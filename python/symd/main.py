from operator import concat
import symd
import numpy as np
import ternviz
import click


def rmsd(p1, p2):
    return np.mean((p1 - p2) ** 2, axis=(1, 2))


def write_xyz(filename, name, traj):
    if len(traj.shape) == 2:
        traj = traj[None]
    with open(filename, "w") as f:
        T, N, _ = traj.shape
        for i in range(T):
            f.write(f"{N}\n")
            f.write(f"{name} t={i}\n")
            for j in range(N):
                f.write("C ")
                for k in range(3):
                    f.write(f"{traj[i,j,k]} ")
                f.write("\n")


def compute_crystal(
    n,
    group,
    w=None,
    ndims=2,
    retries=5,
    steps=10**6,
    steps2=5 * 10**3,
    starting_density=0.2,
):
    # trying to have n be number in UNIT CELL
    # so have to adjust for group size
    m = len(symd.groups.load_group(group, ndims).genpos)
    n = max(2, n // m)
    gname = symd.group_names[ndims][group]
    if w is not None:
        n += sum(w)
        name = f"{gname}-{n}-{sum(w)}"
    else:
        name = f"{gname}-{n}"
    print("Simulating", n, "particles:", name)
    # break out the try/except because we will accept failed NPT (because it jams so hard)
    for i in range(retries):
        np.random.seed(i)
        cell = symd.groups.get_cell(starting_density, group, ndims, n, w)
        # NPT
        md = symd.Symd(
            nparticles=n,
            cell=cell,
            ndims=ndims,
            images=2,
            force="lj",
            wyckoffs=w,
            group=group,
            steps=steps,
            exeDir=f"crystal-{name}",
            pressure=0.25,
            temperature=0.1,
            start_temperature=0.5,
        )
        try:
            md.remove_overlap()
        except RuntimeError as e:
            continue
        md.log_positions(frames=360)
        try:
            md.run()
        except RuntimeError as e:
            d = md.number_density()
            if d < 0.5:
                print("Not dense enough, retrying", d)
                continue

        # NVT
        md.runParams["start_temperature"] = 0.05
        md.runParams["temperature"] = 1e-4
        md.runParams["box_update_period"] = 0
        md.runParams["langevin_gamma"] = 0.5
        md.runParams["steps"] = steps // 4
        md.log_positions(filename="equil.xyz")
        try:
            md.run()
        except RuntimeError as e:
            continue
        config = md.positions[-1]

        # Stability
        fp = np.loadtxt(md.runParams["final_positions"])
        # changing group, so need to read projected cell
        cell = md.read_cell(bravais=True)
        m = fp.shape[0]
        md2 = symd.Symd(
            nparticles=m,
            cell=cell,
            ndims=ndims,
            images=2,
            force="lj",
            wyckoffs=None,
            group=1,
            steps=steps2,
            exeDir=f"melt-{name}",
            temperature=None,
            start_temperature=0.0,
        )
        # run once to get melting traj
        # then again for longer with longer period
        md2.log_positions(period=10)
        md2.runParams["start_positions"] = md.runParams["final_positions"]
        try:
            md2.run()
        except RuntimeError as e:
            continue
        # traj = md2.positions
        traj = md.positions
        csm = rmsd(md2.positions[:, :m], md2.positions[0, :m])
        # csm = []
        # for i in range(md2.positions.shape[0]):
        #    csm.append(compute_symm(md2.positions[i], group, md2.read_cell(), ndims, n))
        return (
            name,
            config,
            md2.positions[-1],
            md2.number_density(),
            csm,
            traj,
            md.number_density(),
            np.arange(0, steps2, 10) * md2.runParams["time_step"],
        )
    return None


@click.command()
@click.argument("n", type=int)
@click.argument("group", type=int)
@click.option("--w", default=None, type=int)
@click.option("--ndims", default=3)
@click.option("--ffmpeg", default="ffmpeg")
@click.option("--vmd", default="vmd")
def crystal(n, group, w, ndims, ffmpeg, vmd):
    if w is not None:
        w = [1] * w
    result = compute_crystal(n, group, w, ndims)
    if result:
        name, config, config2, nd, csm, traj, rho, time = result
        fc = f"{name}-crystal.xyz"
        ft = f"{name}-traj.xyz"
        disp_name = f"{name} (p={rho:.2f})"
        write_xyz(fc, name, config)
        write_xyz(ft, name, traj)
        print("Finished crystallization, now rendering")
        tid = name + "-t"
        cid = name + "-c"
        ternviz.render(
            ft, 800, id=tid, script_name="render-points.vmd", color="white", vmd=vmd
        )
        mt = ternviz.movie(tid, ffmpeg=ffmpeg, short_name=disp_name, color="black")
        ternviz.render(
            fc, 800, id=cid, script_name="render-points.vmd", color="white", vmd=vmd
        )
        mc = ternviz.movie(cid, ffmpeg=ffmpeg, short_name=disp_name, color="black")
        print(ternviz.concat([mt, mc], name, ffmpeg=ffmpeg))
    else:
        print("Failed to compute crystal")
