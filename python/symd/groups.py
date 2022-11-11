from dataclasses import dataclass
import numpy as np
import json
import re
import os
from typing import *
from dataclasses import dataclass, asdict

#: Name of space groups
group_names = [
    None,
    [None, "p1", "p1m"],
    [
        None,
        "p1",
        "p2",
        "pm",
        "pg",
        "cm",
        "pmm",
        "pmg",
        "pgg",
        "cmm",
        "p4",
        "p4m",
        "p4g",
        "p3",
        "p3m1",
        "p31m",
        "p6",
        "p6m",
    ],
    [
        None,
        "P1",
        "P-1",
        "P2",
        "P21",
        "C2",
        "Pm",
        "Pc",
        "Cm",
        "Cc",
        "P2/m",
        "P21/m",
        "C2/m",
        "P2/c",
        "P21/c",
        "C2/c",
        "P222",
        "P2221",
        "P21212",
        "P212121",
        "C2221",
        "C222",
        "F222",
        "I222",
        "I212121",
        "Pmm2",
        "Pmc21",
        "Pcc2",
        "Pma2",
        "Pca21",
        "Pnc2",
        "Pmn21",
        "Pba2",
        "Pna21",
        "Pnn2",
        "Cmm2",
        "Cmc21",
        "Ccc2",
        "Amm2",
        "Aem2",
        "Ama2",
        "Aea2",
        "Fmm2",
        "Fdd2",
        "Imm2",
        "Iba2",
        "Ima2",
        "Pmmm",
        "Pnnn",
        "Pccm",
        "Pban",
        "Pmma",
        "Pnna",
        "Pmna",
        "Pcca",
        "Pbam",
        "Pccn",
        "Pbcm",
        "Pnnm",
        "Pmmn",
        "Pbcn",
        "Pbca",
        "Pnma",
        "Cmcm",
        "Cmce",
        "Cmmm",
        "Cccm",
        "Cmme",
        "Ccce",
        "Fmmm",
        "Fddd",
        "Immm",
        "Ibam",
        "Ibca",
        "Imma",
        "P4",
        "P41",
        "P42",
        "P43",
        "I4",
        "I41",
        "P-4",
        "I-4",
        "P4/m",
        "P42/m",
        "P4/n",
        "P42/n",
        "I4/m",
        "I41/a",
        "P422",
        "P4212",
        "P4122",
        "P41212",
        "P4222",
        "P42212",
        "P4322",
        "P43212",
        "I422",
        "I4122",
        "P4mm",
        "P4bm",
        "P42cm",
        "P42nm",
        "P4cc",
        "P4nc",
        "P42mc",
        "P42bc",
        "I4mm",
        "I4cm",
        "I41md",
        "I41cd",
        "P-42m",
        "P-42c",
        "P-421m",
        "P-421c",
        "P-4m2",
        "P-4c2",
        "P-4b2",
        "P-4n2",
        "I-4m2",
        "I-4c2",
        "I-42m",
        "I-42d",
        "P4/mmm",
        "P4/mcc",
        "P4/nbm",
        "P4/nnc",
        "P4/mbm",
        "P4/mnc",
        "P4/nmm",
        "P4/ncc",
        "P42/mmc",
        "P42/mcm",
        "P42/nbc",
        "P42/nnm",
        "P42/mbc",
        "P42/mnm",
        "P42/nmc",
        "P42/ncm",
        "I4/mmm",
        "I4/mcm",
        "I41/amd",
        "I41/acd",
        "P3",
        "P31",
        "P32",
        "R3",
        "P-3",
        "R-3",
        "P312",
        "P321",
        "P3112",
        "P3121",
        "P3212",
        "P3221",
        "R32",
        "P3m1",
        "P31m",
        "P3c1",
        "P31c",
        "R3m",
        "R3c",
        "P-31m",
        "P-31c",
        "P-3m1",
        "P-3c1",
        "R-3m",
        "R-3c",
        "P6",
        "P61",
        "P65",
        "P62",
        "P64",
        "P63",
        "P-6",
        "P6/m",
        "P63/m",
        "P622",
        "P6122",
        "P6522",
        "P6222",
        "P6422",
        "P6322",
        "P6mm",
        "P6cc",
        "P63cm",
        "P63mc",
        "P-6m2",
        "P-6c2",
        "P-62m",
        "P-62c",
        "P6/mmm",
        "P6/mcc",
        "P63/mcm",
        "P63/mmc",
        "P23",
        "F23",
        "I23",
        "P213",
        "I213",
        "Pm-3",
        "Pn-3",
        "Fm-3",
        "Fd-3",
        "Im-3",
        "Pa-3",
        "Ia-3",
        "P432",
        "P4232",
        "F432",
        "F4132",
        "I432",
        "P4332",
        "P4132",
        "I4132",
        "P-43m",
        "F-43m",
        "I-43m",
        "P-43n",
        "F-43c",
        "I-43d",
        "Pm-3m",
        "Pn-3n",
        "Pm-3n",
        "Pn-3m",
        "Fm-3m",
        "Fm-3c",
        "Fd-3m",
        "Fd-3c",
        "Im-3m",
        "Ia-3d",
    ],
]


@dataclass
class Group:
    """A group with information on genpos, specpos, Bravais lattice"""

    #: Bravais lattice name
    lattice: str
    #: general positions as affine matrix strings (call :func:`str2mat` to convert to matrix)
    genpos: List[str]
    #: list of :obj:`Group` special positions
    specpos: List[object]
    #: name of the group
    name: Optional[str] = None
    #: string indicating asymmetric unit (call :func:`asymm_constraints` to convert to lambda)
    asymm_unit: Optional[str] = None


def _dict2group(d, name=None):
    specpos = []
    if "specpos" in d:
        for s in d["specpos"]:
            g = Group(
                lattice=d["lattice"], genpos=s["sites"], specpos=[], name=s["name"]
            )
            specpos.append(g)
    return Group(
        d["lattice"],
        d["genpos"],
        specpos,
        name,
        d["asymm_unit"] if "asymm_unit" in d else None,
    )


def str2mat(s: str) -> np.ndarray:
    """
    Convert affine matrix specified in xyz notation to matrix. For example, -x, y - x, z.
    Can be 2D or 3D.

    :param s: string in xyz notation
    :return: np.ndarray: matrix
    """
    rows = []
    N = len(s.split(","))
    env = {"x": np.array([1, 0, 0]), "y": np.array([0, 1, 0]), "z": np.array([0, 0, 1])}
    fake_env = {"x": 0, "y": 0, "z": 0}
    for i, si in enumerate(s.split(",")):
        # treat implicit multiplication - 2x = 2 * x
        si = re.sub("(?<=\d)(?=x) | (?<=\d)(?=y) | (?<=\d)(?=z)", "*", si, flags=re.X)
        r = [0] * N
        l = {}
        # use fake ones to get translation
        exec("translation = " + si.strip(), fake_env, l)
        exec("scale = " + si.strip(), env, l)
        # expand
        if type(l["scale"]) != np.ndarray:
            t = np.zeros(N)
            t[i] = l["translation"]
            l["translation"] = t
            l["scale"] = t
        # remove trans and add
        rows.append(
            np.append((l["scale"] - l["translation"])[:N], np.sum(l["translation"]))
        )
    rows.append(np.array(N * [0] + [1]))
    result = np.vstack(rows)
    return result


def asymm_constraints(
    s: str,
) -> Union[Callable[[float, float, float], bool], Callable[[float, float], bool]]:
    """
    Converts inequalities in xyz notation to lambda function. For example,
    0\u2264x\u22642/3;0\u2264y\u22641/3;x\u2264(1+y)/2;y\u2264x/2 is converted to a
    lambda function that returns True if the given point satisfies the constraints. Works in 2D and 3D.

    :param s: string in xyz notation
    :return: lambda function that is ``True`` when the point is in the asymmetric unit
    """
    s = s.replace("≤", "<=")
    env = {}
    in3d = "z" in s
    exec("from math import *", env)
    funcs = []
    for i, si in enumerate(s.split(";")):
        # treat implicit multiplication - 2x = 2 * x
        si = re.sub("(?<=\d)(?=x) | (?<=\d)(?=y) | (?<=\d)(?=z)", "*", si, flags=re.X)
        l = {}
        if in3d:
            exec(f"l{i} = lambda x,y,z:" + si, env, l)
        else:
            exec(f"l{i} = lambda x,y:" + si, env, l)
        funcs.append(l[f"l{i}"])
    if in3d:
        return lambda x, y, z: sum([f(x, y, z) for f in funcs]) == len(funcs)
    else:
        return lambda x, y: sum([f(x, y) for f in funcs]) == len(funcs)


projectors2d = {
    "Square": np.array([4 * [1], 4 * [0], 4 * [0], 4 * [1]]),
    "Rectangular": np.array([[1, 1, 0, 0], 4 * [0], 4 * [0], [0, 0, 1, 1]]),
    "Hexagonal": np.array([4 * [1], 4 * [-1 / 2], 4 * [0], 4 * [np.sqrt(3) / 2]]),
    "Oblique": np.eye(4),
}

projectors3d = {
    "Hexagonal": np.array(
        [
            3 * [1] + 6 * [0],  # ax
            3 * [-1 / 2] + 6 * [0],  # bx
            9 * [0],  # cx
            9 * [0],  # ay
            3 * [np.sqrt(3) / 2] + 6 * [0],  # by
            9 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            6 * [0] + 3 * [1],  # cz
        ]
    ),
    "Cubic": np.array(
        [
            3 * [1] + 6 * [0],  # ax
            9 * [0],  # bx
            9 * [0],  # cx
            9 * [0],  # ay
            3 * [1] + 6 * [0],  # by
            9 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            3 * [1] + 6 * [0],  # cz
        ]
    ),
    "Tetragonal": np.array(
        [
            6 * [1] + 3 * [0],  # ax
            9 * [0],  # bx
            9 * [0],  # cx
            9 * [0],  # ay
            6 * [1] + 3 * [0],  # by
            9 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            6 * [0] + 3 * [1],  # cz
        ]
    ),
    "Triclinic": np.eye(9),
    "Monoclinic": np.array(  # TODO: might be missing potential rotation around z
        [
            [1] + 8 * [0],  # ax
            9 * [0],  # bx
            2 * [0] + [1] + 6 * [0],  # cx
            9 * [0],  # ay
            6 * [1] + 3 * [0],  # by
            5 * [0] + [1] + 3 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            8 * [0] + [1],  # cz
        ]
    ),
    "Orthorhombic": np.array(  # TODO: might be missing potential rotation around z
        [
            3 * [1] + 6 * [0],  # ax
            9 * [0],  # bx
            9 * [0],  # cx
            9 * [0],  # ay
            3 * [0] + 3 * [1] + 3 * [0],  # by
            9 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            6 * [0] + 3 * [1],  # cz
        ]
    ),
    "Trigonal": np.array(  # cubic, plus use b_y as shear
        [
            4 * [1] + 6 * [0],  # ax
            3 * [0] + [1] + 5 * [0],  # bx
            3 * [0] + [1] + 5 * [0],  # cx
            3 * [0] + [1] + 5 * [0],  # ay
            4 * [1] + 6 * [0],  # by
            3 * [0] + [1] + 5 * [0],  # cy
            3 * [0] + [1] + 5 * [0],  # az
            3 * [0] + [1] + 5 * [0],  # bz,
            4 * [1] + 6 * [0],  # cz
        ]
    ),
}


def _projector_key(p):
    # sort so asymmetric unit is valid
    dim = p.shape[0]
    # see which is closest to identity
    s = np.sum((p - np.eye(dim)) ** 2)
    # prefer positive coordinates (in x first)
    # s += np.sum((p < 0) @ np.arange(dim, 0, -1).T)
    return s


def _write_group(f, name, group, dim):
    # tiling code needs to be updated for more than 2D
    def fmt(n):
        if n is None:
            return list(np.round(np.eye(dim**2).astype(float).flatten(), 8))
        return list(np.round(n.astype(float).flatten(), 8))

    try:
        projector = (
            projectors2d[group.lattice] if dim == 2 else projectors3d[group.lattice]
        )
    except KeyError:
        projector = None
    members = [str2mat(s) for s in group.genpos]
    result = {
        "name": name,
        "size": len(members),
        "members": [],
        "projector": fmt(projector),
    }
    # if we do not have identity as zeroth, make one with
    # least translation be zero
    if not np.allclose(members[0][:-1, :-1], np.eye(dim)):
        members.sort(key=_projector_key)
    for m in members:
        result["members"].append(fmt(m))
    # should be same for all, so use last
    dof = np.sum(np.sum(m[:-1, :-1] ** 2, axis=1) > 0)
    result["dof"] = int(dof)
    json.dump(result, f, indent=True)


def load_group(gnum: int, dim: int) -> Group:
    """
    Load one of the 2D planar groups or 3D space groups that tile space. The :obj:`Group`
    contains the name of the Bravais lattice, the general positions,
    and a list of special positions.

    :param gnum: group number (Hall number)
    :param dim: dimensionality of space
    :return: The :obj:`Group`
    """
    gnum = str(gnum)
    from importlib_resources import files
    import symd.data

    fp = files(symd.data).joinpath(f"{dim}dgroups.json")
    with open(fp, "r") as f:
        all_groups = json.load(f)
    if gnum not in all_groups:
        raise KeyError("Could not find group " + gnum)
    group = all_groups[gnum]
    return _dict2group(group)


def prepare_input(
    group: Union[Group, int], dim: int, N: int, name: str, dir="."
) -> List[str]:
    """
    Prepare input files for running symmetry MD.

    :param group: group number (Hall number) or :obj:`Group`
    :param dim: dimensionality of space
    :param N: number of atoms
    :param name: name of the group (used to name output files)
    :param dir: directory to write files to
    :return: list of input files
    """
    if type(group) is int:
        group = load_group(group, dim)
    asymm_unit = asymm_constraints(group.asymm_unit)
    with open(os.path.join(dir, f"{name}.json"), "w") as f:
        _write_group(f, name, group, dim)
    paths = []
    for i, g in enumerate(group.specpos):
        fn = os.path.join(dir, f"{name}-{i:02d}.json")
        paths.append(fn)
        with open(fn, "w") as f:
            _write_group(f, name + f"-{i}", g, dim)
    Ni = N
    with open(os.path.join(dir, f"{name}.dat"), "w") as f:
        while Ni > 0:
            x = np.random.uniform()
            y = np.random.uniform()
            z = np.random.uniform()
            if (dim == 2 and asymm_unit(x, y)) or (dim == 3 and asymm_unit(x, y, z)):
                f.write(f"{x} {y}\n" if dim == 2 else f"{x} {y} {z}\n")
                Ni -= 1
    return paths


def cell_nparticles(group: Group, genpos: int, *specpos: int) -> int:
    """Get number of unit cell particles given genpos and specpos

    :param group: The :obj:`Group`
    :param genpos: The number of particles in the general positions in asymmetric unit
    :param specpos: The numbers of particles in the special positions
    :return: The number of particles in the unit cell
    """
    N = 0

    if specpos is None:
        specpos = []
    for i, w in enumerate(specpos):
        genpos -= w
        if i == len(group.specpos):
            raise ValueError("Too many specpos")
        N += w * len(group.specpos[i].genpos)

    return N + genpos * len(group.genpos)


def _sign(x):
    return bool(x > 0) - bool(x < 0)


def _levi_civta(index):
    p = 1
    d = len(index)
    for i in range(d):
        for j in range(i + 1, d):
            p *= _sign(index[j] - index[i])
    return p


def _rvolume(b, v, i, index):
    d = len(index)
    if i == d:
        return _levi_civta(index) * v
    vi = 0
    for j in range(d):
        index[i] = j
        vi += _rvolume(b, b[j][i] * v, i + 1, index)
    return vi


def cell_volume(b: np.ndarray) -> float:
    """
    Compute volume of given unit cell in arbitrary dimension

    :param b: lattice vectors as columns
    :return: volume of unit cell
    """
    index = [0] * len(b)
    return _rvolume(b, 1, 0, index)


def project_cell(cell: np.ndarray, projector: Union[str, np.ndarray]) -> np.ndarray:
    """
    Project unit cell to constraints of Bravais lattice.

    :param cell: unit cell as columns
    :param projector: projector tensor or str of Bravais lattice
    """
    ndim = cell.shape[0]
    if type(projector) is str:
        projector = projectors2d[projector] if ndim == 2 else projectors3d[projector]
    fub = cell.flatten()
    fb = np.array(projector).reshape(ndim**2, ndim**2) @ fub
    return fb.reshape(ndim, ndim)


def get_cell(
    number_density: float,
    group: Union[int, Group],
    dim: int,
    n: int,
    w: Optional[List[int]] = None,
) -> List[float]:
    """
    Compute unit cell given number density, group, and number of particles in asymmetric unit.

    :param number_density: number density of particles
    :param group: group number (Hall number) or :obj:`Group`
    :param dim: dimensionality of space
    :param n: number of particles in asymmetric unit general positions
    :param w: list of number of particles in special positions
    :return: flattened unit cell (for use in symd MD engine)
    """
    import scipy.optimize as opt

    if w is None:
        w = []
    if type(group) is int:
        group = load_group(group, dim)
    pname = group.lattice
    projector = projectors2d[pname] if dim == 2 else projectors3d[pname]
    N = cell_nparticles(group, n, *w)
    cell = np.eye(dim) * N

    def obj(s):
        c = project_cell(cell * s, projector)
        v = cell_volume(c)
        return (N / v - number_density) ** 2

    result = opt.minimize(obj, x0=1, bounds=[(1e-5, 1e5)])
    return list((result.x * cell).flatten().astype(float))
