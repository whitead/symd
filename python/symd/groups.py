import numpy as np
import json
import re
import os


def str2mat(s):
    rows = []
    N = len(s.split(','))
    env = {'x': np.array([1, 0, 0]), 'y': np.array(
        [0, 1, 0]), 'z': np.array([0, 0, 1])}
    fake_env = {'x': 0, 'y': 0, 'z': 0}
    for i, si in enumerate(s.split(',')):
        # treat implicit multiplication - 2x = 2 * x
        si = re.sub('(?<=\d)(?=x) | (?<=\d)(?=y) | (?<=\d)(?=z)',
                    '*', si, flags=re.X)
        r = [0] * N
        l = {}
        # use fake ones to get translation
        exec('translation = ' + si.strip(), fake_env, l)
        exec('scale = ' + si.strip(), env, l)
        # expand
        if type(l['scale']) != np.ndarray:
            t = np.zeros(N)
            t[i] = l['translation']
            l['translation'] = t
            l['scale'] = t
        # remove trans and add
        rows.append(
            np.append((l['scale'] - l['translation'])[:N], np.sum(l['translation'])))
    rows.append(np.array(N * [0] + [1]))
    result = np.vstack(rows)
    return result


def asymm_constraints(s):
    s = s.replace('≤', '<=')
    env = {}
    in3d = 'z' in s
    exec('from math import *', env)
    funcs = []
    for i, si in enumerate(s.split(';')):
        # treat implicit multiplication - 2x = 2 * x
        si = re.sub('(?<=\d)(?=x) | (?<=\d)(?=y) | (?<=\d)(?=z)',
                    '*', si, flags=re.X)
        l = {}
        if in3d:
            exec(f'l{i} = lambda x,y,z:' + si, env, l)
        else:
            exec(f'l{i} = lambda x,y:' + si, env, l)
        funcs.append(l[f'l{i}'])
    if in3d:
        return lambda x, y, z: sum([f(x, y, z) for f in funcs]) == len(funcs)
    else:
        return lambda x, y: sum([f(x, y) for f in funcs]) == len(funcs)


projectors2d = {
    'Square':
    np.array([
        4 * [1],
        4 * [0],
        4 * [0],
        4 * [1]
    ]),
    'Rectangular':
    np.array([
        [1, 1, 0, 0],
        4 * [0],
        4 * [0],
        [0, 0, 1, 1]
    ]),
    'Hexagonal':
    np.array([
        4 * [1],
        4 * [-1/2],
        4 * [0],
        4 * [np.sqrt(3)/2]
    ]),
    'Oblique': np.eye(4),
}

projectors3d = {
    'Hexagonal':
    np.array([
        3 * [1] + 6 * [0],  # ax
        3 * [-1/2] + 6 * [0],  # bx
        9 * [0],  # cx
        9 * [0],  # ay
        3 * [0] + 3 * [np.sqrt(3)/2] + 3 * [0],  # by
        9 * [0],  # cy
        9 * [0],  # az
        9 * [0],  # bz,
        6 * [0] + 3 * [1],  # cz
    ]),
    'Cubic':
        np.array([
            3 * [1] + 6 * [0],  # ax
            9 * [0],  # bx
            9 * [0],  # cx
            9 * [0],  # ay
            3 * [0] + 3 * [1] + 3 * [0],  # by
            9 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            6 * [0] + 3 * [1],  # cz
        ]),
    'Tetragonal':
    np.array([
        6 * [1] + 3 * [0],  # ax
        9 * [0],  # bx
        9 * [0],  # cx
        9 * [0],  # ay
        6 * [1] + 3 * [0],  # by
        9 * [0],  # cy
        9 * [0],  # az
        9 * [0],  # bz,
        6 * [0] + 3 * [1],  # cz
    ]),
    'Triclinic': np.eye(9),
    'Monoclinic':  # TODO: might be missing potential rotation around z
        np.array([
            6 * [1] + 3 * [0],  # ax
            9 * [0],  # bx
            2 * [0] + [1] + 6 * [0],  # cx
            9 * [0],  # ay
            6 * [1] + 3 * [0],  # by
            5 * [0] + [1] + 3 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            8 * [0] + [1],  # cz
        ]),
    'Orthorhombic':  # TODO: might be missing potential rotation around z
        np.array([
            3 * [1] + 6 * [0],  # ax
            9 * [0],  # bx
            9 * [0],  # cx
            9 * [0],  # ay
            3 * [0] + 3 * [1] + 3 * [0],  # by
            9 * [0],  # cy
            9 * [0],  # az
            9 * [0],  # bz,
            6 * [0] + 3 * [1],  # cz
        ]),
}
projectors3d['Trigonal'] = projectors3d['Hexagonal']


def _projector_key(p):
    # sort so asymmetric unit is valid
    dim = p.shape[0]
    # see which is closest to identity
    s = np.sum((p - np.eye(dim))**2)
    # prefer positive coordinates (in x first)
    #s += np.sum((p < 0) @ np.arange(dim, 0, -1).T)
    return s


def write_group(f, name, group, dim):
    # tiling code needs to be updated for more than 2D
    def fmt(n):
        if n is None:
            return list(np.round(np.eye(dim**2).astype(float).flatten(), 8))
        return list(np.round(n.astype(float).flatten(), 8))
    try:
        projector = projectors2d[group['lattice']
                                 ] if dim == 2 else projectors3d[group['lattice']]
    except KeyError:
        projector = None
    key = 'genpos'
    if key not in group:
        key = 'sites'
    members = [str2mat(s) for s in group[key]]
    result = {'name': name, 'size': len(
        members), 'members': [], 'projector': fmt(projector)}
    # if we do not have identity as zeroth, make one with
    # least translation be zero
    if not np.allclose(members[0][:-1, :-1], np.eye(dim)):
        members.sort(key=_projector_key)
    for m in members:
        result['members'].append(fmt(m))
    # should be same for all, so use last
    dof = np.sum(np.sum(m[:-1, :-1]**2, axis=1) > 0)
    result['dof'] = int(dof)
    json.dump(result, f, indent=True)


def load_group(gnum, dim):
    gnum = str(gnum)
    from importlib_resources import files
    import symd.data
    fp = files(symd.data).joinpath(
        f'{dim}dgroups.json')
    with open(fp, 'r') as f:
        all_groups = json.load(f)
    if gnum not in all_groups:
        raise KeyError('Could not find group ' + gnum)
    group = all_groups[gnum]
    return group


def prepare_input(gnum, dim, N, name, dir='.'):
    group = load_group(gnum, dim)
    asymm_unit = asymm_constraints(group['asymm_unit'])
    with open(os.path.join(dir, f'{name}.json'), 'w') as f:
        write_group(f, name, group, dim)
    paths = []
    for i, g in enumerate(group['specpos']):
        fn = os.path.join(dir, f'{name}-{i:02d}.json')
        paths.append(fn)
        with open(fn, 'w') as f:
            write_group(f, name + f'-{i}', g, dim)
    Ni = N
    with open(os.path.join(dir, f'{name}.dat'), 'w') as f:
        while Ni > 0:
            x = np.random.uniform()
            y = np.random.uniform()
            z = np.random.uniform()
            if (dim == 2 and asymm_unit(x, y)) or (dim == 3 and asymm_unit(x, y, z)):
                f.write(f'{x} {y}\n' if dim == 2 else f'{x} {y} {z}\n')
                Ni -= 1
    return paths
