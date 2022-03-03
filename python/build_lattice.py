#!/usr/bin/python

import copy
from math import *


def main():
    def print_usage():
        print(
            "Usage: build_lattice.py [dims] [box_1_size] [box_..._size] [box_dims_size] [number]"
        )

    if len(sys.argv) < 2:
        print_usage()
        exit()

    n_dims = int(sys.argv[1])

    if len(sys.argv) < 2 + n_dims + 1:
        print_usage()
        exit()

    cell_size = [float(sys.argv[x]) for x in range(2, 2 + n_dims)]
    number = int(sys.argv[2 + n_dims])
    increment = increment_size(n_dims, cell_size, number)
    enumerate_grid(
        lambda x, y: print_grid(x, y, increment, number),
        n_dims - 1,
        [int(ceil(x / increment)) for x in cell_size],
        [],
    )


def increment_size(n_dims, cell_size, number):
    volume = 1
    for x in cell_size:
        volume = volume * x
    return (volume) ** (1.0 / n_dims) / number


def prepend_emit(array, element):
    array_copy = copy.copy(array)
    array_copy.insert(0, element)
    return array_copy


def enumerate_grid(fxn, dim, sizes, indices, count=0):
    if dim > 0:
        for i in range(sizes[dim]):
            count += enumerate_grid(fxn, dim - 1, sizes,
                                    prepend_emit(indices, i), count)

    else:
        for i in range(sizes[dim]):
            fxn(prepend_emit(indices, i), count)
            count += 1
    return count


def print_point(point):
    print(f"{point:8.3f}", end=" ")


def print_grid(indices, count, increment, maxcount=-1):

    # check if to make sure we still need to make particles
    if maxcount > 0 and count < maxcount:
        [print_point(x * increment) for x in indices]
        print("")


if __name__ == "__main__":
    import sys

    main()
