#!/usr/bin/python

import sys, copy

from math import *

def print_usage():
    print "Usage: build_lattice.py [dims] [box_1_size] [box_..._size] [box_dims_size] [number]"

if(len(sys.argv) < 2):
    print_usage()
    exit()

n_dims = int(sys.argv[1])

if(len(sys.argv) < 2 + n_dims + 1):
    print_usage()
    exit()


box_size = [float(sys.argv[x]) for x in range(2,2+n_dims)]
number = int(sys.argv[2 + n_dims])
volume = 1.
for x in box_size:
    volume = volume * x
increment = (volume / number) ** (1. / n_dims)


def prepend_emit(array, element):
    array_copy = copy.copy(array)
    array_copy.insert(0, element)
    return array_copy

def enumerate_grid(fxn, dim, sizes, indices):
    if(dim > 0):
        for i in range(sizes[dim]):
            enumerate_grid(fxn, 
                           dim - 1, 
                           sizes, 
                           prepend_emit(indices, i))

    else:
        for i in range(sizes[dim]):
            fxn(prepend_emit(indices, i))

def print_point(point):
    print point,

def print_grid(indices):
    global number
    
    #check if to make sure we still need to make particles
    if(number > 0):
        number -= 1
        global increment
        [print_point(x * increment) for x in indices]
        print ""
        
enumerate_grid(print_grid, n_dims - 1, 
               [int(ceil(x / increment)) for x in box_size], 
               [])
