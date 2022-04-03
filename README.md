# Symmetric Molecular Dynamics Engine

[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/whitead/symd)
[![PyPI version](https://badge.fury.io/py/symd.svg)](https://badge.fury.io/py/symd)

This is a simple MD engine written in C that can be constraint to obey given symmetry groups at every frame. It can be compiled for any dimensions. It supports the following features:

* NVE, NVT (BAOAB, CSVR thermostats), NPT
* LJ, harmonic, cosine potentials
* Neighbor lists (use `nlj` instead of `lj` for potential)
* Symmetry constraints and all Bravais Lattices

## Purpose

This is not really intended to be a production-ready simulation engine. Instead, it is a reference implementation of symmetry MD. I hope then the algorithm can be incorporated into other engines, or this engine can be extended, for others to use the algorithm. The Python package (see below) should make implementing the method as easy as possible. With that said, the most important files are `box.c` and `group.c`, which contain the key algorithms: `fold_particles` and the NPT algorithms: `make_box`, and `try_rescale`.

## Python Package

There is a companion python package (`symd`) that can be installed via:

```sh
pip install symd
```

This does not include the MD engine, but instead includes necessary details to set-up symmetry constraints like getting affine matrices and Bravais lattice tensors. See [the docs](https://whitead.github.io/symd)

## Compiling


```sh
mkdir build && cd build
cmake ..
make && [sudo] make install
```
The executable is `symdX` where `X` is the number of spatial dimensions

### Dependencies

 * libgsl-dev (GNU Scientific Library headers)
 * CMake


## Python Interface

It is recommended to use the python interface, which simplifies building
the input files and parsing output files. See the examples in the `notebooks` directory.

## Example

An example can be found in the `examples/lj-symm-example` directory.

## Citation
