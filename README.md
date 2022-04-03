Symmetric Molecular Dynamics Engine
=========================

This is a simple MD engine written in C that can be constraint to obey given symmetry groups at every frame. It can be compiled for any dimensions. It supports the following features:

* NVE, NVT (BAOAB, CSVR thermostats), NPT
* LJ, harmonic, cosine potentials
* Neighbor lists (use `nlj` instead of `lj` for potential)
* Symmetry constraints and all Bravais Lattices

Compiling
-------------------------

```sh
mkdir build && cd build
cmake ..
make && [sudo] make install
```
The executable is `symdX` where `X` is the number of spatial dimensions

### Dependencies

 * libgsl-dev (GNU Scientific Library headers)
 * CMake


Python Interface
-----------------

It is recommended to use the python interface, which simplifies building
the input files and parsing output files. See the examples in the `notebooks` directory.

Example
-------------------------

An example can be found in the `examples/lj-symm-example` directory.
