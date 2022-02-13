Symmetric Molecular Dynamics Engine
=========================

TODO:

Find TODO

Nlist broken?
1. Revise to use enclosing cubes and check implementation
2. Investigate - might just be inconsistent when we had jumps but is ok


NPT
1. volume - https://en.wikipedia.org/wiki/Levi-Civita_symbol#Generalization_to_n_dimensions, https://math.stackexchange.com/questions/1605690/how-to-extend-the-parallelepiped-volume-formula-to-higher-dimensions

3D
1. Write out matrix inverse
2. Compiler Flag

Science
1. Wyckoff Monte Carlo
2. BAOAB
3. 3D


Compiling
-------------------------

### Dependencies

 * libgsl-dev (GNU Scientific Library headers)
 * CMake

```sh
mkdir build && cd build
cmake ..
make && make install
```
The executable is `symd`

Python Interface
-----------------
TODO: write this section

Try running `python Tester.py` in the python directory

Example
-------------------------
Examples may be run with:

    cd example_example
    symd run_params_t.json
