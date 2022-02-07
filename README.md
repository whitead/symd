Symmetric Molecular Dynamics Engine
=========================

TODO:

Improve speed for dimensions (low priority)
1. Switch to dynamic arrays (as you go) DONE
2. n_dims becomes pre-processor
3. Use OPenMP in group action

Nlist broken?
1. Investigate - might just be inconsistent when we had jumps but is ok

Groups
1. Better checks on validity!

NPT
1. Finish Python C code generation to work with Bravais Lattices
2.


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
