Symmetric Molecular Dynamics Engine
=========================

TODO:

Improve speed for dimensions (low priority)
1. Switch to dynamic arrays (as you go) DONE
2. n_dims becomes pre-processor

Nlist broken?
1. Investigate - might just be inconsistent when we had jumps but is ok

Groups
1. Better checks on validity!

NPT
1. Need to work with Bravais lattice directly (encode in json), otherwise we violate


PBC
1. Make flag or something so that we can actually have both rect/cubic PBC and our translation groups

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
