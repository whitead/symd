Symmetric Molecular Dynamics Engine
=========================

TODO:

1. Switch to dynamic arrays
2. n_dims becomes pre-processor

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
