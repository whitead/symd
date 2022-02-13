Symmetric Molecular Dynamics Engine
=========================

TODO:

Improve speed for dimensions (low priority)
1. Switch to dynamic arrays (as you go) DONE
2. n_dims becomes pre-processor
3. Use OPenMP in group action

Nlist broken?
1. Investigate - might just be inconsistent when we had jumps but is ok

1. unscale_wrap_coords, unscal_coords, tile_coords,unscale_coords
2. only load asymmetric unit
3. tiling should not be done via groups


volume - https://en.wikipedia.org/wiki/Levi-Civita_symbol#Generalization_to_n_dimensions
https://math.stackexchange.com/questions/1605690/how-to-extend-the-parallelepiped-volume-formula-to-higher-dimensions


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
