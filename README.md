Simple Molecular Dynamics Engine
=========================

This code, written in `C`, is meant to be a simple, accurate, and
easily modified (by me) MD engine for 2D, 3D, ND LJ fluids and 1D, 2D, and
3D, ND harmonic oscillators. The goal is to keep the code
base around 1000 lines at all times. 300 lines are for I/O.


The code has the following features:

* Conserves energy
* Linear scaling wrt atom number with cell/verlet lists
* Multithreading with a few well-placed omp pragmas
* Arbitrary dimension number
* Velocity-Verlet integrator
* Anderson thermostat
* Bussi thermostat
* Lennard-Jones (truncated & shifted), Harmonic Forces, and a cosine soft potential for removing overlap
* Input velocities, positions (XYZ Format)
* Input file format is JSON
* Output forces, temperature, energy, velocities, positions (XYZ Format), etc.
* Periodic boundary conditions (or not for harmonic oscillator)
* g(r) script that corrects for [cube-sphere intersection volume](http://crowsandcats.blogspot.com/2013/05/extending-radial-distributions.html)
* script to generate initial structures


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
The executable is `simple-md`

Python Interface
-----------------
TODO: write this section

Try running `python Tester.py` in the python directory

Example
-------------------------
Examples may be run with:

    cd example_example
    simple-md run_params_t.json
