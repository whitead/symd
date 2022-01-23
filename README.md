Simple Molecular Dynamics Engine
=========================

This code, written in `C`, is meant to be a simple, accurate, and
easily modified (by me) MD engine for 2D, 3D, ND LJ fluids and 1D, 2D, and
3D, ND harmonic oscillators. I use it to test MD methods and sometimes
post my explorations to [my
blog](http://crowsandcats.blogspot.com). The goal is to keep the code
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

Compiling is a little different
than normal. The thermostat, forcefield and integrator are specified
in the `Makefile.` Thus, to compile with a Verlet Integrator, an
Anderson thermostat, and a harmonic force field you use

    cd src
    make harmonic_vverlet_anderson

If you choose `make all`, then all possible versions will be compiled.
Type `make all-single` to use the single threaded versions.

Python Interface
-----------------
TODO: write this section

Try running `python Tester.py` in the python directory

Example
-------------------------
Examples may be run with:

    cd example_example
    ../name_of_executable run_params.txt


TODO/Thoughts/Notes
-------------------------
3. Test 1D lj and 4D lj
4. Write about branching feature
5. checkensemble
6. The cells (used in LJ) do not support open boxes. Not sure whether or not I should implement that. 
7. The g(r) script only supports 3D. I should rewrite it to support 2D and maybe 1D. 
8. If I want to really pursue open boundaries (which I really don't), then I need to remove center-of-mass angular momentum and center-of-mass translation
