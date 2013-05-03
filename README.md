Simple Molecular Dynamics Engine
=========================

This code, written in `C`, is meant to be a simple, accurate, and
easily modified (by me) MD engine for 2D, 3D LJ fluids and 1D, 2D,
and 3D harmonic oscillators. I use it to test MD methods and sometimes
post my explorations to [my
blog](http://crowsandcats.blogspot.com). The goal is to keep the code
base around 1500 lines at all times.


The code has the following features:

* Conserves energy
* Linear scaling wrt atom number with cell/verlet lists
* Arbitrary dimension number (implementing cell/verlet lists in arbitrary dimensions was fun)
* Velocity-Verlet integrator
* Anderson thermostat
* Bussi thermostat (conserves conjugate energy)
* Lennary-Jones (truncated & shifted) or Harmonic Forces
* Input velocities, positions (XYZ Format)
* Output forces, temperature, energy, velocities, positions (XYZ Format), etc.
* Periodic boundary conditions (or not)
* g(r) script that corrects for cube-sphere intersection volume
* script to generate initial structures
* Removes COM motion


Compiling
-------------------------
Compiling is a little different
than normal. The thermostat, forcefield and integrator are specified
in the `Makefile.` Thus, to compile with a Verlet Integrator, an
Anderson thermostat, and a harmonic force field you use

    cd src
    make harmonic_vverlet_anderson

If you choose `make all`, then all possible versions will be compiled.

Example
-------------------------
Examples may be run with:

    cd example_example
    ../name_of_executable run_params.txt

TODO/Thoughts/Notes
-------------------------
1. Make a few plots showing conservation, linear scaling
2. Write down g(r) cube-sphere intersection information
3. Test 1D lj and 4D lj
4. Write about branching feature
5. Get checkensemble to work
6. The cells (used in LJ) do not support open boxes. Not sure whether or not I should implement that. 
7. The g(r) script only supports 3D. I should rewrite it to support 2D and maybe 1D. 3D cube-sphere intersection sounds like I would die.