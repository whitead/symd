Simple Molecular Dynamics Engine
=========================

This code, written in `C`, provides a very simple minimal molecular
dynamics engine. I wouldn't recommend using it. I use it to test MD
methods and sometimes post my explorations to [my
blog](http://crowsandcats.blogspot.com).


The code has the following features:

* Conserves energy
* Conserves conjugate energy thing in Bussi thermostat
* Arbitrary dimension number (1-3 I would recommend)
* Velocity-Verlet integrator
* Anderson thermostat
* Bussi thermostat
* LJ or Harmonic Forces
* Input velocites, positions (XYZ Format)
* Output forces, temperature, energy, velocities, positions (XYZ Format)
* Periodic boundary conditions


Compiling
-------------------------

If you really want to use this code, compiling is a little different
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

TODO/Thoughts
-------------------------
1. Remove COM motion
2. I hate it, but I should go ahead and implement a non `n^2` force calculation at some point
