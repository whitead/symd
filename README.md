Simple Molecular Dynamics Engine
=========================

This code, written in `C`, provides a very simple minimal molecular
dynamics engine. I wouldn't recommend using it. I use it to test MD
methods and sometimes post my explorations to [my
blog](http://crowsandcats.blogspot.com).


The code has the following features:

* Conserves energy
* Arbitrary dimension number (1-5 I would recommend)
* Verlet integrator
* Anderson thermostat
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
    make harmonic_verlet_anderson

If you choose `make all`, then all possible versions will be compiled.

Example
-------------------------
Examples may be run with:

    cd example_example
    ../name_of_executable run_params.txt

TODO/Thoughts
-------------------------
1. Alright, tried to implement Bussi-Parrinello with GSL. Haven't compiled yet. It's bed time.
1. Remove COM motion
2. Implement the Bussi-Parrinello thermostat. The code and implementation details may be found in resamplekin.f90. The only detail to add is that the random number distributions should be discarded and  the GNU scilibrary should be used. The center-of-mass motion needs to removed as well for the whole system and its corresponding degrees of freedom removed for thermostat calculations. Also, note that some of the implemenations floating around for the Bussi thermostat draw N random numbers, where N is the number of particles. Don't do that.
2. I hate it, but I should go ahead and implement a non `n^2` force calculation at some point


