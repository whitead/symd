Simple Molecular Dynamics Engine
=========================

This code, written in `C`, provides a very simple minimal molecular
dynamics engine. I wouldn't recommend using it. I use it to test MD
methods and sometimes post my explorations to [my
blog](http://crowsandcats.blogspot.com).

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
To run the example, try the following commands after compiling:
    cd example
    ../harmonic_verlet_anderson run_params.txt

TODO/Thoughts
-------------------------
1. Implement the Bussi-Parrinello thermostat
2. I hate it, but I should go ahead and implement a non `n^2` force calculation
3. Perhaps I should try implement that Hilbert curve dynamic R-Tree decomposition method for fun
4. Try those prior belief molecular dynamics equations I wrote down on that paper I lost in my office
5. Use actual xyz format
6. I'm out of chunky peanut butter, I should buy more