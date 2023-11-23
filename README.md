# tides
Software under development for numerically investigating the tidal evolution of celestial bodies based on the physical theory developed by Clodoaldo Ragazzo and Lucas Ruiz.

This work is part of a postdoctoral project from the São Paulo Research Foundation (FAPESP - Grants 2021/11306-0 and 2022/12785-1), carried out by Vitor M. de Oliveira at the Institute of Mathematics and Statistics of the University of São Paulo (Brazil) and at the Department of Physics of the University of Coimbra (Portugal), under the supervision of prof. Clodoaldo Ragazzo (IME/USP) and with the collaboration of prof. Alexandre Correia (CFisUC).

Author: V. M. de Oliveira

Last update on this file: November 21th 2023

## Important notes
The software simulates the tidal evolution of any number of celestial bodies interacting gravitationally with each other. The rheology model adopted here is the generalized Voigt one, which can be reduced to the Maxwell model. The equations of motion are numerically integrated using a Prince-Dormand Runge-Kutta scheme of 8th and 9th order and fixed stepsize from the GNU Scientific Library.

Up to this point, the simulation assumes that there is no prestress and that the angular velocity of each body is parallel to its biggest moment of inertia.

## How to compile
```makefile
make
```

This command will compile the code using the compiler at ``cc``. If one wants to use another compiler, set the environment variable ``CC`` alongside the command.

## How to run
```makefile
make run INPUT=configuration_file.dat
```

The ``configuration_file.dat`` should include the simulation name under ``name``, the location of the input files under ``input_folder``, the location of the output files under ``output_folder``, the name of the files containing the system and the integrator specifications, under ``system_specs`` and ``integrator_specs``, respectively. The number of bodies used can also be given in this file under ``number_of_bodies``. If this value is not passed to the program, the total number of bodies in ``system_specs`` is assumed. The results are written in the ``output_folder`` as ``results_name``.

## Additional features
After running the program, it is possible to calculate the orbital elements of every body using

```makefile
make orbit INPUT=configuration_file.dat
```

And also to calculate the spin variables of every body using

```makefile
make spin INPUT=configuration_file.dat
```

After that, it is possible to plot the time evolution of the main variables in the system via ``Gnuplot`` using

```makefile
make plot INPUT=configuration_file.dat
```

## All
If one wants to run the code and calculate all of the additional features, it is possible to do so by typing

```makefile
make all INPUT=configuration_file.dat
```

## Examples
There are also some example files which can be run using

```makefile
make examples
```

This command runs the program with the input files in the folder ``examples``, which contains 3 example files: the specs for the Earth-Moon system with no obliquity, zero orbital eccentricity, and the Moon being treated as a point mass; the specs for the Earth-Moon-Sun system with all objects being treated as point masses; and the specs for rigid Earth.

## Main references
Ragazzo C., Ruiz L. S. Viscoelastic tides: models for use in Celestial Mechanics, _Celestial Mechanics and Dynamical Astronomy_, Springer, v. 128, n.1., p. 19--59, 2017.

## Disclaimer
We are not offering any license yet, but it will be done under a General Public License in the future. Users are allowed to view and fork this repository under the Terms of Service on GitHub.