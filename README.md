# tides
Software under development for numerically investigating the tidal evolution of celestial bodies based on the physical theory developed by Clodoaldo Ragazzo and Lucas Ruiz.

This work is part of a postdoctoral project from the São Paulo Research Foundation (FAPESP - Grants 2021/11306-0 and 2022/12785-1), carried out at the Institute of Mathematics and Statistics of the University of São Paulo (Brazil) and at the Department of Physics of the University of Coimbra (Portugal) with the collaboration of prof. Alexandre Correia.

Author: V. M. de Oliveira

Last update on this file: October 4th 2023

## Important notes
The software simulates the tidal evolution of any number of celestial bodies interacting gravitationally with each other. The rheology model adopted here is the generalized Voigt one, which can be reduced to the Maxwell model. The equations of motion are numerically integrated using a Prince-Dormand Runge-Kutta scheme of 8th and 9th order and fixed stepsize from the GNU Scientific Library.

Up to this point, the simulation assumes that there is no prestress and that the angular velocity of each body is parallel to its biggest moment of inertia.

## How to compile

```sh
make
```

## How to run

```sh
./TIDES configuration_file.dat  or  make run input=configuration_file.dat
```

The ``configuration_file.dat`` should include the simulation name under ``name``, the location of the input files under ``input_folder``, the location of the output files under ``output_folder``, the name of the files containing the system and the integrator specifications, under ``system_specs`` and ``integrator_specs``, respectively. The number of bodies used can also be given in this file under ``number_of_bodies``. If this value is not passed to the program, the total number of bodies in ``system_specs`` is assumed. The results are written in the ``output_folder`` as ``results_name``.

## Examples

```sh
make && make examples
```

This command compiles the code and runs it with the input files in the folder ``example``, which contains 2 examples: the specs for the Earth-Moon system with no obliquity, zero orbital eccentricity, and the Moon being treated as a point mass, and the specs for the Earth-Moon-Sun system with all objects being treated as point masses.

## Main references
Ragazzo C., Ruiz L. S. Viscoelastic tides: models for use in Celestial Mechanics, _Celestial Mechanics and Dynamical Astronomy_, Springer, v. 128, n.1., p. 19--59, 2017.

## Disclaimer
We are not offering any license yet, but it will be done under a General Public License in the future. Users are allowed to view and fork this repository under the Terms of Service on GitHub.