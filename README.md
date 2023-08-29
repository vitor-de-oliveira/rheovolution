# tides
Software under development for numerically investigating the tidal evolution of celestial bodies based on the physical theory developed by Clodoaldo Ragazzo and Lucas Ruiz.

This work is part of a postdoctoral project from the São Paulo Research Foundation (FAPESP - Grants 2021/11306-0 and 2022/12785-1), carried out at the Institute of Mathematics and Statistics of the University of São Paulo (Brazil) and at the Department of Physics of the University of Coimbra (Portugal) with the collaboration of prof. Alexandre Correia.

Author: V. M. de Oliveira

Last update on this file: August 29th 2023

## Important notes
The software simulates the tidal evolution of any number of celestial bodies interacting gravitationally with each other. The rheology model adopted here is the generalized Voigt one, which can be reduced to the Maxwell model. The equations of motion are numerically integrated using a Prince-Dormand Runge-Kutta scheme of 8th and 9th order and fixed stepsize from the GNU Scientific Library.

Up to this point, the simulation assumes that there is no prestress and that the angular velocity of each body is parallel to its biggest moment of inertia.

## How to compile

```sh
make local
```

## How to run

```sh
./TIDES configuration_file.dat
```

The ``configuration_file.dat`` is composed by a few information, such as the name of the input and the output dir and the path to the ``system specs`` and also the ``integrator specs`` file.

It also contains the type of system specs file used. ``1`` stands for variables directly used in the theory, while ``2`` stands for values given in eliptical elements.

An example of all input files can be found in ``example``.

## Main references
Ragazzo C., Ruiz L. S. Viscoelastic tides: models for use in Celestial Mechancis, _Celestial Mechanics and Dynamical Astronomy_, Springer, v. 128, n.1., p. 19--59, 2017.

## Disclaimer
We are not offering any license yet, but it will be done under a General Public License in the future. Users are allowed to view and fork this repository under the Terms of Service on GitHub.