# tides
Software under development for numerically investigating the tidal evolution of celestial bodies based on the physical theory developed by Clodoaldo Ragazzo and Lucas Ruiz.

This work is part of a postdoctoral project from the São Paulo Research Foundation (FAPESP - Grants 2021/11306-0 and 2022/12785-1), carried out at the Institute of Mathematics and Statistics of the University of São Paulo (Brazil) and at the Department of Physics of the University of Coimbra (Portugal).

Author: V. M. de Oliveira

Last update on this file: June 7th 2023

## Important notes
Up to this point, the software simulates the tidal evolution of one extended body moving under the gravitational field of a point mass. The rheology model for teh extended body implemented here is the generalized Voigt one, which can be reduced to the Maxwell model.
The equations of motion are numerically integrated using a Runge-Kutta scheme of 4th order from the Gnu Scientific Library (link).

## How to compile

As prerequisites...

The makefile was written to compile the program on the group's cluster which uses modules...

there is an option of...uses gcc..

```sh
make local
```

## How to run
The software needs two files ...

```sh
./TIDES system_pararameters_file.txt integrator_pararameters_file.txt
```

A directory named ``output`` will be created at root and a copy of the input files will be written there.

## Main references
Ragazzo (2023) ...

## Disclaimer
We are not offering any license yet, but it will be done under a General Public License in the future. Users are allowed to view and fork this repository under the Terms of Service on GitHub.