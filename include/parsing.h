/**
 * Library for parsing the input file
*/

#ifndef PSG_H
#define PSG_H

#define _GNU_SOURCE // getline
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h> // ssize_t

#include "algelin3d.h"
#include "celmec.h"
#include "dynamical_system.h"

int
count_columns(const char *s);

int
read_system_type_2	(cltbdy **body,
					 int number_of_bodies,
					 const char file[]);

int
convert_input	(double *m1, double *m2, double *I0, double *R,
				 double *kf, double omega[],
				 double *alpha, double *eta,
				 double tilde_x[], double tilde_x_dot[],
				 bool *centrifugal, bool *tidal,
				 const double G,
				 const char file[],
				 const char units[],
				 const int number_of_bodies);

#endif