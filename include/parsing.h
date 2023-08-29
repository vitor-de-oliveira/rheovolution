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
convert_input	(cltbdy	**body,
				 const int number_of_bodies,
				 const double G,
				 const char file[],
				 const char units[]);

#endif
