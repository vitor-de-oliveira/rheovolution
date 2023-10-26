#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include "linear_algebra.h"
#include "celestial_mechanics.h"
#include "tidal_theory.h"

/* field for the Generalized-Voigt rheology */
int
field_GV(double t, 
		 const double y[],
		 double f[],
       	 void *params);

#endif
