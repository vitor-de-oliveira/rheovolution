#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include "linear_algebra.h"
#include "celestial_mechanics.h"
#include "tidal_theory.h"
#include "data_processing.h"

/* struct for field parameters */

typedef struct FieldParameters {

	siminf	simulation;
	cltbdy	*bodies;

} fldpar;

/* state vector */

int
mount_state_vector (double **state_vector,
					size_t *dim_state_vector,
					const fldpar params);

int
retrieve_state_vector	(cltbdy **bodies,
						 const double *state_vector,
						 const siminf simulation);

/* field for the generalised Voigt rheology */

int
field(double t, 
	  const double y[],
	  double f[],
      void *params);

#endif
