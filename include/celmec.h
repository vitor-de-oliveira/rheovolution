/**
 * Auxiliar functions in Celestial Mechanics
 * 
 * Author: Vitor M. de Oliveira
 * Date: 11 july 2023
**/

#ifndef CEL_MEC_H
#define CEL_MEC_H

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

/* Kepler 3rd law */
double
kepler_period(double m1, double m2, double G, double a);

/* functions for solving Kepler equation */
struct root_params_kepler
{
	double e, M;
};

double
root_function_kepler(double E,
                     void *params);

double
root_derivative_kepler	(double E,
                         void *params);

void
root_fdf_kepler	(double E,
                 void *params,
                 double *y,
                 double *dy);

double
kepler_equation	(const double e,
                 const double M);

#endif
