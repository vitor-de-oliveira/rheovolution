/**
 * Auxiliar functions in Celestial Mechanics
 * 
 * Author: Vitor M. de Oliveira
 * Date: 11 july 2023
**/

#ifndef CEL_MEC_H
#define CEL_MEC_H

#include <math.h>
#include <algelin3d.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

/* Kepler 3rd law */
double
kepler_period(double m1, double m2, double G, double a);

double
kepler_period_only_m1(double m1, double G, double a);

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

/* orbital elements from state vectors */
/* taken from text based on Vallado, 2007 */
/* Fundamentals of Astrodynamics and Applications */
/* Also based on transformations on lib orbital */
/* by RazerM (github) */

double
calculate_semi_major_axis   (const double G,
                             const double m1,
                             const double m2,
                             const double x[],
                             const double v[]);

double
calculate_eccentricity  (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[]);
                     
double
calculate_inclination   (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[]);

double
calculate_true_anomaly  (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[]);

double
calculate_eccentric_anomaly (const double G,
                             const double m1,
                             const double m2,
                             const double x[],
                             const double v[]);

double
calculate_mean_anomaly  (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[]);

double
calculate_argument_of_periapsis (const double G,
                                 const double m1,
                                 const double m2,
                                 const double x[],
                                 const double v[]);

double
calculate_longitude_of_the_ascending_node   (const double G,
                                             const double m1,
                                             const double m2,
                                             const double x[],
                                             const double v[]);

#endif
