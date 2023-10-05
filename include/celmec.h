/**
 * Auxiliar functions in Celestial Mechanics
 * 
 * Author: Vitor M. de Oliveira
 * Date: 11 july 2023
**/

#ifndef CEL_MEC_H
#define CEL_MEC_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "algelin3d.h"

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

/* Struct for celestial bodies (Generalized-Maxwell) */
typedef struct CelestialBody{

	/* body identification */
    char    name[100];  		// name

	/* bulk, orbital and rheology parameters */
	double  mass;       		// mass (sun mass)
	double	R;		    		// radius (km)
	double	lod;	    		// length of day (day)

	double	rg;		    		// moment of inertia factor
	double	J2;		    		// gravity field coefficient
	double	C22;        		// gravity field coefficient
	double	I0;					// mean moment of inertia
	
	double 	obl;	    		// obliquity (deg)
	double	psi;	    		// longitude of the ascending node 
								// of the body equator or precession angle (deg)	
	double	lib;	    		// angle between the ascending node and the
								// lowest principal moment of inertia (deg)
	
	double	a;		    		// semi-major axis (AU)
	double	e;		    		// orbit eccentricity
	double	I;		    		// inclination (deg)
	double	M;		    		// mean anomaly (deg)
	double	w;		    		// argument of periapsis (deg)
	double	Omega;	    		// longitude of the ascending node (deg)

	double	kf;		    		// fluid Love number
	double	Dt;		    		// tidal lag (s)
	double	tau;	    		// Maxwell relaxation time plus tidal lag (yr)

	double	gamma;				// gravitational modulus
	double	alpha;				// elastic modulus
	double	eta;				// viscosity
	double	alpha_0;			// prestress elastic modulus
	int		elements;			// number of voigt elements
	double	*alpha_elements;	// elastic modulus for Voigt elements
	double	*eta_elements;		// viscosity for Voigt elements

	/* orbital options */
	bool	keplerian;			// sets a fixed Keplerian orbit
	bool	orbit_2body;		// considers only the gravitational interaction
								// with the central body

	/* deformation options */
	bool	point_mass;			// point mass
	bool	centrifugal;		// centrifugal force
	bool	tidal;				// tidal force

	/* state variables */
	double	x[3];				// position
	double	x_dot[3];			// velocity
	double	l[3];				// angular momentum
	double	b0_me[5];			// main elements of prestress matrix
	double	u_me[5];			// main elements of u matrix
	double	*bk_me;				// main elements of Voigt elements matrix

	/* body frame */
	double	Y[9];				// transformation (rotation matrix)
								// to the body frame

	/* non-state variables */
	double	omega[3];			// angular velocity
	double	b[9];				// deformation matrix

} cltbdy;

/* prints CelestialBody */
int
print_CelestialBody(cltbdy body);

/* all orbital elements from state vectors */
int
calculate_orbital_elements  (cltbdy *body, 
                             const cltbdy body_ref,
                             const double G);

#endif
