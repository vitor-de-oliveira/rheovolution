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
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "linear_algebra.h"

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

/* gravitational potential terms */
/* should be calculated using an inertia tensor I[] */
/* measured in the body's reference frame */

double
calculate_rg(const double m, const double R, const double I[9]);

double
calculate_J2(const double m, const double R, const double I[9]);

double
calculate_C22(const double m, const double R, const double I[9]);

double
calculate_S22(const double m, const double R, const double I[9]);

double
calculate_C21(const double m, const double R, const double I[9]);

double
calculate_S21(const double m, const double R, const double I[9]);

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

/* struct for celestial bodies (Generalized-Voigt) */

typedef struct CelestialBody {

	/* body identification */
    char    name[100];  		// name

	/* bulk, orbital and rheology parameters */
	double  mass;       		// mass (sun mass)
	double	R;		    		// radius (km)
	double	lod;	    		// length of day (day)
	double	azi;				// azimuthal angle of angular velocity (deg)
	double	pol;				// polar angle of angular velocity (deg)

	double	I0;					// mean moment of inertia
	double	rg;		    		// moment of inertia factor
	double	J2;		    		// gravity field coefficient
	double	C22;        		// gravity field coefficient
	double	S22;        		// gravity field coefficient
	double	C21;        		// gravity field coefficient
	double	S21;        		// gravity field coefficient

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

	double	k0;					// Love number at zero frequency
	double	Dt;		    		// tidal lag (s)
	double	tau;	    		// Maxwell relaxation time plus tidal lag (yr)

	double	gamma_0;			// gravitational modulus + prestress elastic modulus
	double	alpha;				// elastic modulus
	double	eta;				// viscosity
	int		elements;			// number of Voigt elements
	double	*alpha_elements;	// elastic modulus for Voigt elements
	double	*eta_elements;		// viscosity for Voigt elements

	/* deformation options */
	bool	point_mass;			// point mass
	bool	deformable;			// deformable
	bool	prestress;			// permanent deformation
	bool	centrifugal;		// centrifugal force
	bool	tidal;				// tidal force

	/* state variables (inertial frame) */
	double	x[3];				// position
	double	x_dot[3];			// velocity
	double	l[3];				// angular momentum
	double	p_me[5];			// main elements of prestress matrix
	double	b_eta_me[5];		// main elements of b_eta matrix
	double	*bk_me;				// main elements of Voigt elements matrices

	/* real deformation */
	double	bs_me[5];			// deformation calculated from
								// the stokes coefficients

	/* rotation quaternion */
	double	q[4];				// quaternion which transforms
								// from the body to the inertial frame 

	/* non-state variables (inertial frame) */
	double	omega[3];			// angular velocity
	double	b[9];				// deformation matrix

	/* relative motion (arbitraty frame) */
	double	relative_x[3];
	double	relative_x_dot[3];

} cltbdy;

/* copies CelestialBody */
int
copy_CelestialBody	(cltbdy *body_dest,
					 const cltbdy body_src);

/* creates copied CelestialBody */
cltbdy
create_and_copy_CelestialBody(const cltbdy body_src);

/* prints CelestialBody */
int
print_CelestialBody(cltbdy body);

/* prints state variables CelestialBody */
int
print_state_variables(cltbdy body);

/* all orbital elements from state vectors */
int
calculate_orbital_elements  (cltbdy *body, 
                             const cltbdy body_ref,
                             const double G);

/* initialization of angular velocity vector */

int
initialize_angular_velocity_on_z_axis(cltbdy *body);

int
initialize_angular_velocity(cltbdy *body);

/* body orientation */

int
calculate_obliquity_free_body_from_angular_velocity(cltbdy *body);

int
calculate_obliquity_on_orbit_from_angular_velocity(cltbdy *body);

int
calculate_obliquity_free_body_from_figure_axis_of_solid_frame(cltbdy *body);

int
calculate_obliquity_on_orbit_from_figure_axis_of_solid_frame(cltbdy *body);

double
angle_between_spin_axis_and_figure_axis_of_solid_frame(const cltbdy body);

double
angle_between_spin_axis_and_figure_axis(const cltbdy body);

double
angle_between_spin_axis_and_angular_momentum(const cltbdy body);

/* auxiliary functions */

int
calculate_center_of_mass(double center_of_mass[3],
				 		 const cltbdy *bodies,
			 	 		 const int number_of_bodies,
			 	 		 const double G);

double
largest_time_scale	(const cltbdy *bodies,
			 	 	 const int number_of_bodies,
			 	 	 const double G);

#endif
