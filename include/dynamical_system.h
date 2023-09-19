#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include "algelin3d.h"

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

	/* non-state variables */
	double	omega[3];			// angular velocity
	double	b[9];				// deformation matrix

} cltbdy;

/* prints CelestialBody */
int
print_CelestialBody(cltbdy body);

/* field for Generalized-Voigt rheology */
int
field_GV(double t, 
		 const double y[],
		 double f[],
       	 void *params);

int
hat_map(double x_hat[9], const double x[3]);

int
construct_traceless_symmetric_matrix(double M[9], 
	const double M_main_elements[5]);

int
get_main_elements_traceless_symmetric_matrix(double M_main_elements[5], 
	const double M[9]);

double
parameter_gamma(const double G,
				const double I0, 
				const double R,
				const double kf);

double
calculate_c(const cltbdy body);

int
calculate_f_tide(double f_tide[9],
				 const int id,
			 	 const cltbdy *body,
			 	 const int number_of_bodies,
			 	 const double G);

int
calculate_g	(double g[9], 
			 const int id,
			 const cltbdy *body,
			 const int number_of_bodies,
			 const double G);

int
calculate_f_cent(double f_cent[9],
				 const double omega[3]);

int
calculate_b	(const int id,
			 cltbdy *body,
			 const int number_of_bodies,
			 const double G);
int
calculate_inertia_tensor(double I[9],
						 const double I0,
						 const double b[9]);

int
calculate_l	(const int id,
			 cltbdy *body,
			 const int number_of_bodies,
			 const double G);

int
calculate_omega	(const int id,
				 cltbdy *body,
			 	 const int number_of_bodies,
			 	 const double G);

int
calculate_total_angular_momentum(double l_total[3],
				 		 		 const cltbdy *bodies,
			 	 		 		 const int number_of_bodies,
			 	 		 		 const double G);

int
calculate_center_of_mass(double center_of_mass[3],
				 		 const cltbdy *bodies,
			 	 		 const int number_of_bodies,
			 	 		 const double G);

double
calculate_J2(const double m, const double R, const double I[9]);

double
calculate_C22(const double m, const double R, const double I[9]);

double
calibrate_Imk2(const double rate, const double dist, 
	const double m1, const double m2, const double I0, 
	const double R,	const double omega_z, const double G);

int
calculate_tau_v_and_tau(double tau_v_pair[2], double tau_pair[2],
	const double nu, const double Imk2, const double dist,
	const double m1, const double m2, const double kf, 
	const double omega_z, const double G);

int
calculate_k2(double *re, double *im, const double sigma, 
	const double kf, const double tau_v, const double tau);

#endif
