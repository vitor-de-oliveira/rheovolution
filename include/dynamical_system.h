#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include "algelin3d.h"

/* Struct for celestial bodies */
typedef struct CelestialBody{
    char    name[100];  // name
	double  mass;       // mass (sun mass)
	double	lod;	    // length of day (day)
	double 	obl;	    // obliquity (deg)
	double	psi;	    // longitude of the ascending node 
						// of the body equator or precession angle (deg)
	double	R;		    // radius (km)
	double	rg;		    // moment of inertia factor
	double	J2;		    // gravity field coefficient
	double	C22;        // gravity field coefficient
	double	lib;	    // angle between the ascending node and the
						// lowest principal moment of inertia (deg)
	double	kf;		    // fluid Love number
	double	Dt;		    // tidal lag (s)
	double	tau;	    // Maxwell relaxation time plus tidal lag (yr)
	double	a;		    // semi-major axis (AU)
	double	e;		    // orbit eccentricity
	double	I;		    // inclination (deg)
	double	M;		    // mean anomaly (deg)
	double	w;		    // argument of periapsis (deg)
	double	Omega;	    // longitude of the ascending node (deg)

	bool	centrifugal; // centrifugal force
	bool	tidal;		 // tidal force
} cltbdy;

/* field for 1 extended body and 1 point mass */
int
field_1EB1PM(double t, const double y[], double f[],
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
parameter_gamma(const double G,	const double I0, 
	const double R, const double kf);

double
calculate_c(const double gamma, const double alpha_0, const double alpha);

int
calculate_f_tide(double f_tide[9], const double G, const double m2,
	const double tilde_x[3]);

int
calculate_g(double g[9], const double G, const double m2, 
	const double alpha_0, const double alpha, const double tilde_x[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool tidal);

int
calculate_f_cent(double f_cent[9], const double omega[3]);

int
calculate_b(double b[9], const double G, const double m2, 
	const double gamma, const double alpha_0, const double alpha,
	const double tilde_x[3], const double omega[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool centrifugal, const bool tidal);

int
calculate_inertia_tensor(double I[9], const double I0, const double b[9]);

int
calculate_l(double l[3], const double I0, 
	const double b[9], const double omega[3]);

int
calculate_omega(double omega[3], const double omega_seed[3], const double G, 
	const double m2, const double I0, const double gamma, const double alpha_0, 
	const double alpha, const double tilde_x[3], const double l[3],
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool centrifugal, const bool tidal);

int
total_angular_momentum(double l_total[3],
	const double m1, const double m2,
	const double tilde_x[3], const double tilde_x_dot[3],
	const double l[3]);

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
