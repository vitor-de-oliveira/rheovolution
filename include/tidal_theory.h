#ifndef TIDAL_H
#define TIDAL_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include "linear_algebra.h"
#include "celestial_mechanics.h"

int
calculate_inertia_tensor(double I[9],
						 const double I0,
						 const double b[9]);

int
body_frame_deformation_from_stokes_coefficients	(double B[9],
												 const cltbdy body);

double
parameter_gamma(const double G,
				const double I0, 
				const double R,
				const double kf);

double
parameter_alpha_0	(const double G,
					 const double I0, 
					 const double R,
					 const double kf,
					 const double ks);

double
calculate_c(const cltbdy body);

int
calculate_f_tide(double f_tide[9],
				 const int id,
			 	 const cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G);

int
calculate_f_ps	(double f_ps[9],
			 	 const cltbdy body);

int
calculate_f_rheo(double f_rheo[9],
			 	 const cltbdy body);

int
calculate_g	(double g[9], 
			 const int id,
			 const cltbdy *bodies,
			 const int number_of_bodies,
			 const double G);

int
calculate_f_cent(double f_cent[9],
				 const double omega[3]);

int
calculate_b	(const int id,
			 cltbdy *bodies,
			 const int number_of_bodies,
			 const double G);

int
calculate_l	(cltbdy *body,
			 const int number_of_bodies,
			 const double G);

int
calculate_omega	(const int id,
				 cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G);

int
calculate_total_angular_momentum(double l_total[3],
				 		 		 const cltbdy *bodies,
			 	 		 		 const int number_of_bodies,
			 	 		 		 const double G);

/* semi-major axis growth in the Earth-Moon system */

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
calculate_tau_v_and_tau_from_Rek2(double *tau_v, 
	double *tau, const double Rek2, const double nu,
	const double kf, const double sigma);

int
calculate_k2(double *re, double *im, const double sigma, 
	const double kf, const double tau_v, const double tau);

#endif
