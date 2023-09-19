#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include "algelin3d.h"
#include "celmec.h"

/* field for Generalized-Voigt rheology */
int
field_GV(double t, 
		 const double y[],
		 double f[],
       	 void *params);

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
			 	 const cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G);

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
calculate_inertia_tensor(double I[9],
						 const double I0,
						 const double b[9]);

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
