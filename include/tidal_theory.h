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
parameter_gamma_0	(const double G,
					 const double I0, 
					 const double R,
					 const double k0);

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
calculate_F_cent_mean	(double F_cent_mean[9],
				 	  	 const double mean_omega);

int
calculate_b	(const int id,
			 cltbdy *bodies,
			 const int number_of_bodies,
			 const double G);

int
calculate_l	(cltbdy *body);

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

#endif
