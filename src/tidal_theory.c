#include "tidal_theory.h"

int
calculate_inertia_tensor(double I[9], const double I0, const double b[9])
{
	double Id[9];
	identity_matrix(Id);
	linear_combination_square_matrix(I, I0, Id, -1.0 * I0, b);
	return 0;
}

double
parameter_gamma(const double G,	const double I0, 
	const double R, const double kf)
{
	return 3.0 * I0 * G / (pow(R, 5.0) * kf);
}

double
calculate_c(const cltbdy body)
{
	return body.gamma + body.alpha_0 + body.alpha;
}

int
calculate_f_tide(double f_tide[9],
				 const int id,
			 	 const cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G)
{
	/* initializing tidal force on body id as null */
	null_matrix(f_tide);

	/* loop over all bodies different from body id */
	for (int i = 0; i < number_of_bodies; i++)
	{
		if (i != id)
		{
			double relative_x[3];
			linear_combination_vector(relative_x,
				1.0, bodies[id].x,
				-1.0, bodies[i].x);
			double x_tensor_x[9];
			tensor_product(x_tensor_x, relative_x, relative_x);
			double Id[9];
			identity_matrix(Id);
			double scaled_id[9];
			scale_square_matrix(scaled_id, 
				norm_squared_vector(relative_x) / 3.0, Id);
			double x_norm_fifth = pow(norm_vector(relative_x), 5.0);
			double f_tide_component[9];
			linear_combination_square_matrix(f_tide_component, 
				 3.0 * G * bodies[i].mass / x_norm_fifth, x_tensor_x,
				-3.0 * G * bodies[i].mass / x_norm_fifth, scaled_id);
			linear_combination_square_matrix(f_tide, 
				1.0, f_tide,
				1.0, f_tide_component);
		}
	}

	return 0;
}

int
calculate_g	(double g[9], 
			 const int id,
			 const cltbdy *bodies,
			 const int number_of_bodies,
			 const double G)
{
	/* for testing */
	// null_matrix(g);
	// return 0;

	/* construct b0, and u matrices */
	double b0[9], u[9];
	construct_traceless_symmetric_matrix(b0, bodies[id].b0_me);
	construct_traceless_symmetric_matrix(u, bodies[id].u_me);

	/* prestress contribution */
	double alpha_0_b0[9];
	scale_square_matrix(alpha_0_b0, bodies[id].alpha_0, b0);

	/* calculate g without Voigt elements */
	linear_combination_square_matrix(g,
		 1.0, alpha_0_b0,
		-1.0, u);

	/* add tidal force if chosen */
	if (bodies[id].tidal == true)
	{
		/* calculate f_tide */
		double f_tide[9];
		calculate_f_tide(f_tide, id, bodies, number_of_bodies, G);
		linear_combination_square_matrix(g,
			1.0, g,
			1.0, f_tide);
	}

	/* add Voigt elements to g */
	double **bk_me_2d_array, **bk;
	if (bodies[id].elements > 0)
	{
		bk_me_2d_array 	= (double **) malloc(bodies[id].elements * sizeof(double *));
		bk 				= (double **) malloc(bodies[id].elements * sizeof(double *));
		for (int i = 0; i < bodies[id].elements; i++)
		{
			bk_me_2d_array[i] = (double *) malloc(5 * sizeof(double));
			for (int j = 0; j < 5; j++)
			{
				bk_me_2d_array[i][j] = bodies[id].bk_me[j + (i*5)];

				/* for testing */
				// printf("bk_me = %f\n",bk_me[i][j]);
			}
			bk[i] = (double *) malloc(9 * sizeof(double));
			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);

			/* for testing */
			// print_square_matrix(bk[i]);

			linear_combination_square_matrix(g,
							1.0, g,
				 bodies[id].alpha, bk[i]);
		}
	}

	/* for testing */
	// for (int i  = 0; i < elements; i++)
	// {
	// 	printf("\nbk_%d = \n", i+1);
	// 	print_square_matrix(bk[i]);
	// }

	/* freeing Voigt elements */
	if (bodies[id].elements > 0)
	{
		for (int i = 0; i < bodies[id].elements; i++) free(bk_me_2d_array[i]);
		free(bk_me_2d_array);
		for (int i = 0; i < bodies[id].elements; i++) free(bk[i]);
		free(bk);
	}

	return 0;
}

int
calculate_f_cent(double f_cent[9], const double omega[3])
{
	double omega_hat[9];
	hat_map(omega_hat, omega);
	double omega_hat_squared[9];
	square_matrix_times_square_matrix(omega_hat_squared,
		omega_hat, omega_hat);
	double trace_omega_hat_squared 
		= trace_square_matrix(omega_hat_squared);
	double Id[9];
	identity_matrix(Id);
	linear_combination_square_matrix(f_cent,
		-1.0, omega_hat_squared,
		trace_omega_hat_squared / 3.0, Id);

	return 0;
}

int
calculate_b	(const int id,
			 cltbdy *bodies,
			 const int number_of_bodies,
			 const double G)
{
	/* for testing */
	// null_matrix(b);
	// return 0;

	if (bodies[id].point_mass == true) /* if body is a point mass, b equals 0 */
	{
		null_matrix(bodies[id].b);
	}
	else if (bodies[id].centrifugal == false && bodies[id].tidal == false) /* if body is not deformable, b equals to b0 */
	{
		construct_traceless_symmetric_matrix(bodies[id].b, bodies[id].b0_me);
	}
	else
	{
		/* calculate g */
		double g[9];
		calculate_g(g, id, bodies, number_of_bodies, G);

		/* calculate c */
		double c = calculate_c(bodies[id]);

		/* calculate b */
		scale_square_matrix(bodies[id].b, 1.0 / c, g);

		/* add centrifugal force if chosen */
		if (bodies[id].centrifugal == true)
		{
			/* calculate centrifugal force */
			double f_cent[9];
			calculate_f_cent(f_cent, bodies[id].omega);
			linear_combination_square_matrix(bodies[id].b,
				1.0, bodies[id].b,
				1.0 / c, f_cent);
		}
	}

	/* for testing */
	// printf("\nb = \n");
	// print_square_matrix(b);
	// printf("\ntilde_x = \n");
	// print_vector(tilde_x);
	// printf("\nb0 = \n");
	// print_square_matrix(b0);
	// printf("\nu = \n");
	// print_square_matrix(u);
	// printf("\nG = %f\n", G);
	// printf("\nm2 = %f\n", m2);
	// printf("\ngamma = %f\n", gamma);
	// printf("\nalpha = %f\n", alpha);
	// printf("\nalpha_0 = %f\n", alpha_0);
	// exit(42);

	// printf("%d\n", id);
	// print_CelestialBody(bodies[id]);

	return 0;
}

int
calculate_l	(cltbdy *body,
			 const int number_of_bodies,
			 const double G)
{
	double I[9];
	calculate_inertia_tensor(I, (*body).I0, (*body).b);
	square_matrix_times_vector((*body).l, I, (*body).omega);
	
	return 0;
}

int
calculate_omega	(const int id,
				 cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G)
{
	/* for testing */
	// scale_vector(omega, 1.0 / I0, l);
	// return 0;

	if (bodies[id].point_mass == true)
	{
		null_vector(bodies[id].omega);
		return 0;
	}

	/* method parameters */
	int 	number_of_iterates = 5;
	double 	max_error = 1e-8;
	double	error = 1.0;
	double 	omega[3], previous_omega[3];
	copy_vector(omega, bodies[id].omega);

	for (int i = 0; i < number_of_iterates; i++)
	// while (error > max_error)
	{
		/* for testing */
		// printf("%d\n", id);
		// print_vector(omega);

		/* store omega previous value */
		copy_vector(previous_omega, omega);

		/* calculate g */
		double g[9];
		calculate_g(g, id, bodies, number_of_bodies, G);

		/* calculate c */
		double c = calculate_c(bodies[id]);

		/* calculate H = 0 */
		double Id[9];
		identity_matrix(Id);
		double aux_H_first_term[9];
		linear_combination_square_matrix(aux_H_first_term,
			 1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (bodies[id].centrifugal == true)
		{
			linear_combination_square_matrix(aux_H_first_term,
				1.0, aux_H_first_term,
				2.0 * norm_squared_vector(omega) / (3.0 * c), Id);
		}
		double H_first_term[3];
		square_matrix_times_vector(H_first_term, aux_H_first_term, omega);
		double H[3];
		linear_combination_vector(H, 
							1.0, H_first_term,
			 -1.0 / bodies[id].I0, bodies[id].l);
		double minus_H[3];
		scale_vector(minus_H, -1.0, H);

		/* calculate DH */
		double DH[9];
		linear_combination_square_matrix(DH,
			1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (bodies[id].centrifugal == true)
		{
			double omega_tensor_omega[9];
			tensor_product(omega_tensor_omega, omega, omega);
			double DH_third_term[9];
			linear_combination_square_matrix(DH_third_term,
				2.0 * norm_squared_vector(omega) / (3.0 * c), Id,
				4.0 / (3.0 * c), omega_tensor_omega);
			linear_combination_square_matrix(DH,
				1.0, DH,
				1.0, DH_third_term);
		}

		/* solving linear equation m*x=b using LU decomposition */
		gsl_matrix_view m
			= gsl_matrix_view_array (DH, 3, 3);
		gsl_vector_view b
			= gsl_vector_view_array (minus_H, 3);
		double omega_minus_previous_omega[] = { 0.0, 0.0, 0.0 };
		gsl_vector_view x
			= gsl_vector_view_array (omega_minus_previous_omega, 3);

		int s;
		gsl_permutation * p = gsl_permutation_alloc (3);
		gsl_linalg_LU_decomp (&m.matrix, p, &s);
		gsl_linalg_LU_solve (&m.matrix, p, &b.vector, &x.vector);

		linear_combination_vector(omega, 
			1.0, omega_minus_previous_omega,
			1.0, previous_omega);

		// error 
		// 	= norm_vector(omega_minus_previous_omega) / norm_vector(previous_omega);
		
		error = norm_vector(omega_minus_previous_omega);

		/* for testing */
		// printf("iter = %d error = %1.5e\n", i+1, error);

		gsl_permutation_free (p);

	}

	if (error > max_error)
	{
		fprintf(stderr, "Error: error higher than the max allowed\n");
		fprintf(stderr, "for omega calculation.\n");
		exit(99);
	}

	copy_vector(bodies[id].omega, omega);

	/* for testing */
	// printf("\nomega = \n");
	// print_vector(omega);
	// exit(42);

	return 0;
}

int
calculate_total_angular_momentum(double l_total[3],
				 		 		 const cltbdy *bodies,
			 	 		 		 const int number_of_bodies,
			 	 		 		 const double G)
{
	null_vector(l_total);

	for (int i = 0; i < number_of_bodies; i++)
	{
		double x_cross_x_dot[3];
		cross_product(x_cross_x_dot, bodies[i].x, bodies[i].x_dot);
		double l_total_component[3];
		linear_combination_vector(l_total_component,
			1.0, bodies[i].l,
			bodies[i].mass, x_cross_x_dot);
		linear_combination_vector(l_total,
			1.0, l_total,
			1.0, l_total_component);
	}

	return 0;
}

double
calibrate_Imk2(const double rate, const double dist, 
	const double m1, const double m2, const double I0, 
	const double R,	const double omega_z, const double G)
{
	double Imk2;

	double alpha = rate;
	double r = dist;
	double a = R;
	double Omega = omega_z;

	double M = m1 + m2;
	double m = (m1 * m2) / M;

	double r2 = pow(r, 2.0);
	double r3 = pow(r, 3.0);
	double r5 = pow(r, 5.0);

	double n = sqrt((G * M) / r3);

	double n2 = pow(n, 2.0);

	double n_p = -(3.0 / 2.0) * sqrt((G * M) / r5);
	double Omega_p = -(m / (2.0 * I0)) * sqrt((G * M) / r);

	double h_1 = m * n * n_p * r2;
	double h_2 = m * n2 * r;
	double h_3 = I0 * Omega * Omega_p;
	double h_4 = (G * m1 * m2) / r2;

	double h = h_1 + h_2 + h_3 + h_4;

	double N_2 = sqrt(5.0 / (4.0 * M_PI * 24.0));

	double a_hat = (G * m2) / (2.0 * N_2 * r3);

	double a5 = pow(a, 5.0);
	double a_hat2 = pow(a_hat, 2.0);

	double omega_SD = 2.0 * (Omega - n); // Semi-diurnal

	double beta = -(5.0 * a5 * a_hat2 * omega_SD) / (32.0 * M_PI * G);

	Imk2 = -(alpha * h) / beta;

	return Imk2;
}

int
calculate_tau_v_and_tau(double tau_v_pair[2], double tau_pair[2],
	const double nu, const double Imk2, const double dist,
	const double m1, const double m2, const double kf, 
	const double omega_z, const double G)
{
	double tau_local_minus;
	double tau_local_plus;
	double tau_v_local_minus;
	double tau_v_local_plus;

	double r = dist;
	double b = Imk2;
	double Omega = omega_z;

	double M = m1 + m2;

	double r3 = pow(r, 3.0);

	double n = sqrt((G * M) / r3);

	double omega_tilde = 2.0 * (Omega - n); // Semi-diurnal

	double kf2 = pow(kf, 2.0);
	double nu2 = pow(nu, 2.0);
	double b2 = pow(b, 2.0);

	tau_local_minus = (-kf * nu) / (2.0 * b * omega_tilde) 
		- (1.0 / (2.0 * omega_tilde)) * sqrt(((kf2 * nu2) / b2) - 4.0);
	tau_local_plus = (-kf * nu) / (2.0 * b * omega_tilde) 
		+ (1.0 / (2.0 * omega_tilde)) * sqrt(((kf2 * nu2) / b2) - 4.0);

	tau_v_local_minus = nu * tau_local_minus; 
	tau_v_local_plus = nu * tau_local_plus; 

	tau_v_pair[0] = tau_v_local_minus;
	tau_v_pair[1] = tau_v_local_plus;
	tau_pair[0] = tau_local_minus;
	tau_pair[1] = tau_local_plus;

	return 0;
}

int
calculate_k2(double *re, double *im, const double sigma, 
	const double kf, const double tau_v, const double tau)
{
	double tau_e = tau - tau_v;
	double sigma2 = sigma * sigma;
	double tau2 = tau * tau;

	*re =  kf * ((1.0 + sigma2 * tau_e * tau) / (1.0 + sigma2 * tau2));
	*im = -kf * ((sigma * tau_v) / (1.0 + sigma2 * tau2));

	return 0;
}
