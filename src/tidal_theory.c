#include "tidal_theory.h"

int
calculate_inertia_tensor(double I[9], const double I0, const double b[9])
{
	double Id[9];
	identity_matrix(Id);
	linear_combination_square_matrix(I, I0, Id, -1.0 * I0, b);
	return 0;
}

int
body_frame_deformation_from_stokes_coefficients	(double B[9],
												 const cltbdy body)
{
	double div = 3.0 * body.rg - 2.0 * body.J2;

	B[0] = (body.J2 + 6.0 * body.C22) / div;
	B[1] = (6.0 * body.S22) / div;
	B[2] = (3.0 * body.C21) / div;
	B[3] = B[1];
	B[4] = (body.J2 - 6.0 * body.C22) / div;
	B[5] = (3.0 * body.S21) / div;
	B[6] = B[2];
	B[7] = B[5];
	B[8] = -1.0 * (B[0] + B[4]); // B[8] = (-2.0 * body.J2) / div

	return 0;
}

double
parameter_gamma(const double G,
				const double I0, 
				const double R,
				const double kf)
{
	return 3.0 * I0 * G / (pow(R, 5.0) * kf);
}

double
parameter_alpha_0	(const double G,
					 const double I0, 
					 const double R,
					 const double kf,
					 const double ks)
{
	return (3.0 * I0 * G / pow(R, 5.0)) * (1.0 / ks - 1.0 / kf);
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
			double scaled_Id[9];
			scale_square_matrix(scaled_Id, 
				norm_squared_vector(relative_x) / 3.0, Id);
			double x_norm_fifth = pow(norm_vector(relative_x), 5.0);
			double f_tide_component[9];
			linear_combination_square_matrix(f_tide_component, 
				 3.0 * G * bodies[i].mass / x_norm_fifth, x_tensor_x,
				-3.0 * G * bodies[i].mass / x_norm_fifth, scaled_Id);
			linear_combination_square_matrix(f_tide, 
				1.0, f_tide,
				1.0, f_tide_component);
		}
	}

	return 0;
}

int
calculate_f_ps	(double f_ps[9],
			 	 const cltbdy body)
{
	double b0[9];
	construct_traceless_symmetric_matrix(b0, body.b0_me);
	scale_square_matrix(f_ps, body.alpha_0, b0);

	return 0;
}

int
calculate_f_rheo(double f_rheo[9],
			 	 const cltbdy body)
{
	double u[9];
	construct_traceless_symmetric_matrix(u, body.u_me);
	scale_square_matrix(f_rheo, -1.0, u);
	if (body.elements > 0)
	{
		double bk_me_2d_array[body.elements][5];
		double bk[body.elements][9];
		for (int i = 0; i < body.elements; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				bk_me_2d_array[i][j] = body.bk_me[j + (i*5)];
			}
			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);
			linear_combination_square_matrix(f_rheo,
				1.0, f_rheo,
				body.alpha, bk[i]);
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
	double f_rheo[9], f_ps[9], f_tide[9];
	null_matrix(f_rheo);
	null_matrix(f_ps);
	null_matrix(f_tide);
	
	/* f_rheo (rheology) */
	if (bodies[id].deformable == true)
	{
		calculate_f_rheo(f_rheo, bodies[id]);
	}

	/* f_ps (prestress) */
	if (bodies[id].prestress == true)
	{
		calculate_f_ps(f_ps, bodies[id]);
	}

	/* f_tide (tides) */
	if (bodies[id].tidal == true)
	{
		calculate_f_tide(f_tide, id, bodies, number_of_bodies, G);
	}

	/* calculate g */
	linear_combination_three_square_matrix(g,
		1.0, f_rheo,
		1.0, f_ps,
		1.0, f_tide);

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
calculate_f_cent_static(double f_cent_static[9], const double mean_omega)
{
	double a = (mean_omega * mean_omega) / 3.0;
	double M[] = {1.0, 0.0, 0.0,
				  0.0, 1.0, 0.0,
				  0.0, 0.0, -2.0};
	scale_square_matrix(f_cent_static, a, M);

	return 0;
}

int
calculate_b	(const int id,
			 cltbdy *bodies,
			 const int number_of_bodies,
			 const double G)
{
	/**
	 * if body is a point mass, b equals 0 
	 * if body is not deformable, b equals to b0
	**/
	if (bodies[id].point_mass == true)
	{
		null_matrix(bodies[id].b);
	}
	else if (bodies[id].deformable == false)
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
	/**
	 * if body is a point mass, omega is constant
	**/
	if (bodies[id].point_mass == true)
	{
		return 0;
	}

	/* method parameters */
	int 	number_of_iterates = 5;
	double 	max_error = 1e-8;
	double	error = 1.0;
	double 	omega[3], previous_omega[3];
	copy_vector(omega, bodies[id].omega);

	/* auxiliary variable */
	double Id[9];
	identity_matrix(Id);

	for (int i = 0; i < number_of_iterates; i++)
	// while (error > max_error) // an alternative
	{
		/* store previous value of omega */
		copy_vector(previous_omega, omega);

		double	H[3], DH[9];

		if (bodies[id].deformable == false)
		{
			double b0_i[9];
			construct_traceless_symmetric_matrix(b0_i, 
				bodies[id].b0_me);

			/* calculate H = 0 */
			double aux_H_first_term[9];
			linear_combination_square_matrix(aux_H_first_term,
				 1.0, Id,
				-1.0, b0_i);
			double H_first_term[3];
			square_matrix_times_vector(H_first_term, 
				aux_H_first_term, omega);
			linear_combination_vector(H, 
								 1.0, H_first_term,
				-1.0 / bodies[id].I0, bodies[id].l);

			/* calculate DH */
			linear_combination_square_matrix(DH,
				 1.0, Id,
				-1.0, b0_i);
		}
		else
		{
			/* calculate g */
			double g[9];
			calculate_g(g, id, bodies, number_of_bodies, G);

			/* calculate c */
			double c = calculate_c(bodies[id]);

			/* calculate H = 0 */
			double aux_H_first_term[9];
			linear_combination_square_matrix(aux_H_first_term,
				1.0, Id,
				-1.0 / c, g);

			/* add centrifugal term in H if chosen */
			if (bodies[id].centrifugal == true)
			{
				linear_combination_square_matrix(aux_H_first_term,
					1.0, aux_H_first_term,
					2.0 * norm_squared_vector(omega) / (3.0 * c), Id);
			}
			double H_first_term[3];
			square_matrix_times_vector(H_first_term, 
				aux_H_first_term, omega);
			linear_combination_vector(H, 
								1.0, H_first_term,
				-1.0 / bodies[id].I0, bodies[id].l);

			/* calculate DH */
			linear_combination_square_matrix(DH,
				1.0, Id,
				-1.0 / c, g);

			/* add centrifugal term to DH if chosen */
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
		}

		double minus_H[3];
		scale_vector(minus_H, -1.0, H);

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

		// error = norm_vector(omega_minus_previous_omega) / 
		//			norm_vector(previous_omega);
		
		error = norm_vector(omega_minus_previous_omega);

		/* for testing */
		// printf("iter = %d error = %1.5e\n", i+1, error);

		gsl_permutation_free (p);

	}

	if (error > max_error)
	{
		fprintf(stderr, "Error: error higher than the max allowed\n");
		fprintf(stderr, "for omega calculation in body %d.\n", id + 1);
		exit(99);
	}

	copy_vector(bodies[id].omega, omega);

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
calculate_tau_v_and_tau_from_Rek2(double *tau_v, 
	double *tau, const double Rek2, const double nu,
	const double kf, const double sigma)
{
	*tau = sqrt((kf-Rek2)/(sigma*sigma*(Rek2-kf*(1.0-nu))));
	*tau_v = nu * (*tau); 

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

/***** Implementing gV ******/

double*
convert_parameters_gV_summations_C_and_D(int m,
								  		 double sigma,
								 		 double alpha_elements[],
								 		 double eta_elements[])
{
	double C = 0.0, D = 0.0;
	for (int i = 0; i < m + 1; i++)
	{
		double denominator 
			= alpha_elements[i] * alpha_elements[i]
			+ (sigma * eta_elements[i]) * (sigma * eta_elements[i]);
		C += alpha_elements[i] / denominator;
		D += eta_elements[i] / denominator;
	}

	double *C_and_D = (double *) malloc (2 * sizeof(double));
	C_and_D[0] = C;
	C_and_D[1] = D;

	return C_and_D;
}

double
convert_parameters_gV_components_g (double sigma,
									double tau_b,
								 	double gamma,
								 	double alpha,
								 	double eta,
									double C_sigma,
								 	double D_sigma)
{
	double g = (1.0 + sigma * sigma * eta * D_sigma) * tau_b
				- eta / alpha - eta / gamma - eta * C_sigma;

	return g;
}

int
rosenbrock_f (const gsl_vector * x, void *params,
              gsl_vector * f)
{
  double a = ((struct rparams *) params)->a;
  double b = ((struct rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);

  const double y0 = a * (1 - x0);
  const double y1 = b * (x1 - x0 * x0);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}

int
convert_parameters_gV_f (const gsl_vector * x, void *params,
              gsl_vector * f)
{
	double gamma = ((struct gV_conversion_params *) params)->gamma;
	int m = ((struct gV_conversion_params *) params)->m;
	double *sigma = ((struct gV_conversion_params *) params)->sigma;
	double *tau_a = ((struct gV_conversion_params *) params)->tau_a;
	double *tau_b = ((struct gV_conversion_params *) params)->tau_b;

	double alpha = gsl_vector_get (x, 0);
	double eta = gsl_vector_get (x, 1);

	double *alpha_elements;
	double *eta_elements;

	double *C = (double *) malloc((m + 1) * sizeof(double));
	double *D = (double *) malloc((m + 1) * sizeof(double));

	if (m > 0)
	{
		alpha_elements = (double *) malloc(m * sizeof(double));
		eta_elements = (double *) malloc(m * sizeof(double));

		for (int i = 0; i < 2 * m; i = i + 2)
		{
			alpha_elements[i] = gsl_vector_get (x, i);
			eta_elements[i] = gsl_vector_get (x, i+1);
		}

		for (int i = 0; i < m + 1; i++)
		{
			double *C_and_D;
			
			C_and_D = convert_parameters_gV_summations_C_and_D(m,
					sigma[i], alpha_elements, eta_elements);

			C[i] = C_and_D[0];
			D[i] = C_and_D[1];

			free(C_and_D);
		}
	}
	else
	{
		C[0] = 0.0;
		D[0] = 0.0;
	}

	for (int i = 0; i < 2 * (m + 1); i = i + 2)
	{
		double component_f = 
			(1.0 + sigma[i] * sigma[i] * eta * D[i]) * tau_a[i]
				- eta / alpha - eta * C[i];

		double component_g = 
			(1.0 + sigma[i] * sigma[i] * eta * D[i]) * tau_b[i]
				- eta / alpha - eta / gamma - eta * C[i];

		gsl_vector_set (f, i, component_f);
		gsl_vector_set (f, i + (m + 1), component_g);
	}

	if (m > 0)
	{
		free(alpha_elements);
		free(eta_elements);
	}
	free(C);
	free(D);

	return GSL_SUCCESS;
}

double*
convert_parameters_gV_partial_derivatives_C_D_alpha_k_eta_k	(int m,
															 double sigma,
															 double alpha_k,
															 double eta_k)
{
	double del_C_del_alpha_k;
	double del_C_del_eta_k;
	double del_D_del_alpha_k;
	double del_D_del_eta_k;

	double denominator 
		= pow(sigma * sigma * eta_k * eta_k + alpha_k * alpha_k, 2.0);
	double numerator_1 
		= sigma * sigma * eta_k * eta_k - alpha_k * alpha_k;
	double numerator_2 
		= - 2.0 * sigma * sigma * alpha_k * eta_k;
	double numerator_3 
		= - 2.0 * alpha_k * eta_k;
	double numerator_4 
		= -1.0 * numerator_1;

	del_C_del_alpha_k = numerator_1 / denominator;
	del_C_del_eta_k = numerator_2 / denominator;
	del_D_del_alpha_k = numerator_3 / denominator;
	del_D_del_eta_k = numerator_4 / denominator;

	double *derivatives_C_D_alphak_etak = (double *) malloc (4 * sizeof(double));
	derivatives_C_D_alphak_etak[0] = del_C_del_alpha_k;
	derivatives_C_D_alphak_etak[1] = del_C_del_eta_k;
	derivatives_C_D_alphak_etak[2] = del_D_del_alpha_k;
	derivatives_C_D_alphak_etak[3] = del_D_del_eta_k;

	return derivatives_C_D_alphak_etak;
}

int
rosenbrock_df (const gsl_vector * x, void *params,
               gsl_matrix * J)
{
  const double a = ((struct rparams *) params)->a;
  const double b = ((struct rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);

  const double df00 = -a;
  const double df01 = 0;
  const double df10 = -2 * b  * x0;
  const double df11 = b;

  gsl_matrix_set (J, 0, 0, df00);
  gsl_matrix_set (J, 0, 1, df01);
  gsl_matrix_set (J, 1, 0, df10);
  gsl_matrix_set (J, 1, 1, df11);

  return GSL_SUCCESS;
}

int
convert_parameters_gV_df (const gsl_vector * x, void *params,
               gsl_matrix * J)
{
	double gamma = ((struct gV_conversion_params *) params)->gamma;
	int m = ((struct gV_conversion_params *) params)->m;
	double *sigma = ((struct gV_conversion_params *) params)->sigma;
	double *tau_a = ((struct gV_conversion_params *) params)->tau_a;
	double *tau_b = ((struct gV_conversion_params *) params)->tau_b;

	double alpha = gsl_vector_get (x, 0);
	double eta = gsl_vector_get (x, 1);

	double *C = (double *) malloc((m + 1) * sizeof(double));
	double *D = (double *) malloc((m + 1) * sizeof(double));

	double *alpha_elements;
	double *eta_elements;

	double **del_C_del_alpha_k, **del_C_del_eta_k;
	double **del_D_del_alpha_k, **del_D_del_eta_k;

	if (m > 0)
	{
		alpha_elements = (double *) malloc(m * sizeof(double));
		eta_elements = (double *) malloc(m * sizeof(double));

		for (int i = 0; i < 2 * m; i = i + 2)
		{
			alpha_elements[i] = gsl_vector_get (x, i);
			eta_elements[i] = gsl_vector_get (x, i+1);
		}

		for (int i = 0; i < m + 1; i++)
		{
			double *C_and_D;
			
			C_and_D = convert_parameters_gV_summations_C_and_D(m,
					sigma[i], alpha_elements, eta_elements);

			C[i] = C_and_D[0];
			D[i] = C_and_D[1];

			free(C_and_D);
		}

		del_C_del_alpha_k = (double **) malloc((m + 1) * sizeof(double *));
		del_C_del_eta_k = (double **) malloc((m + 1) * sizeof(double *));
		del_D_del_alpha_k = (double **) malloc((m + 1) * sizeof(double *));
		del_D_del_eta_k = (double **) malloc((m + 1) * sizeof(double *));
		for (int i = 0; i < m + 1; i++)
		{
			del_C_del_alpha_k[i] = (double *) malloc(m * sizeof(double));
			del_C_del_eta_k[i] = (double *) malloc(m * sizeof(double));
			del_D_del_alpha_k[i] = (double *) malloc(m * sizeof(double));
			del_D_del_eta_k[i] = (double *) malloc(m * sizeof(double));
		}


		for (int i = 0; i < m + 1; i++)
		{
			for (int j = 0; j < m; j++)
			{
				double *partial_C_D_alpha_k_eta_k;

				partial_C_D_alpha_k_eta_k 
					= convert_parameters_gV_partial_derivatives_C_D_alpha_k_eta_k(m,
						sigma[i], alpha_elements[j], eta_elements[j]);

				del_C_del_alpha_k[i][j]	= partial_C_D_alpha_k_eta_k[0];
				del_C_del_eta_k[i][j] 	= partial_C_D_alpha_k_eta_k[1];
				del_D_del_alpha_k[i][j] = partial_C_D_alpha_k_eta_k[2];
				del_D_del_eta_k[i][j] 	= partial_C_D_alpha_k_eta_k[3];

				free(partial_C_D_alpha_k_eta_k);
			}
		}

	}
	else
	{
		C[0] = 0.0;
		D[0] = 0.0;
	}

	for (int i = 0; i < 2 * (m + 1); i = i + 2)
	{
		double component_del_f_del_alpha
			= eta / (alpha * alpha);

		double component_del_f_del_eta
			= sigma[i] * sigma[i] * D[i] * tau_a[i]
				- 1.0 / alpha - C[i];

		double component_del_g_del_alpha
			= component_del_f_del_alpha;

		double component_del_g_del_eta
			= sigma[i] * sigma[i] * D[i] * tau_b[i]
				- 1.0 / alpha - 1.0 / gamma - C[i];		

		gsl_matrix_set (J, i, 0, component_del_f_del_alpha);
		gsl_matrix_set (J, i, 1, component_del_f_del_eta);
		gsl_matrix_set (J, i + (m + 1), 0, component_del_g_del_alpha);
		gsl_matrix_set (J, i + (m + 1), 1, component_del_g_del_eta);

		if (m > 0)
		{
			for (int j = 0; j < 2 * m; j = j + 2)
			{
				double component_del_f_del_alpha_k
					= sigma[i] * sigma[i] * eta * tau_a[i] 
						* del_D_del_alpha_k[i][j]
						- eta * del_C_del_alpha_k[i][j];

				double component_del_f_del_eta_k
					= sigma[i] * sigma[i] * eta * tau_a[i] 
						* del_D_del_eta_k[i][j]
						- eta * del_C_del_eta_k[i][j];

				double component_del_g_del_alpha_k
					= sigma[i] * sigma[i] * eta * tau_b[i] 
						* del_D_del_alpha_k[i][j]
						- eta * del_C_del_alpha_k[i][j];

				double component_del_g_del_eta_k
					= sigma[i] * sigma[i] * eta * tau_b[i]
						* del_D_del_eta_k[i][j]
						- eta * del_C_del_eta_k[i][j];

				gsl_matrix_set (J, i, 2 + j, component_del_f_del_alpha_k);
				gsl_matrix_set (J, i, 2 + j + 1, component_del_f_del_eta_k);
				gsl_matrix_set (J, i + (m + 1), 2 + j, component_del_g_del_alpha_k);
				gsl_matrix_set (J, i + (m + 1), 2 + j + 1, component_del_g_del_eta_k);

			}
		}
	}

	if (m > 0)
	{
		free(alpha_elements);
		free(eta_elements);
		for (int i = 0; i < m + 1; i++)
		{
			free(del_C_del_alpha_k[i]);
			free(del_C_del_eta_k[i]);
			free(del_D_del_alpha_k[i]);
			free(del_D_del_eta_k[i]);		
		}
		free(del_C_del_alpha_k);
		free(del_C_del_eta_k);
		free(del_D_del_alpha_k);
		free(del_D_del_eta_k);
	}
	free(C);
	free(D);

	return GSL_SUCCESS;
}

int
rosenbrock_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * J)
{
  rosenbrock_f (x, params, f);
  rosenbrock_df (x, params, J);

  return GSL_SUCCESS;
}

int
print_state (size_t iter, gsl_multiroot_fdfsolver * s)
{
	printf ("iter = %3lu x = % .3f % .3f "
			"f(x) = % .3e % .3e\n",
			iter,
			gsl_vector_get (s->x, 0),
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->f, 0),
			gsl_vector_get (s->f, 1));

	return 0;
}

int
convert_parameters_gV()
{
	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;

	int status;
	size_t iter = 0;

	const size_t n = 2;
	struct rparams p = {1.0, 10.0};
	gsl_multiroot_function_fdf f = {&rosenbrock_f,
									&rosenbrock_df,
									&rosenbrock_fdf,
									n, &p};

	double x_init[2] = {-10.0, -5.0};
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set (x, 0, x_init[0]);
	gsl_vector_set (x, 1, x_init[1]);

	T = gsl_multiroot_fdfsolver_gnewton;
	s = gsl_multiroot_fdfsolver_alloc (T, n);
	gsl_multiroot_fdfsolver_set (s, &f, x);

	print_state (iter, s);

	do
		{
		iter++;

		status = gsl_multiroot_fdfsolver_iterate (s);

		print_state (iter, s);

		if (status)
			break;

		status = gsl_multiroot_test_residual (s->f, 1e-7);
		}
	while (status == GSL_CONTINUE && iter < 1000);

	printf ("status = %s\n", gsl_strerror (status));

	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free (x);

	return 0;
}
