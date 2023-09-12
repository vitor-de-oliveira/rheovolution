#include "dynamical_system.h"

int
field_1EB1PM(double t, const double y[], double f[],
       		 void *params)
{
	(void)(t);

	/* preparing parameters */

	double 	*par = (double *)params;

	double	omega_seed[] 	= { par[0], par[1], par[2] };
	double 	G				= par[3];
	double 	m1		 		= par[4];
	double 	m2		 		= par[5];
	double 	I0 		 		= par[6];
	double	gamma	 		= par[7];
	double 	alpha 	 		= par[8];
	double 	eta 	 		= par[9];
	double	alpha_0	 		= par[10];
	int 	elements 		= (int) par[11];
	double	*alpha_elements, *eta_elements;
	if (elements > 0)
	{
		alpha_elements 	= (double *) malloc(elements * sizeof(double));
		eta_elements 	= (double *) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			alpha_elements[i] 	= par[12 + (2*i)];
			eta_elements[i]		= par[13 + (2*i)];

			/* for testing */
			// printf("alpha_%d = %f\n", i, alpha_elements[i]);
			// printf("eta_%d = %f\n", i, eta_elements[i]);
		}
	}
	bool	centrifugal		= (bool) par[12 + (elements * 2)];
	bool	tidal			= (bool) par[13 + (elements * 2)];

	/* preparing variables */

	double tilde_x[3], tilde_x_dot[3], l[3];
	double b0_me[9], u_me[9];
	double *bk_me = *(&bk_me), **bk_me_2d_array;

	for (int i = 0; i < 3; i++)	tilde_x[i] 		= y[0 + i];
	for (int i = 0; i < 3; i++) tilde_x_dot[i] 	= y[3 + i];
	for (int i = 0; i < 3; i++) l[i] 			= y[6 + i];
	for (int i = 0; i < 5; i++) b0_me[i] 		= y[9 + i];
	for (int i = 0; i < 5; i++) u_me[i] 		= y[14 + i];
	if (elements > 0)
	{
		bk_me = (double *) malloc(elements * 5 * sizeof(double));
		for (int i  = 0; i < elements * 5; i++)
		{
			bk_me[i] = y[19 + i];
		}
		bk_me_2d_array = (double **) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			bk_me_2d_array[i] = (double *) malloc(5 * sizeof(double));
			for (int j = 0; j < 5; j++)
			{
				bk_me_2d_array[i][j] = y[19 + j + (i*5)];

				/* for testing */
				// printf("bk_me_2d_array = %f\n",bk_me_2d_array[i][j]);
			}
		}
	}

	double b0[9], u[9];
	construct_traceless_symmetric_matrix(b0, b0_me);
	construct_traceless_symmetric_matrix(u, u_me);

	double **bk = *(&(*(&bk)));
	if (elements > 0)
	{
		bk = (double **) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			bk[i] = (double *) malloc(9 * sizeof(double));
			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);

			/* for testing */
			// print_square_matrix(bk[i]);
		}
	}

	/* calculate omega and b */
	double omega[3], b[9];
	// copy_vector(omega, omega_seed); // for testing
	calculate_omega(omega, omega_seed, G, m2, I0, gamma, alpha_0, 
		alpha, tilde_x, l, b0_me, u_me, elements, bk_me, centrifugal, tidal);
	calculate_b(b, G, m2, gamma, alpha_0, alpha,
		tilde_x, omega, b0_me, u_me, elements, bk_me,
		centrifugal, tidal);
	
	/* for testing */
	// printf("omega inside = \n");
	// print_vector(omega);
	// printf("b = \n");
	// print_square_matrix(b);
	// exit(42);
	// null_matrix(b);
	// null_matrix(omega);
	// double b_me[5];
	// for (int i = 0; i < 5; i++)
	// {
	// 	b_me[i] = ((double) i) * 0.00000000001;
	// }
	// construct_traceless_symmetric_matrix(b, b_me);

	double omega_hat[9];
	hat_map(omega_hat, omega);

	/* useful definitions */

	double minus_G_times_total_mass = -1.0 * G * (m1 + m2);

	double tilde_x_norm 		= norm_vector(tilde_x);
	double tilde_x_norm_cube 	= pow(tilde_x_norm, 3.0);
	double tilde_x_norm_fifth 	= pow(tilde_x_norm, 5.0);

	double tau = eta / alpha;
	double *tau_elements;
	if (elements > 0)
	{
		tau_elements = (double *) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			tau_elements[i] = eta_elements[i] / alpha_elements[i];
		}			
	}

	double bx[3];
	square_matrix_times_vector(bx, b, tilde_x);
	double x_cross_bx[3];
	cross_product(x_cross_bx, tilde_x, bx);

	double lambda[9];
	linear_combination_square_matrix(lambda, 1.0, u, alpha, b);
	for (int i = 0; i < elements; i++)
	{
		linear_combination_square_matrix(lambda, 1.0, lambda, -1.0*alpha, bk[i]);
	}

	/* calculating components */

	// tilde_x component

	double component_tilde_x[] = { 0.0, 0.0, 0.0 };
	copy_vector (component_tilde_x, tilde_x_dot);

	// tilde_x_dot component

	double component_tilde_x_dot_1st_term[] = { 0.0, 0.0, 0.0 };
	scale_vector (component_tilde_x_dot_1st_term, 
		1.0 / tilde_x_norm_cube, tilde_x);

	double component_tilde_x_dot_2nd_term[] = { 0.0, 0.0, 0.0 };
	double tilde_x_norm_seventh = pow(tilde_x_norm, 7.0);
	double bx_dot_x = dot_product(bx, tilde_x);
	scale_vector (component_tilde_x_dot_2nd_term, 
		(15. * I0 * bx_dot_x) / (2. * m1 * tilde_x_norm_seventh), tilde_x);

	double component_tilde_x_dot_3rd_term[] = { 0.0, 0.0, 0.0 };
	scale_vector (component_tilde_x_dot_3rd_term, 
		(-3.0 * I0) / (m1 * tilde_x_norm_fifth), bx);

	double component_tilde_x_dot[] = { 0.0, 0.0, 0.0 };
	linear_combination_three_vector(component_tilde_x_dot,
		minus_G_times_total_mass, component_tilde_x_dot_1st_term, 
		minus_G_times_total_mass, component_tilde_x_dot_2nd_term, 
		minus_G_times_total_mass, component_tilde_x_dot_3rd_term);

	// l component

	double component_l[] = { 0.0, 0.0, 0.0 };
	scale_vector (component_l, 
		-3.0 * G * m2 * I0 / tilde_x_norm_fifth, x_cross_bx);

	// b0 component

	double component_b0[] = { 0.0, 0.0, 0.0,
							  0.0, 0.0, 0.0,
							  0.0, 0.0, 0.0 };
	commutator(component_b0, omega_hat, b0);
	double component_b0_me[] = { 0.0, 0.0, 0.0,
							          0.0, 0.0 };
	get_main_elements_traceless_symmetric_matrix(component_b0_me, component_b0);

	// u component

	double omega_hat_comm_u[9];
	commutator(omega_hat_comm_u, omega_hat, u);
	double component_u[] = { 0.0, 0.0, 0.0,
							 0.0, 0.0, 0.0,
							 0.0, 0.0, 0.0 };
	linear_combination_square_matrix(component_u, 
		1.0, omega_hat_comm_u, -1.0 / tau, lambda);
	double component_u_me[] = { 0.0, 0.0, 0.0,
							         0.0, 0.0 };
	get_main_elements_traceless_symmetric_matrix(component_u_me, component_u);

	/* for testing */
	// printf("\nomega_hat_comm_u = \n");
	// print_square_matrix(omega_hat_comm_u);
	// printf("\ntau = \n");
	// printf("%f\n", tau);
	// printf("\nlambda = \n");
	// print_square_matrix(lambda);
	// printf("\ncomponent_u = \n");
	// print_square_matrix(component_u);	
	// printf("\ncomponent_u_me = \n");
	// printf("%1.10e %1.10e %1.10e %1.10e %1.10e\n",
	// 	component_u_me[0], component_u_me[1], component_u_me[2],
	// 	component_u_me[3], component_u_me[4]);	
	// printf("\nu_me = \n");
	// printf("%1.10e %1.10e %1.10e %1.10e %1.10e\n",
	// 	u_me[0], u_me[1], u_me[2], u_me[3], u_me[4]);	
	// exit(42);

	// bk components

	double **component_bk_me = *(&(*(&component_bk_me)));
	if (elements > 0)
	{
		component_bk_me = (double **) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			component_bk_me[i] = (double *) malloc(5 * sizeof(double));

			double omega_hat_comm_bk[9];
			commutator(omega_hat_comm_bk, omega_hat, bk[i]);
			double minus_bk_over_tau_elements[9];
			scale_square_matrix(minus_bk_over_tau_elements, 
				-1.0 / tau_elements[i], bk[i]);
			double lambda_over_eta_elements[9];
			scale_square_matrix(lambda_over_eta_elements,
				1.0 / eta_elements[i], lambda);

			/* for testing */
			// printf("\ntau_%d = %f\n", i, tau_elements[i]);
			// printf("\nminus_1_over_tau_%d = %f\n", i, -1.0 / tau_elements[i]);
			// printf("\nomega_hat_comm_bk = \n");
			// print_square_matrix(omega_hat_comm_bk);
			// printf("\nminus_bk_over_tau_elements = \n");
			// print_square_matrix(minus_bk_over_tau_elements);
			// printf("\nlambda_over_eta_elements = \n");
			// print_square_matrix(lambda_over_eta_elements);

			double component_bk[] = { 0.0, 0.0, 0.0,
							 		  0.0, 0.0, 0.0,
							          0.0, 0.0, 0.0 };
			linear_combination_three_square_matrix(component_bk,
				1.0, omega_hat_comm_bk,
				1.0, minus_bk_over_tau_elements,
				1.0, lambda_over_eta_elements);

			get_main_elements_traceless_symmetric_matrix(component_bk_me[i],
				component_bk);

			/* for testing */
			// printf("\nalpha_%d = %f\n", i, alpha_elements[i]);
			// printf("\ntau_%d = %f\n", i, tau_elements[i]);
			// printf("\ncomponent_bk = \n");
			// print_square_matrix(component_bk);
		}
	}

	/* writing components */	

	for (int i = 0; i < 3; i++) f[i] 		= component_tilde_x[i];
	for (int i = 0; i < 3; i++) f[3 + i] 	= component_tilde_x_dot[i];
	for (int i = 0; i < 3; i++) f[6 + i] 	= component_l[i];
	for (int i = 0; i < 5; i++) f[9 + i]	= component_b0_me[i];
	for (int i = 0; i < 5; i++) f[14 + i]	= component_u_me[i];
	for (int i = 0; i < elements; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			f[19 + j + (i*5)] = component_bk_me[i][j];
		}
	}

	/* for testing */
	// printf("tilde_x = \n");
	// print_vector(tilde_x);
	// printf("\ntilde_x_dot = \n");
	// print_vector(tilde_x_dot);
	// printf("\nl = \n");
	// print_vector(l);
	// printf("\nb0 = \n");
	// print_square_matrix(b0);
	// printf("\nb = \n");
	// print_square_matrix(b);
	// printf("\nu = \n");
	// print_square_matrix(u);
	// for (int i  = 0; i < elements; i++)
	// {
	// 	printf("\nbk_%d = \n", i+1);
	// 	print_square_matrix(bk[i]);
	// }
	// printf("\nomega = \n");
	// print_vector(omega);
	// printf("\nG = %e\n", G);
	// printf("\nm1 = %e\n", m1);
	// printf("\nm2 = %e\n", m2);
	// printf("\nI0 = %e\n", I0);
	// printf("\ngamma = %e\n", gamma);
	// printf("\nalpha = %e\n", alpha);
	// printf("\neta = %e\n", eta);
	// printf("\nalpha_0 = %e\n", alpha_0);
	// for (int i  = 0; i < elements; i++)
	// {
	// 	printf("\nalpha_%d = %e\n", i+1, alpha_elements[i]);
	// 	printf("\neta%d = %e\n", i+1, eta_elements[i]);
	// }
	// printf("\ncomponent_tilde_x = \n");
	// print_vector(component_tilde_x);
	// printf("\ncomponent_tilde_x_dot = \n");
	// print_vector(component_tilde_x_dot);
	// printf("\ncomponent_l = \n");
	// print_vector(component_l);
	// printf("\ncomponent_u = \n");
	// print_square_matrix(component_u);	
	// printf("\ncomponent_b0 = \n");
	// print_square_matrix(component_b0);
	// printf("\nlambda = \n");
	// print_square_matrix(lambda);
	// exit(42);
	
	/* freeing Voigt elements */
	if (elements > 0)
	{
		free(alpha_elements);
		free(eta_elements);
		free(tau_elements);
		free(bk_me);
		for (int i = 0; i < elements; i++) free(bk_me_2d_array[i]);
		free(bk_me_2d_array);
		for (int i = 0; i < elements; i++) free(bk[i]);
		free(bk);
		for (int i = 0; i < elements; i++) free(component_bk_me[i]);
		free(component_bk_me);
	}

	return GSL_SUCCESS;
}

int
hat_map(double x_hat[9], const double x[3])
{
	x_hat[0] = 0.0;
	x_hat[1] = -1.0 * x[2];
	x_hat[2] = x[1];
	x_hat[3] = x[2];
	x_hat[4] = 0.0;
	x_hat[5] = -1.0 * x[0];
	x_hat[6] = -1.0 * x[1];
	x_hat[7] = x[0];
	x_hat[8] = 0.0;

	return 0;
}

int
check_map(double x[3], const double x_hat[9])
{
	x[0] = -1.0 * x_hat[5];
	x[1] = x_hat[2];
	x[2] = -1.0 * x_hat[1];

	return 0;
}

int
construct_traceless_symmetric_matrix(double M[9], 
	const double M_main_elements[5])
{
	double M_11 = M_main_elements[0];
	double M_12 = M_main_elements[1];
	double M_13 = M_main_elements[2];
	double M_22 = M_main_elements[3];
	double M_23 = M_main_elements[4];

	M[0] = M_11;
	M[1] = M_12;
	M[2] = M_13;
	M[3] = M_12;
	M[4] = M_22;
	M[5] = M_23;
	M[6] = M_13;
	M[7] = M_23;
	M[8] = -1.0 * (M_11 + M_22);

	return 0;
}

int
get_main_elements_traceless_symmetric_matrix(double M_main_elements[5], 
	const double M[9])
{
	double M_11 = M[0];
	double M_12 = M[1];
	double M_13 = M[2];
	double M_22 = M[4];
	double M_23 = M[5];

	M_main_elements[0] = M_11;
	M_main_elements[1] = M_12;
	M_main_elements[2] = M_13;
	M_main_elements[3] = M_22;
	M_main_elements[4] = M_23;

	return 0;
}

double
parameter_gamma(const double G,	const double I0, 
	const double R, const double kf)
{
	return 3.0 * I0 * G / (pow(R, 5.0) * kf);
}

double
calculate_c(const double gamma, const double alpha_0, const double alpha)
{
	return gamma + alpha_0 + alpha;
}

int
calculate_f_tide(double f_tide[9], const double G, const double m2,
	const double tilde_x[3])
{
	double x_tensor_x[9];
	tensor_product(x_tensor_x, tilde_x, tilde_x);
	double Id[9];
	identity_matrix(Id);
	double scaled_id[9];
	scale_square_matrix(scaled_id, 
		norm_squared_vector(tilde_x) / 3.0, Id);
	double tilde_x_norm_fifth = pow(norm_vector(tilde_x), 5.0);
	linear_combination_square_matrix(f_tide, 
		 3.0 * G * m2 / tilde_x_norm_fifth, x_tensor_x,
		-3.0 * G * m2 / tilde_x_norm_fifth, scaled_id);

	return 0;
}

int
calculate_g(double g[9], const double G, const double m2, 
	const double alpha_0, const double alpha, const double tilde_x[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool tidal)
{
	/* for testing */
	// null_matrix(g);
	// return 0;

	/* construct b0, and u matrices */
	double b0[9], u[9];
	construct_traceless_symmetric_matrix(b0, b0_me);
	construct_traceless_symmetric_matrix(u, u_me);

	/* prestress contribution */
	double alpha_0_b0[9];
	scale_square_matrix(alpha_0_b0, alpha_0, b0);

	/* calculate g without Voigt elements */
	linear_combination_square_matrix(g,
		1.0, alpha_0_b0,
		-1.0, u);

	/* add tidal force if chosen */
	if (tidal == true)
	{
		/* calculate f_tide */
		double f_tide[9];
		calculate_f_tide(f_tide, G, m2, tilde_x);
		linear_combination_square_matrix(g,
			1.0, g,
			1.0, f_tide);
	}

	/* add Voigt elements to g */
	double **bk_me_2d_array, **bk;
	if (elements > 0)
	{
		bk_me_2d_array 	= (double **) malloc(elements * sizeof(double));
		bk 				= (double **) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			bk_me_2d_array[i] = (double *) malloc(5 * sizeof(double));
			for (int j = 0; j < 5; j++)
			{
				bk_me_2d_array[i][j] = bk_me[j + (i*5)];

				/* for testing */
				// printf("bk_me = %f\n",bk_me[i][j]);
			}
			bk[i] = (double *) malloc(9 * sizeof(double));
			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);

			/* for testing */
			// print_square_matrix(bk[i]);

			linear_combination_square_matrix(g,
				1.0,   g,
				alpha, bk[i]);
		}
	}

	/* for testing */
	// for (int i  = 0; i < elements; i++)
	// {
	// 	printf("\nbk_%d = \n", i+1);
	// 	print_square_matrix(bk[i]);
	// }

	/* freeing Voigt elements */
	if (elements > 0)
	{
		for (int i = 0; i < elements; i++) free(bk_me_2d_array[i]);
		free(bk_me_2d_array);
		for (int i = 0; i < elements; i++) free(bk[i]);
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
calculate_b(double b[9], const double G, const double m2, 
	const double gamma, const double alpha_0, const double alpha,
	const double tilde_x[3], const double omega[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool centrifugal, const bool tidal)
{
	/* for testing */
	// null_matrix(b);
	// return 0;

	/* if body is not deformable, b equals to b0 */
	if (centrifugal == false && tidal == false)
	{
		construct_traceless_symmetric_matrix(b, b0_me);
		return 0;
	}

	/* calculate g */
	double g[9];
	calculate_g(g, G, m2, alpha_0, alpha, tilde_x, 
				b0_me, u_me, elements, bk_me,
				tidal);

	/* calculate c */
	double c = calculate_c(gamma, alpha_0, alpha);

	/* calculate b */
	scale_square_matrix(b, 1.0 / c, g);

	/* add centrifugal force if chosen */
	if (centrifugal == true)
	{
		/* calculate centrifugal force */
		double f_cent[9];
		calculate_f_cent(f_cent, omega);
		linear_combination_square_matrix(b,
			1.0, b,
			1.0 / c, f_cent);
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

	return 0;
}

int
calculate_inertia_tensor(double I[9], const double I0, const double b[9])
{
	double Id[9];
	identity_matrix(Id);
	linear_combination_square_matrix(I, I0, Id, -1.0 * I0, b);
	return 0;
}

int
calculate_l(double l[3], const double I0, 
	const double b[9], const double omega[3])
{
	double Id[9];
	identity_matrix(Id);
	double I[9];
	calculate_inertia_tensor(I, I0, b);
	square_matrix_times_vector(l, I, omega);
	return 0;
}

int
calculate_omega(double omega[3], const double omega_seed[3], const double G, 
	const double m2, const double I0, const double gamma, const double alpha_0, 
	const double alpha, const double tilde_x[3], const double l[3],
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[],
	const bool centrifugal, const bool tidal)
{
	/* for testing */
	// scale_vector(omega, 1.0 / I0, l);
	// return 0;

	/* method parameters */
	int 	number_of_iterates = 5;
	double 	max_error = 1e-8;
	double	error = 1.0;
	double 	previous_omega[3];
	copy_vector(omega, omega_seed);

	for (int i = 0; i < number_of_iterates; i++)
	{
		/* for testing */
		// print_vector(omega);

		/* store omega previous value */
		copy_vector(previous_omega, omega);

		/* calculate g */
		double g[9];
		calculate_g(g, G, m2, alpha_0, alpha, tilde_x, 
					b0_me, u_me, elements, bk_me,
					tidal);

		/* calculate c */
		double c = calculate_c(gamma, alpha_0, alpha);

		/* calculate H = 0 */
		double Id[9];
		identity_matrix(Id);
		double aux_H_first_term[9];
		linear_combination_square_matrix(aux_H_first_term,
			1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (centrifugal == true)
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
			-1.0 / I0, l);
		double minus_H[3];
		scale_vector(minus_H, -1.0, H);

		/* calculate DH */
		double DH[9];
		linear_combination_square_matrix(DH,
			1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (centrifugal == true)
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

	/* for testing */
	// printf("\nomega = \n");
	// print_vector(omega);
	// exit(42);

	return 0;
}

int
total_angular_momentum(double l_total[3],
	const double m1, const double m2,
	const double tilde_x[3], const double tilde_x_dot[3],
	const double l[3])
{
	double reduced_mass = (m1 * m2) / (m1 + m2);
	double x_cross_x_dot[3];
	cross_product(x_cross_x_dot, tilde_x, tilde_x_dot);
	double l_center_of_mass[3];
	scale_vector(l_center_of_mass, reduced_mass, x_cross_x_dot);

	linear_combination_vector(l_total, 1.0, l_center_of_mass, 1.0, l);

	return 0;
}

double
calculate_J2(const double m, const double R, const double I[9])
{
	double J2;

	double I_11 = I[0];
	double I_22 = I[4];
	double I_33 = I[8];

	J2 = (2.0 * I_33 - I_11 - I_22) / (2.0 * m * R * R);

	return J2;
}

double
calculate_C22(const double m, const double R, const double I[9])
{
	double C22;

	double I_11 = I[0];
	double I_22 = I[4];

	C22 = (I_22 - I_11) / (4.0 * m * R * R);

	return C22;
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