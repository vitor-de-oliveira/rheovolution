#include "dynamical_system.h"

int
field_1EB1PM(double t, const double y[], double f[],
       		 void *params)
{
	(void)(t);

	/* preparing parameters */

	double 	*par = (double *)params;

	double 	G		 = par[0];
	double 	m1		 = par[1];
	double 	m2		 = par[2];
	double 	I0 		 = par[3];
	double	gamma	 = par[4];
	double 	alpha 	 = par[5];
	double 	eta 	 = par[6];
	double	alpha_0	 = par[7];
	int 	elements = par[8];
	double	*alpha_elements, *eta_elements;
	if (elements > 0)
	{
		alpha_elements 	= (double *) malloc(elements * sizeof(double));
		eta_elements 	= (double *) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			alpha_elements[i] 	= par[9 + (2*i)];
			eta_elements[i]		= par[10 + (2*i)];

			/* for testing */
			// printf("alpha_%d = %f\n", i, alpha_elements[i]);
			// printf("eta_%d = %f\n", i, eta_elements[i]);
		}
	}

	/* preparing variables */

	double tilde_x[3], tilde_x_dot[3], l[3];
	double b0_me[9], u_me[9], *bk_me, **bk_me_2d_array;

	for (int i = 0; i < 3; i++)	tilde_x[i] 		= y[0 + i];
	for (int i = 0; i < 3; i++) tilde_x_dot[i] 	= y[3 + i];
	for (int i = 0; i < 3; i++) l[i] 			= y[6 + i];
	for (int i = 0; i < 5; i++) b0_me[i] 		= y[9 + i];
	for (int i = 0; i < 5; i++) u_me[i] 		= y[14 + i];
	if (elements > 0)
	{
		bk_me = (double *) malloc(elements * 5 * sizeof(double));
		bk_me_2d_array = (double **) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			bk_me[i] = y[19 + i];
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

	double **bk;
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
	calculate_omega(omega);	// tricky part
	null_matrix(b);
	// calculate_b(b, G, m2, gamma, alpha_0, alpha,
	// 	tilde_x, omega, b0_me, u_me, elements, bk_me);

	double omega_hat[9];
	hat_map(omega_hat, omega);

	/* useful definitions */

	double minus_G_times_total_mass = -1.0 * G * (m1 + m2);

	double tilde_x_norm 		= norm_vector(tilde_x);
	double tilde_x_norm_cube 	= pow(tilde_x_norm, 3.0);
	double tilde_x_norm_fifth 	= pow(tilde_x_norm, 5.0);
	double tilde_x_norm_seventh = pow(tilde_x_norm, 7.0);

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
	double bx_dot_x;
	bx_dot_x = dot_product(bx, tilde_x);
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
		minus_G_times_total_mass / tilde_x_norm_cube, tilde_x);

	double component_tilde_x_dot_2nd_term[] = { 0.0, 0.0, 0.0 };
	scale_vector (component_tilde_x_dot_2nd_term, 
		(15. * I0 * bx_dot_x) / (2. * m1 * tilde_x_norm_seventh), tilde_x);

	double component_tilde_x_dot_3rd_term[] = { 0.0, 0.0, 0.0 };
	scale_vector (component_tilde_x_dot_3rd_term, 
		(-3.0 * I0) / (m1 * tilde_x_norm_fifth), bx);

	double component_tilde_x_dot[] = { 0.0, 0.0, 0.0};
	linear_combination_three_vector(component_tilde_x_dot,
		1.0, component_tilde_x_dot_1st_term, 
		1.0, component_tilde_x_dot_2nd_term, 
		1.0, component_tilde_x_dot_3rd_term);

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
	linear_combination_square_matrix(component_u, 1.0, omega_hat_comm_u, -1.0 / tau, lambda);
	double component_u_me[] = { 0.0, 0.0, 0.0,
							         0.0, 0.0 };
	get_main_elements_traceless_symmetric_matrix(component_u_me, component_u);

	// bk components

	double **component_bk_me;
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

			double component_bk[] = { 0.0, 0.0, 0.0,
							 		  0.0, 0.0, 0.0,
							          0.0, 0.0, 0.0 };
			linear_combination_three_square_matrix(component_bk,
				1.0, omega_hat_comm_bk,
				1.0, minus_bk_over_tau_elements,
				1.0, lambda_over_eta_elements);

			get_main_elements_traceless_symmetric_matrix(component_bk_me[i],
				component_bk);
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
parameter_gamma_homogeneous_body(const double G,
	const double I0, const double R)
{
	return 2.0 * I0 * G / pow(R, 5.0);
}

int
calculate_b(double b[9], const double G, const double m2, 
	const double gamma, const double alpha_0, const double alpha,
	const double tilde_x[3], const double omega[3], 
	const double b0_me[5], const double u_me[5],
	const int elements, const double bk_me[])
{
	double b0[9], u[9];
	construct_traceless_symmetric_matrix(b0, b0_me);
	construct_traceless_symmetric_matrix(u, u_me);
	double omega_hat[9];
	hat_map(omega_hat, omega);
	
	double x_tensor_x[9];
	tensor_product(x_tensor_x, tilde_x, tilde_x);
	double Id[9];
	identity_matrix(Id);
	double scaled_id[9];
	scale_square_matrix(scaled_id, 
		norm_squared_vector(tilde_x) / 3.0, Id);
	double tilde_x_norm_fifth = pow(norm_vector(tilde_x), 5.0);
	double f_tide[9];
	linear_combination_square_matrix(f_tide, 
		 3.0 * G * m2 / tilde_x_norm_fifth, x_tensor_x,
		-3.0 * G * m2 / tilde_x_norm_fifth, scaled_id);

	double alpha_0_b0[9];
	scale_square_matrix(alpha_0_b0, alpha_0, b0);
	double g[9];
	linear_combination_three_square_matrix(g,
		 1.0, f_tide,
		 1.0, alpha_0_b0,
		-1.0, u);

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
	
	double c = gamma + alpha_0 + alpha;

	double omega_hat_squared[9];
	square_matrix_times_square_matrix(omega_hat_squared,
		omega_hat, omega_hat);
	double trace_omega_hat_squared;
	trace_omega_hat_squared = trace_square_matrix(omega_hat_squared);
	double scaled_id_2[9];
	scale_square_matrix(scaled_id_2, 
		trace_omega_hat_squared / 3.0, Id);

	linear_combination_three_square_matrix(b,
		-1.0 / c, omega_hat_squared,
		 1.0 / c, scaled_id_2,
		 1.0 / c, g);

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
calculate_l(double l[3], const double I0, 
	const double b[9], const double omega[3])
{
	double Id[9];
	identity_matrix(Id);
	double I[9];
	linear_combination_square_matrix(I, I0, Id, -1.0 * I0, b);
	square_matrix_times_vector(l, I, omega);
	return 0;
}

int
calculate_omega(double omega[3])
{
	double null_v[3];
	null_vector(null_v);
	copy_vector(omega, null_v);
	return 0;
}