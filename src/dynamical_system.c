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
	double 	alpha 	 = par[4];
	double 	eta 	 = par[5];
	int 	elements = par[6];
	double	*alpha_elements, *eta_elements;
	if (elements > 0)
	{
		alpha_elements 	= (double *) malloc(elements * sizeof(double));
		eta_elements 	= (double *) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			alpha_elements[i] 	= par[7 + (2*i)];
			eta_elements[i]		= par[8 + (2*i)];

			/* for testing */
			// printf("alpha_%d = %f\n", i, alpha_elements[i]);
			// printf("eta_%d = %f\n", i, eta_elements[i]);
		}
	}

	/* preparing variables */

	double tilde_x[3], tilde_x_dot[3], l[3];
	double b0_me[9], u_me[9], **bk_me;

	for (int i = 0; i < 3; i++)	tilde_x[i] 		= y[0 + i];
	for (int i = 0; i < 3; i++) tilde_x_dot[i] 	= y[3 + i];
	for (int i = 0; i < 3; i++) l[i] 			= y[6 + i];
	for (int i = 0; i < 5; i++) b0_me[i] 		= y[9 + i];
	for (int i = 0; i < 5; i++) u_me[i] 		= y[14 + i];
	if (elements > 0)
	{
		bk_me = (double **) malloc(elements * sizeof(double));
		for (int i = 0; i < elements; i++)
		{
			bk_me[i] = (double *) malloc(5 * sizeof(double));
			for (int j = 0; j < 5; j++)
			{
				bk_me[i][j] = y[19 + j + (i*5)];

				/* for testing */
				// printf("bk_me = %f\n",bk_me[i][j]);
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
			construct_traceless_symmetric_matrix(bk[i], bk_me[i]);

			/* for testing */
			// print_square_matrix(bk[i]);
		}
	}

	/* calculate omega and b */
	double omega[3], b[9]; 
	calculate_omega(omega);	// tricky part
	calculate_b(b);

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
		for (int i = 0; i < elements; i++) free(bk_me[i]);
		free(bk_me);
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

int
calculate_b(double b[9])
{
	double null_M[9];
	null_matrix(null_M);
	copy_square_matrix(b, null_M);
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