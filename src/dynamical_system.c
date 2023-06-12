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
		}			
	}

	/* preparing variables */

	double tilde_x[3], tilde_x_dot[3], l[3];
	double b0_me[9], u_me[9], *bk_me;

	for (int i = 0; i < 3; i++)	tilde_x[i] 		= y[0 + i];
	for (int i = 0; i < 3; i++) tilde_x_dot[i] 	= y[3 + i];
	for (int i = 0; i < 3; i++) l[i] 			= y[6 + i];
	for (int i = 0; i < 5; i++) b0_me[i] 		= y[9 + i];
	for (int i = 0; i < 5; i++) u_me[i] 		= y[14 + i];
	if (elements > 0)
	{
		bk_me = (double *) malloc(elements * 5 * sizeof(double));
		for (int i = 0; i < (elements * 5); i++)
		{
			bk_me[i] = y[19 + i];
		}
	}

	double b0[9], u[9];
	construct_traceless_symmetric_matrix(b0, b0_me);
	construct_traceless_symmetric_matrix(u, u_me);

	double *bk;
	if (elements > 0)
	{
		bk = (double *) malloc(elements * 9 * sizeof(double));
		construct_traceless_symmetric_matrix(bk, bk_me);
	}

	/* calculate omega and b */
	double omega[3], b[9]; 
	calculate_omega();	// tricky part
	calculate_b(b);

	/* useful definitions */

	double minus_G_times_total_mass = -1.0 * G * (m1 + m2);

	double tilde_x_norm 		= norm_vector(tilde_x);
	double tilde_x_norm_cube 	= pow(tilde_x_norm, 3.0);
	double tilde_x_norm_fifth 	= pow(tilde_x_norm, 5.0);
	double tilde_x_norm_seventh = pow(tilde_x_norm, 7.0);

	double bx[3];
	square_matrix_times_vector(bx, b, tilde_x);
	double bx_dot_x;
	bx_dot_x = dot_product(bx, tilde_x);

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

	/* writing components */	

	for (int i = 0; i < 3; i++) f[i] 		= component_tilde_x[i];
	for (int i = 0; i < 3; i++) f[3 + i] 	= component_tilde_x_dot[i];

	f[6]  = 0.0;
  	f[7]  = 0.0;
  	f[8]  = 0.0;
	f[9]  = 0.0;
  	f[10] = 0.0;
  	f[11] = 0.0;
	f[12] = 0.0;
  	f[13] = 0.0;
  	f[14] = 0.0;
  	f[15] = 0.0;
  	f[16] = 0.0;
  	f[17] = 0.0;
  	f[18] = 0.0;

	for (int i = 0; i < (elements * 5); i++) f[19 + i] = 0.0;
	
	/* freeing Voigt elements */
	if (elements > 0)
	{
		free(alpha_elements);
		free(eta_elements);
		free(bk_me);
		free(bk);
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
	copy_vector(b, null_M);
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
calculate_omega()
{
	return 0;
}