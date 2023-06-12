#include "dynamical_system.h"

int
field_1EB1PM(double t, const double y[], double f[],
       		 void *params)
{
	(void)(t);

	/* preparing parameters */

	double 	*par = (double *)params;
	double	*alpha_elements, *eta_elements, *bk;
	double	omega_hat[9], b[9];

	double 	G		 = par[0];
	double 	m1		 = par[1];
	double 	m2		 = par[2];
	double 	I0 		 = par[3];
	double 	alpha 	 = par[4];
	double 	eta 	 = par[5];
	int 	elements = par[6];

	if (elements > 0)
	{
		alpha_elements 	= (double *) malloc(elements * sizeof(double));
		eta_elements 	= (double *) malloc(elements * sizeof(double));
		bk				= (double *) calloc(elements * 9, sizeof(double));
	}

	for (int i = 0; i < elements; i++) 	alpha_elements[i] = par[7 + (2*i)];
	for (int i = 0; i < elements; i++) 	eta_elements[i] = par[8 + (2*i)];
	for (int i = 0; i < 9; i++) 		omega_hat[i] = par[7 + (2*elements) + i];
	for (int i = 0; i < 9; i++)			b[i] = par[16 + (2*elements) + i];
	
  	double minus_G_times_total_mass = -1.0 * G * (m1 + m2);

	/* preparing variables */

	double tilde_x[] 		= { y[0], y[1], y[2] };
	double tilde_x_dot[] 	= { y[3], y[4], y[5] };
	double l_hat[] 			= { y[6], y[7], y[8],
								y[9], y[10], y[11],
								y[12], y[13], y[14] };
	double b0[]				= { y[15], y[16], y[17],
								y[18], y[19], y[20],
								y[21], y[22], y[23] };
	double u[]				= { y[24], y[25], y[26],
								y[27], y[28], y[29],
								y[30], y[31], y[32] };
	for (int i = 0; i < (elements * 9); i++)
	{
		bk[i] = y[33 + i];
	}

	double tilde_x_norm 		= norm_vector (tilde_x);
	double tilde_x_norm_cube 	= pow(tilde_x_norm, 3.0);
	double tilde_x_norm_fifth 	= pow(tilde_x_norm, 5.0);
	double tilde_x_norm_seventh = pow(tilde_x_norm, 7.0);

	/* calculating components */

	// tilde_x component

	double component_tilde_x[] = { tilde_x_dot[0], tilde_x_dot[1], tilde_x_dot[2] };

	// 1st term of tilde_x_dot component

	double component_tilde_x_dot_1st_term[] = { 0.0, 0.0, 0.0};
	scale_vector (component_tilde_x_dot_1st_term, 
		minus_G_times_total_mass / tilde_x_norm_cube, tilde_x);

	// 2nd term of tilde_x_dot component

	double component_tilde_x_dot_2nd_term[] = { 0.0, 0.0, 0.0};


	// 3rd term of tilde_x_dot component

	double component_tilde_x_dot_3rd_term[] = { 0.0, 0.0, 0.0};


	// tilde_x_dot component

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
  	f[19] = 0.0;
  	f[20] = 0.0;
  	f[21] = 0.0;
  	f[22] = 0.0;
  	f[23] = 0.0;
  	f[24] = 0.0;
  	f[25] = 0.0;
  	f[26] = 0.0;
  	f[27] = 0.0;
  	f[28] = 0.0;
  	f[29] = 0.0;
  	f[30] = 0.0;
  	f[31] = 0.0;
  	f[32] = 0.0;

	// calculating bk components

	for (int i = 0; i < (elements * 9); i++) f[33 + i] = bk[i];
	
	/* freeing Voigt elements */

	if (elements > 0)
	{
		free(alpha_elements);
		free(eta_elements);
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

double
parameter_gamma()
{
	double gamma = 0.0;
	return gamma;
}

int
calculate_f  (double f[9], double omega_hat[9],
              double x[3], double G, double m2)
{
	return 0;
}

int
calculate_lambda(double lambda[9], double f[9],
                 double b0_matrix[9], double u[9],
                 double gamma, double alpha_0,
                 double alpha)
{
	return 0;
}

int
calculate_b (double b[9], double f[9], double lambda[9],
             double b0_matrix[9],
             double gamma, double alpha_0)
{
	return 0;
}

int
calculate_l_hat (double l_hat[9], double omega_hat[9], double b[9],
		     	 double I0)
{
	return 0;
}

int
calculate_omega_hat (double omega_hat[9], double b[9], double l_hat[9],
	             	 double I0)
{
	return 0;
}