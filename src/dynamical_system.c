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

	gsl_matrix_view omega_hat_gsl	= gsl_matrix_view_array(omega_hat, 3, 3);
	gsl_matrix_view b_gsl			= gsl_matrix_view_array(b, 3, 3);
	
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

	gsl_vector_view tilde_x_gsl		= gsl_vector_view_array(tilde_x, 3);
	gsl_vector_view tilde_X_dot_gsl = gsl_vector_view_array(tilde_x_dot, 3);
	gsl_matrix_view l_hat_gsl 		= gsl_matrix_view_array(l_hat, 3, 3);
	gsl_matrix_view b0_gsl 			= gsl_matrix_view_array(b0, 3, 3);
	gsl_matrix_view u_gsl 			= gsl_matrix_view_array(u, 3, 3);

	double tilde_x_norm 		= gsl_blas_dnrm2 (&tilde_x_gsl.vector);
	double tilde_x_norm_cube 	= pow(tilde_x_norm, 3.0);
	double tilde_x_norm_fifth 	= pow(tilde_x_norm, 5.0);
	double tilde_x_norm_seventh = pow(tilde_x_norm, 7.0);

	/* preparing and writting components */

	// calculating tilde_x component

	double component_tilde_x[] = { tilde_x_dot[0], tilde_x_dot[1], tilde_x_dot[2] };

	for (int i = 0; i < 3; i++) f[i] = component_tilde_x[i];

	// calculating the 1st term of tilde_x_dot component

	double component_tilde_x_dot_1st_term[] = { 0.0, 0.0, 0.0};

	gsl_vector_view component_tilde_x_dot_1st_term_gsl 
		= gsl_vector_view_array(component_tilde_x_dot_1st_term, 3);

	gsl_blas_daxpy ((1.0 / tilde_x_norm_cube), &tilde_x_gsl.vector, 
		&component_tilde_x_dot_1st_term_gsl.vector);

	// calculating the 2nd term of tilde_x_dot component

	double component_tilde_x_dot_2nd_term[] = { 0.0, 0.0, 0.0};

	gsl_vector_view component_tilde_x_dot_2nd_term_gsl 
		= gsl_vector_view_array(component_tilde_x_dot_2nd_term, 3);
		
	gsl_blas_dgemv(CblasNoTrans, 1.0, &b_gsl.matrix, &tilde_x_gsl.vector, 
		0.0, &component_tilde_x_dot_2nd_term_gsl.vector);

	double result_dot_product;
	gsl_blas_ddot(&component_tilde_x_dot_2nd_term_gsl.vector, 
		&tilde_x_gsl.vector, &result_dot_product);
	

	
	// calculating the 3rd term of tilde_x_dot component

	double component_tilde_x_dot_3rd_term[] = { 0.0, 0.0, 0.0};

	gsl_vector_view component_tilde_x_dot_3rd_term_gsl 
		= gsl_vector_view_array(component_tilde_x_dot_3rd_term, 3);

	// calculating tilde_x_dot component

	double component_tilde_x_dot[] = { 0.0, 0.0, 0.0};

	gsl_vector_view component_tilde_x_dot_gsl 
		= gsl_vector_view_array(component_tilde_x_dot, 3);

	gsl_blas_daxpy (minus_G_times_total_mass, &component_tilde_x_dot_1st_term_gsl.vector, 
		&component_tilde_x_dot_gsl.vector);
	gsl_blas_daxpy (minus_G_times_total_mass, &component_tilde_x_dot_2nd_term_gsl.vector, 
		&component_tilde_x_dot_gsl.vector);
	gsl_blas_daxpy (minus_G_times_total_mass, &component_tilde_x_dot_3rd_term_gsl.vector, 
		&component_tilde_x_dot_gsl.vector);

	for (int i = 0; i < 3; i++) f[3 + i] = component_tilde_x_dot[i];

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