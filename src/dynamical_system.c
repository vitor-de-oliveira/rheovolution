#include "dynamical_system.h"

int
field (double t, const double y[], double f[],
       void *params)
{
	(void)(t);

	double *par = (double *)params;

	double I0 		= par[0];
	double alpha 	= par[1];
	double eta 		= par[2];
	double G		= par[3];
	double m1		= par[4];
	double m2		= par[5];
	
  	double total_mass = m1 + m2;
  	
  	double r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
  	double r_cube = r * r * r;
  
  	/**
  	 * y[0]  = x
  	 * y[1]  = y
  	 * y[2]  = z
  	 * y[3]  = x_dot
  	 * y[4]  = y_dot
  	 * y[5]  = z_dot
  	 * y[6]  = l_11
  	 * y[7]  = l_12
  	 * y[8]  = l_13
  	 * y[9]  = l_21
  	 * y[10] = l_22
  	 * y[11] = l_23
  	 * y[12] = l_31
  	 * y[13] = l_32
  	 * y[14] = l_33
  	 * y[15] = b_0_11
  	 * y[16] = b_0_12
  	 * y[17] = b_0_13
  	 * y[18] = b_0_21
  	 * y[19] = b_0_22
  	 * y[20] = b_0_23
  	 * y[21] = b_0_31
  	 * y[22] = b_0_32
  	 * y[23] = b_0_33
  	 * y[24] = u_11
  	 * y[25] = u_12
  	 * y[26] = u_13
  	 * y[27] = u_21
  	 * y[28] = u_22
  	 * y[29] = u_23
  	 * y[30] = u_31
  	 * y[31] = u_32
  	 * y[32] = u_33
	*/
  
	f[0]  = y[3];
	f[1]  = y[4];
	f[2]  = y[5];
	f[3]  = -1.0 * G * total_mass * y[0] / r_cube;
	f[4]  = -1.0 * G * total_mass * y[1] / r_cube;
	f[5]  = -1.0 * G * total_mass * y[2] / r_cube;

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