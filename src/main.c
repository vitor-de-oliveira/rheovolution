#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "dynamical_system.h"

int
main(int argc, char *argv[]) 
{
	/* check number of input files */
	if( argc < 3 )
	{
	   	printf("Please, provide both the system parameters");
		printf(" and the integrator parameters files\n");
	   	exit(2);
   	}
   	else if ( argc > 3 )
   	{
	   	printf("Too many arguments\n");
	   	exit(4);
   	}

	/* make a copy of the input files */
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}
	FILE *in1 = fopen(argv[1], "r");
	FILE *in1_copy = fopen("output/input_pars_system_copy.txt" , "w");
	char ch = fgetc(in1);
    while(ch != EOF)
    {
        fputc(ch, in1_copy);
        ch = fgetc(in1);
    }
	fclose(in1);
	fclose(in1_copy);
	FILE *in2 = fopen(argv[2], "r");
	FILE *in2_copy = fopen("output/input_pars_integrator_copy.txt" , "w");
	char ch2 = fgetc(in2);
    while(ch2 != EOF)
    {
        fputc(ch2, in2_copy);
        ch2 = fgetc(in2);
    }
	fclose(in2);
	fclose(in2_copy);

	/* auxiliary variables for fscanf */
	char var_name[100];
	double var_value;

	/* state variables given by user*/
	double x[3], x_dot[3], b0[3];
	/* non-state variables given by user */
	double omega[3];
	/* system parameters given by user */
	double I0, alpha_0, alpha, eta;
	double G, m1, m2;

	/* reading system specs from user */
	FILE *in1_again = fopen(argv[1], "r");
	while(fscanf(in1_again, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "x") == 0)					x[0] = var_value;
		else if (strcmp(var_name, "y") == 0)			x[1] = var_value;
		else if (strcmp(var_name, "z") == 0)			x[2] = var_value;
		else if (strcmp(var_name, "x_dot") == 0)		x_dot[0] = var_value;
		else if (strcmp(var_name, "y_dot") == 0)		x_dot[1] = var_value;
		else if (strcmp(var_name, "z_dot") == 0)		x_dot[2] = var_value;
		else if (strcmp(var_name, "b0_x") == 0)			b0[0] = var_value;
		else if (strcmp(var_name, "b0_y") == 0)			b0[1] = var_value;
		else if (strcmp(var_name, "b0_z") == 0)			b0[2] = var_value;
		else if (strcmp(var_name, "omega_x") == 0)		omega[0] = var_value;
		else if (strcmp(var_name, "omega_y") == 0)		omega[1] = var_value;
		else if (strcmp(var_name, "omega_z") == 0)		omega[2] = var_value;
		else if (strcmp(var_name, "I0") == 0)			I0 = var_value;
		else if (strcmp(var_name, "alpha_0") == 0)		alpha_0 = var_value;
		else if (strcmp(var_name, "alpha") == 0)		alpha = var_value;
		else if (strcmp(var_name, "eta") == 0)			eta = var_value;
		else if (strcmp(var_name, "G") == 0)			G = var_value;
		else if (strcmp(var_name, "m1") == 0)			m1 = var_value;
		else if (strcmp(var_name, "m2") == 0)			m2 = var_value;
	}
	fclose(in1_again);

	/* integrator parameters */
	double 	h, t0, t1, eps_abs, eps_rel;
	int		data_step;

	/* reading integrator specs from user */
	FILE *in2_again = fopen(argv[2], "r");
	while(fscanf(in2_again, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "h") == 0)					h = var_value;
		else if (strcmp(var_name, "t0") == 0)			t0 = var_value;
		else if (strcmp(var_name, "t1") == 0)			t1 = var_value;
		else if (strcmp(var_name, "eps_abs") == 0)		eps_abs = var_value;
		else if (strcmp(var_name, "eps_rel") == 0)		eps_rel = var_value;
		else if (strcmp(var_name, "data_step") == 0)	data_step = (int) var_value;
	}
	fclose(in2_again);

	double	omega_hat[9];
	hat_map(omega_hat, omega);

	double	b0_matrix[9];
	for (int i = 0; i < 9; i++) b0_matrix[i] = 0.0;
	b0_matrix[0] = b0[0];
	b0_matrix[3] = b0[1];
	b0_matrix[6] = b0[2];

	/* state variables not given by user */
	double l_hat[9], u[9];
	/* non-state variables not given by user */
	double b[9], f[9], lambda[9];
	/* system parameters not given by user */
	double gamma = parameter_gamma();

	/* calculation of non-state variables values */
	calculate_f(f, omega_hat, x, G, m2);
	calculate_lambda(lambda, f, b0_matrix, u, gamma, alpha_0, alpha);
	calculate_b(b, f, lambda, b0_matrix, gamma, alpha_0);

	/* initialization of remainder state variables values */
	calculate_l_hat(l_hat, omega_hat, b, I0);
	for (int i = 0; i < 9; i++) u[i] = 0.0;

	/* variables and parameters passed as field parameters */
	int		dim_params = 24; // optmize this
	double	params[dim_params]; 
	params[0] = I0;
	params[1] = alpha;
	params[2] = eta;
	params[3] = G;
	params[4] = m1;
	params[5] = m2;
	for (int i = 0; i < 9; i++) params[6 + i] = b[i];
	for (int i = 0; i < 9; i++) params[15 + i] = omega_hat[i];

	/* integration loop variables */
	int		dim = 33;	// optmize this
	double 	t = t0;
	double 	y[dim];
	y[0]  = x[0];
	y[1]  = x[1];
	y[2]  = x[2];
	y[3]  = x_dot[0];
	y[4]  = x_dot[1];
	y[5]  = x_dot[2];
	y[6]  = l_hat[0];
  	y[7]  = l_hat[1];
  	y[8]  = l_hat[2];
	y[9]  = l_hat[3];
  	y[10] = l_hat[4];
  	y[11] = l_hat[5];
	y[12] = l_hat[6];
  	y[13] = l_hat[7];
  	y[14] = l_hat[8];
  	y[15] = b0_matrix[0];
  	y[16] = b0_matrix[1];
  	y[17] = b0_matrix[2];
  	y[18] = b0_matrix[3];
  	y[19] = b0_matrix[4];
  	y[20] = b0_matrix[5];
  	y[21] = b0_matrix[6];
  	y[22] = b0_matrix[7];
  	y[23] = b0_matrix[8];
  	y[24] = u[0];
  	y[25] = u[1];
  	y[26] = u[2];
  	y[27] = u[3];
  	y[28] = u[4];
  	y[29] = u[5];
  	y[30] = u[6];
  	y[31] = u[7];
  	y[32] = u[8];

	/* GSL variables */
	const gsl_odeiv2_step_type * T
		= gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * s
    	= gsl_odeiv2_step_alloc (T, dim);
  	gsl_odeiv2_control * c
    	= gsl_odeiv2_control_y_new (eps_abs, eps_rel);
  	gsl_odeiv2_evolve * e
    	= gsl_odeiv2_evolve_alloc (dim);
  	gsl_odeiv2_system sys = {field, NULL, dim, params};

	/* integration loop */
	int counter = 0;
	while (t < t1)
	{
		int status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &t, h, y);
		
		if (status != GSL_SUCCESS) break;
		
		if (counter % data_step == 0)
		{
			printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
		}

		/* update b and omega_hat */
		for (int i = 0; i < 1; i++)	// can be modified to add more loops
		{
			calculate_omega_hat(omega_hat, b, l_hat, I0);
			calculate_f(f, omega_hat, x, G, m2);
			calculate_lambda(lambda, f, b0_matrix, u, gamma, alpha_0, alpha);
			calculate_b(b, f, lambda, b0_matrix, gamma, alpha_0);
		}
		calculate_omega_hat(omega_hat, b, l_hat, I0);

		/* update params array */
		for (int i = 0; i < 9; i++) params[6 + i] = b[i];
		for (int i = 0; i < 9; i++) params[15 + i] = omega_hat[i];
		
		/* update sys */
		sys.params = params;

		counter++;
	}

	/* free GSL variables */
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	return 0;
}
