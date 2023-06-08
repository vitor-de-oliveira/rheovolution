#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "dynamical_system.h"

#define PI 3.14159265358979323846

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
	char 	var_name[100];
	double 	var_value;

	/* orbital parameters given by user */
	double 	e, a, T;
	/* state variables given by user*/
	double 	b0[3];
	/* non-state variables given by user */
	double 	omega[3];
	/* system parameters given by user */
	int		elements = 0; //number of Voigt elements
	double  I0, alpha_0, alpha, eta;
	double 	G, m1, m2;

	/* reading system specs from user */
	FILE *in1_again = fopen(argv[1], "r");
	while(fscanf(in1_again, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "e") == 0)					e = var_value;
		else if (strcmp(var_name, "a") == 0)			a = var_value;
		else if (strcmp(var_name, "T") == 0)			T = var_value;
		else if (strcmp(var_name, "G") == 0)			G = var_value;
		else if (strcmp(var_name, "m1") == 0)			m1 = var_value;
		else if (strcmp(var_name, "m2") == 0)			m2 = var_value;
		else if (strcmp(var_name, "I0") == 0)			I0 = var_value;
		else if (strcmp(var_name, "b0_x") == 0)			b0[0] = var_value;
		else if (strcmp(var_name, "b0_y") == 0)			b0[1] = var_value;
		else if (strcmp(var_name, "b0_z") == 0)			b0[2] = var_value;
		else if (strcmp(var_name, "omega_x") == 0)		omega[0] = var_value;
		else if (strcmp(var_name, "omega_y") == 0)		omega[1] = var_value;
		else if (strcmp(var_name, "omega_z") == 0)		omega[2] = var_value;
		else if (strcmp(var_name, "alpha_0") == 0)		alpha_0 = var_value;
		else if (strcmp(var_name, "alpha") == 0)		alpha = var_value;
		else if (strcmp(var_name, "eta") == 0)			eta = var_value;
		else if (strcmp(var_name, "elements") == 0)		elements = (int) var_value;
	}
	fclose(in1_again);

	/* setting Voigt elements */
	double	*alpha_elements, *eta_elements, *bk;
	if (elements > 0)
	{
		alpha_elements 	= (double *) malloc(elements * sizeof(double));
		eta_elements 	= (double *) malloc(elements * sizeof(double));
		bk				= (double *) calloc(elements * 9, sizeof(double));

		char name_element_alpha[20], name_element_eta[20];
		bool element_alpha_found[elements];
		bool element_eta_found[elements];
		for (int i = 0; i < elements; i++)
		{
			element_alpha_found[i] = false;
			element_eta_found[i] = false;
		}

		FILE *in1_elements = fopen(argv[1], "r");
		while(fscanf(in1_again, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
		{
			for (int i = 0; i < elements; i++)
			{
				sprintf(name_element_alpha, "alpha_%d", i+1);
				sprintf(name_element_eta, "eta_%d", i+1);
				if (strcmp(var_name, name_element_alpha) == 0)
				{
					alpha_elements[i] = var_value;
					element_alpha_found[i] = true;
				}
				else if (strcmp(var_name, name_element_eta) == 0)
				{
					eta_elements[i] = var_value;
					element_eta_found[i] = true;
				}
			}
		}
		fclose(in1_elements);
		for (int i = 0; i < elements; i++)
		{
			if (element_alpha_found[i] == false || 
				element_eta_found[i] == false)
			{
				printf("Error: parameters missing for Voigt elements.\n");
				exit(10);
			}
		}
	}

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

	/* position and velocity at periapsis by Murray */
	double	tilde_x[3], tilde_x_dot[3];
	double	n = (2.0 * PI) / T;
	tilde_x[0] = a * (1.0 - e);
	tilde_x[1] = 0.0;
	tilde_x[2] = 0.0;	
	tilde_x_dot[0] = 0.0;
	tilde_x_dot[1] = n * a * sqrt((1.0 + e)/(1.0 - e));
	tilde_x_dot[2] = 0.0;

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
	calculate_f(f, omega_hat, tilde_x, G, m2);
	calculate_lambda(lambda, f, b0_matrix, u, gamma, alpha_0, alpha);
	calculate_b(b, f, lambda, b0_matrix, gamma, alpha_0);

	/* initialization of remainder state variables values */
	calculate_l_hat(l_hat, omega_hat, b, I0);
	for (int i = 0; i < 9; i++) u[i] = 0.0;

	/* variables and parameters passed as field parameters */
	int		dim_params = 25 + (elements * 2); // optmize this
	double	params[dim_params]; 
	params[0] = G;
	params[1] = m1;
	params[2] = m2;
	params[3] = I0;
	params[4] = alpha;
	params[5] = eta;
	params[6] = elements;
	for (int i = 0; i < elements; i++) 	params[7 + (2*i)] = alpha_elements[i];
	for (int i = 0; i < elements; i++) 	params[8 + (2*i)] = eta_elements[i];
	for (int i = 0; i < 9; i++) 		params[7 + (2*elements) + i] = omega_hat[i];
	for (int i = 0; i < 9; i++)			params[16 + (2*elements) + i] = b[i];

	/* integration loop variables */
	int		dim = 33 + (elements * 9);	// optmize this
	double 	t = t0;
	double 	y[dim];
	for (int i = 0; i < 3; i++) 				y[0 + i] = tilde_x[i];
	for (int i = 0; i < 3; i++) 				y[3 + i] = tilde_x_dot[i];
	for (int i = 0; i < 9; i++) 				y[6 + i] = l_hat[i];
	for (int i = 0; i < 9; i++) 				y[15 + i] = b0_matrix[i];
	for (int i = 0; i < 9; i++) 				y[24 + i] = u[i];
	for (int i = 0; i < (elements * 9); i++) 	y[33 + i] = bk[i];

	/* GSL variables */
	const gsl_odeiv2_step_type * ode_type
		= gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * ode_step
    	= gsl_odeiv2_step_alloc (ode_type, dim);
  	gsl_odeiv2_control * ode_control
    	= gsl_odeiv2_control_y_new (eps_abs, eps_rel);
  	gsl_odeiv2_evolve * ode_evolve
    	= gsl_odeiv2_evolve_alloc (dim);
  	gsl_odeiv2_system sys = {field_1EB1PM, NULL, dim, params};

	/* integration loop */
	int counter = 0;
	while (t < t1)
	{
		int status = 
			gsl_odeiv2_evolve_apply_fixed_step (ode_evolve, 
				ode_control, ode_step, &sys, &t, h, y);
		
		if (status != GSL_SUCCESS) break;
		
		if (counter % data_step == 0)
		{
			printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
		}

		/* update b and omega_hat */
		for (int i = 0; i < 1; i++)	// can be modified to add more loops
		{
			calculate_omega_hat(omega_hat, b, l_hat, I0);
			calculate_f(f, omega_hat, tilde_x, G, m2);
			calculate_lambda(lambda, f, b0_matrix, u, gamma, alpha_0, alpha);
			calculate_b(b, f, lambda, b0_matrix, gamma, alpha_0);
		}
		calculate_omega_hat(omega_hat, b, l_hat, I0);

		/* update params array */
		for (int i = 0; i < 9; i++) params[7 + (2*elements) + i]  = omega_hat[i];
		for (int i = 0; i < 9; i++) params[16 + (2*elements) + i] = b[i];
		
		/* update sys */
		sys.params = params;

		counter++;
	}

	/* free Voigt elements */
	if (elements > 0)
	{
		free(alpha_elements);
		free(eta_elements);	
		free(bk);
	}

	/* free GSL variables */
	gsl_odeiv2_evolve_free (ode_evolve);
	gsl_odeiv2_control_free (ode_control);
	gsl_odeiv2_step_free (ode_step);

	return 0;
}
