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
	double 	b0_diag[3];
	/* non-state variables given by user */
	double 	omega_vec[3];
	/* system parameters given by user */
	int		elements = 0; //number of Voigt elements
	double  G, m1, m2, I0;
	double	alpha_0, alpha, eta;

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
		else if (strcmp(var_name, "b0_x") == 0)			b0_diag[0] = var_value;
		else if (strcmp(var_name, "b0_y") == 0)			b0_diag[1] = var_value;
		else if (strcmp(var_name, "b0_z") == 0)			b0_diag[2] = var_value;
		else if (strcmp(var_name, "omega_x") == 0)		omega_vec[0] = var_value;
		else if (strcmp(var_name, "omega_y") == 0)		omega_vec[1] = var_value;
		else if (strcmp(var_name, "omega_z") == 0)		omega_vec[2] = var_value;
		else if (strcmp(var_name, "alpha_0") == 0)		alpha_0 = var_value;
		else if (strcmp(var_name, "alpha") == 0)		alpha = var_value;
		else if (strcmp(var_name, "eta") == 0)			eta = var_value;
		else if (strcmp(var_name, "elements") == 0)		elements = (int) var_value;
	}
	fclose(in1_again);

	/* setting up Voigt elements */
	double	*alpha_elements, *eta_elements;
	if (elements > 0)
	{
		alpha_elements 	= (double *) malloc(elements * sizeof(double));
		eta_elements 	= (double *) malloc(elements * sizeof(double));

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

	/* position and velocity at periapsis given by Murray */
	double	tilde_x[3], tilde_x_dot[3];
	double	n = (2.0 * PI) / T;
	tilde_x[0] 		= a * (1.0 - e);
	tilde_x[1] 		= 0.0;
	tilde_x[2] 		= 0.0;	
	tilde_x_dot[0]	= 0.0;
	tilde_x_dot[1] 	= n * a * sqrt((1.0 + e)/(1.0 - e));
	tilde_x_dot[2] 	= 0.0;

	/* complete b0_me */
	double b0_me[5];
	b0_me[0] = b0_diag[0];
	b0_me[1] = 0.0;
	b0_me[2] = 0.0;
	b0_me[3] = b0_diag[1];
	b0_me[4] = 0.0;

	/* variables and parameters passed as field parameters */
	int		dim_params = 7 + (elements * 2); // optmize this
	double	params[dim_params]; 
	params[0] = G;
	params[1] = m1;
	params[2] = m2;
	params[3] = I0;
	params[4] = alpha;
	params[5] = eta;
	params[6] = elements;
	for (int i = 0; i < elements; i++)
	{
		params[7 + (2*i)] = alpha_elements[i];
		params[8 + (2*i)] = eta_elements[i];
	}

	/* variables not given by user */
	double u_me[5], *bk_me;
	for (int i = 0; i < 5; i++) u_me[i] = 0.0;
	if (elements > 0)
	{
		bk_me = (double *) calloc(elements * 5, sizeof(double));
	}

	/* calculate first l */
	double 	b[9], l[3];
	calculate_b(b);
	calculate_l(l, I0, b, omega_vec);

	/* integration loop variables */
	int		dim = 19 + (elements * 5);	// optmize this
	double 	t = t0;
	double 	y[dim];
	for (int i = 0; i < 3; i++) 				y[0 + i] = tilde_x[i];
	for (int i = 0; i < 3; i++) 				y[3 + i] = tilde_x_dot[i];
	for (int i = 0; i < 3; i++) 				y[6 + i] = l[i];
	for (int i = 0; i < 5; i++) 				y[9 + i] = b0_me[i];
	for (int i = 0; i < 5; i++) 				y[14 + i] = u_me[i];
	for (int i = 0; i < (elements * 5); i++) 	y[19 + i] = bk_me[i];

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

		counter++;
	}

	/* free Voigt elements */
	if (elements > 0)
	{
		free(alpha_elements);
		free(eta_elements);	
		free(bk_me);
	}

	/* free GSL variables */
	gsl_odeiv2_evolve_free (ode_evolve);
	gsl_odeiv2_control_free (ode_control);
	gsl_odeiv2_step_free (ode_step);

	return 0;
}
