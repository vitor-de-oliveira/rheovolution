#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "dynamical_system.h"
#include "convert.h"

#define PI 3.14159265358979323846

int
main(int argc, char *argv[]) 
{
	/* check number of input files */
	if( argc < 3 )
	{
	   	fprintf(stderr, "Please, provide both the system parameters");
		fprintf(stderr, " and the integrator parameters files\n");
	   	exit(2);
   	}
   	else if ( argc > 3 )
   	{
	   	fprintf(stderr, "Too many arguments\n");
	   	exit(4);
   	}

	/* convert input */
	convert_input(argv[1]);
	exit(99);

	/* auxiliary variables for fscanf */
	char 	var_name[100];
	double 	var_value;
	
	/* orbital parameters given by user */
	double 	e = 0.0, a = 0.0, T = 0.0;
	/* state variables given by user*/
	double 	b0_diag[] = {0.0, 0.0, 0.0};
	/* non-state variables given by user */
	double 	omega[] = {0.0, 0.0, 0.0};
	/* system parameters given by user */
	int		elements = 0; //number of Voigt elements
	double  G = 0.0, m1 = 0.0, m2 = 0.0, I0 = 0.0, R = 0.0;
	double	alpha_0 = 0.0, alpha = 0.0, eta = 0.0;

	/* verification variables for system input */
	int 	number_system_inputs = 18;
	bool	input_system_received[number_system_inputs];
	for (int i = 0; i < number_system_inputs; i++)
	{
		input_system_received[i] = false;
	}

	/* reading system specs from user */
	FILE *in1 = fopen(argv[1], "r");
	if	(in1 == NULL)
	{
		fprintf(stderr, "Warning: could not read system input file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	while(fscanf(in1, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "e") == 0)
		{
			e = var_value;
			input_system_received[0] = true;
		}					
		else if (strcmp(var_name, "a") == 0)
		{
			a = var_value;
			input_system_received[1] = true;
		}
		else if (strcmp(var_name, "T") == 0)
		{
			T = var_value;
			input_system_received[2] = true;
		}
		else if (strcmp(var_name, "G") == 0)
		{
			G = var_value;
			input_system_received[3] = true;
		}
		else if (strcmp(var_name, "m1") == 0)
		{
			m1 = var_value;
			input_system_received[4] = true;
		}
		else if (strcmp(var_name, "m2") == 0)
		{
			m2 = var_value;
			input_system_received[5] = true;
		}
		else if (strcmp(var_name, "I0") == 0)
		{
			I0 = var_value;
			input_system_received[6] = true;
		}
		else if (strcmp(var_name, "R") == 0)
		{
			R = var_value;
			input_system_received[7] = true;
		}
		else if (strcmp(var_name, "b0_x") == 0)
		{
			b0_diag[0] = var_value;
			input_system_received[8] = true;
		}
		else if (strcmp(var_name, "b0_y") == 0)
		{
			b0_diag[1] = var_value;
			input_system_received[9] = true;
		}
		else if (strcmp(var_name, "b0_z") == 0)
		{
			b0_diag[2] = var_value;
			input_system_received[10] = true;
		}
		else if (strcmp(var_name, "omega_x") == 0)
		{
			omega[0] = var_value;
			input_system_received[11] = true;
		}
		else if (strcmp(var_name, "omega_y") == 0)
		{
			omega[1] = var_value;
			input_system_received[12] = true;
		}
		else if (strcmp(var_name, "omega_z") == 0)
		{
			omega[2] = var_value;
			input_system_received[13] = true;
		}
		else if (strcmp(var_name, "alpha_0") == 0)
		{
			alpha_0 = var_value;
			input_system_received[14] = true;
		}
		else if (strcmp(var_name, "alpha") == 0)
		{
			alpha = var_value;
			input_system_received[15] = true;
		}
		else if (strcmp(var_name, "eta") == 0)
		{
			eta = var_value;
			input_system_received[16] = true;
		}
		else if (strcmp(var_name, "elements") == 0)
		{
			elements = (int) var_value;
			input_system_received[17] = true;
		}
	}
	fclose(in1);

	/* system input verification */
	for (int i = 0; i < number_system_inputs; i++)
	{
		if(input_system_received[i] == false)
		{
			fprintf(stderr, "Error: there is at least one missing input ");
			fprintf(stderr, "from %s.\n", argv[1]);
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}
	}
	if (fabs(alpha) < 1e-14 || fabs(eta) < 1e-14)
	{
		fprintf(stderr, "Warning: nor alpha nor eta should be zero.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}

	/* setting up Voigt elements */
	double	*alpha_elements = *(&alpha_elements);
	double	*eta_elements = *(&eta_elements);
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
		while(fscanf(in1_elements, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
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
				fprintf(stderr, "Error: parameters missing for Voigt elements.\n");
				exit(10);
			}
			if (alpha_elements[i] < 1e-13 || eta_elements[i] < 1e-13)
			{
				fprintf(stderr, "Warning: nor alpha nor eta should be zero.\n");
				exit(13);
			}
		}
	}

	/* integrator parameters */
	double 	h = 0.0, t0 = 0.0, t1 = 0.0;
	double	eps_abs = 0.0, eps_rel = 0.0;
	int		data_step = 0;

	/* verification variables for integrator input */
	int 	number_integrator_inputs = 6;
	bool	input_integrator_received[number_integrator_inputs];
	for (int i = 0; i < number_integrator_inputs; i++)
	{
		input_integrator_received[i] = false;
	}

	/* reading integrator specs from user */
	FILE *in2 = fopen(argv[2], "r");
	if	(in2 == NULL)
	{
		fprintf(stderr, "Warning: could not read integrator input file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	while(fscanf(in2, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "h") == 0)
		{
			h = var_value;
			input_integrator_received[0] = true;
		}
		else if (strcmp(var_name, "t0") == 0)
		{
			t0 = var_value;
			input_integrator_received[1] = true;
		}
		else if (strcmp(var_name, "t1") == 0)
		{
			t1 = var_value;
			input_integrator_received[2] = true;
		}
		else if (strcmp(var_name, "eps_abs") == 0)
		{
			eps_abs = var_value;
			input_integrator_received[3] = true;
		}
		else if (strcmp(var_name, "eps_rel") == 0)
		{
			eps_rel = var_value;
			input_integrator_received[4] = true;
		}
		else if (strcmp(var_name, "data_step") == 0)
		{
			data_step = (int) var_value;
			input_integrator_received[5] = true;
		}
	}
	fclose(in2);

	/* parameter input verification */
	for (int i = 0; i < number_integrator_inputs; i++)
	{
		if(input_integrator_received[i] == false)
		{
			fprintf(stderr, "Error: there is at least one missing input ");
			fprintf(stderr, "from %s.\n", argv[2]);
			exit(13);
		}
	}

	/* make a copy of the input files */
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}
	FILE *in1_to_copy = fopen(argv[1], "r");
	FILE *in1_copy = fopen("output/input_pars_system_copy.txt" , "w");
	char ch = fgetc(in1_to_copy);
    while(ch != EOF)
    {
        fputc(ch, in1_copy);
        ch = fgetc(in1_to_copy);
    }
	fclose(in1_to_copy);
	fclose(in1_copy);
	FILE *in2_to_copy = fopen(argv[2], "r");
	FILE *in2_copy = fopen("output/input_pars_integrator_copy.txt" , "w");
	char ch2 = fgetc(in2_to_copy);
    while(ch2 != EOF)
    {
        fputc(ch2, in2_copy);
        ch2 = fgetc(in2_to_copy);
    }
	fclose(in2_to_copy);
	fclose(in2_copy);

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

	/* variables not given by user */
	double gamma = parameter_gamma_homogeneous_body(G, I0, R);
	double u_me[5], *bk_me = *(&bk_me);
	for (int i = 0; i < 5; i++) u_me[i] = 0.0;
	if (elements > 0)
	{
		bk_me = (double *) calloc(elements * 5, sizeof(double));
	}

	/* for testing */
	// for (int i = 0; i < 5; i++)
	// {
	// 	u_me[i] = (double) (-1 * i);
	// }
	// for (int i = 0; i < elements * 5; i++)
	// {
	// 	bk_me[i] = (double) i;
	// }

	/* calculate first l */
	double 	b[9], l[3];
	calculate_b(b, G, m2, gamma, alpha_0, alpha,
		tilde_x, omega, b0_me, u_me, elements, bk_me);
	calculate_l(l, I0, b, omega);

	/* for testing */
	// printf("b = \n");
	// print_square_matrix(b);
	// printf("\nomega = \n");
	// print_vector(omega);
	// exit(42);
	// null_matrix(b);

	/* variables and parameters passed as field parameters */
	int		dim_params = 12 + (elements * 2); // optmize this
	double	params[dim_params]; 
	params[0] = omega[0];
	params[1] = omega[1];
	params[2] = omega[2];
	params[3] = G;
	params[4] = m1;
	params[5] = m2;
	params[6] = I0;
	params[7] = gamma;
	params[8] = alpha;
	params[9] = eta;
	params[10] = alpha_0;
	params[11] = elements;
	for (int i = 0; i < elements; i++)
	{
		params[12 + (2*i)] = alpha_elements[i];
		params[13 + (2*i)] = eta_elements[i];
	}

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

	/* integration loop */
	int counter = 0;
	while (t < t1)
	// while (counter < 1) // for testing
	{
		/* for testing */
		// printf("omega 1 = \n");
		// print_vector(omega);

		if (fabs(t1-t) < h) h = t1-t;	//smaller last step

	  	gsl_odeiv2_system sys = {field_1EB1PM, NULL, dim, params};
	
		int status = 
			gsl_odeiv2_evolve_apply_fixed_step (ode_evolve, 
				ode_control, ode_step, &sys, &t, h, y);

		if (status != GSL_SUCCESS)
		{
			fprintf(stderr, "Error: GSL odeiv2 status = %d\n", status);
			break;
		}

		/* update variables */
		for (int i = 0; i < 3; i++) 				tilde_x[i] = y[0 + i];
		for (int i = 0; i < 3; i++) 				tilde_x_dot[i] = y[3 + i];
		for (int i = 0; i < 3; i++) 				l[i] = y[6 + i];
		for (int i = 0; i < 5; i++) 				b0_me[i] = y[9 + i];
		for (int i = 0; i < 5; i++) 				u_me[i] = y[14 + i];
		for (int i = 0; i < (elements * 5); i++) 	bk_me[i] = y[19 + i];

		/* update omega */
		double omega_seed[3];
		copy_vector(omega_seed, omega);
		calculate_omega(omega, omega_seed, G, m2, I0, gamma, alpha_0, 
			alpha, tilde_x, l, b0_me, u_me, elements, bk_me);
		for (int i = 0; i < 3; i++) params[i] = omega[i];

		/* calculate total angular momentum */
		double l_total[3];
		total_angular_momentum(l_total, m1, m2, 
			tilde_x, tilde_x_dot, l);

		/* writes output */
		if (counter % data_step == 0)
		{
			printf ("%.5e %.5e %.5e %.5e %.5e %.5e\n", 
				t, y[0], y[1],
				l_total[0], l_total[1], l_total[2]);
		}

		/* for testing */
		// printf("omega 2 = \n");
		// print_vector(omega);

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
