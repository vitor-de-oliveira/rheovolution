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
#include "celmec.h"

#define PI 3.14159265358979323846
#define G 1.0
// #define G 6.6743e-11

int
main(int argc, char *argv[]) 
{
	/* check number of command line arguments */
	if( argc < 4 )
	{
	   	fprintf(stderr, "Please, provide the type of system parameters file,");
		fprintf(stderr, " and both the system parameters");
		fprintf(stderr, " and the integrator parameters files\n");
	   	exit(2);
   	}
   	else if ( argc > 4 )
   	{
	   	fprintf(stderr, "Too many arguments\n");
	   	exit(4);
   	}

	/* auxiliary variables for fscanf */
	char 	var_name[100];
	double 	var_value;
	
	/* orbital parameters given by user */
	double 	e = 0.0, a = 0.0;
	/* state variables given by user*/
	double 	b0_diag[] = {0.0, 0.0, 0.0};
	/* non-state variables given by user */
	double 	omega[] = {0.0, 0.0, 0.0};
	/* system parameters given by user */
	int		elements = 0; //number of Voigt elements
	double  m1 = 0.0, m2 = 0.0;
	double	I0 = 0.0, R = 0.0, kf = 0.0;
	double	alpha = 0.0, eta = 0.0, alpha_0 = 0.0;
	/* Voigt elements for Maxwell generalized rheology */
	double	*alpha_elements = *(&alpha_elements);
	double	*eta_elements = *(&eta_elements);
	/* position and velocity */
	double	tilde_x[] = {0.0, 0.0, 0.0};
	double 	tilde_x_dot[] = {0.0, 0.0, 0.0};
	
	if (atoi(argv[1]) == 1)
	{
		/* verification variables for system input */
		int 	number_system_inputs = 17;
		bool	input_system_received[number_system_inputs];
		for (int i = 0; i < number_system_inputs; i++)
		{
			input_system_received[i] = false;
		}

		/* reading system specs from user */
		FILE *in1 = fopen(argv[2], "r");
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
			else if (strcmp(var_name, "m1") == 0)
			{
				m1 = var_value;
				input_system_received[2] = true;
			}
			else if (strcmp(var_name, "m2") == 0)
			{
				m2 = var_value;
				input_system_received[3] = true;
			}
			else if (strcmp(var_name, "I0") == 0)
			{
				I0 = var_value;
				input_system_received[4] = true;
			}
			else if (strcmp(var_name, "R") == 0)
			{
				R = var_value;
				input_system_received[5] = true;
			}
			else if (strcmp(var_name, "kf") == 0)
			{
				kf = var_value;
				input_system_received[6] = true;
			}
			else if (strcmp(var_name, "b0_x") == 0)
			{
				b0_diag[0] = var_value;
				input_system_received[7] = true;
			}
			else if (strcmp(var_name, "b0_y") == 0)
			{
				b0_diag[1] = var_value;
				input_system_received[8] = true;
			}
			else if (strcmp(var_name, "b0_z") == 0)
			{
				b0_diag[2] = var_value;
				input_system_received[9] = true;
			}
			else if (strcmp(var_name, "omega_x") == 0)
			{
				omega[0] = var_value;
				input_system_received[10] = true;
			}
			else if (strcmp(var_name, "omega_y") == 0)
			{
				omega[1] = var_value;
				input_system_received[11] = true;
			}
			else if (strcmp(var_name, "omega_z") == 0)
			{
				omega[2] = var_value;
				input_system_received[12] = true;
			}
			else if (strcmp(var_name, "alpha_0") == 0)
			{
				alpha_0 = var_value;
				input_system_received[13] = true;
			}
			else if (strcmp(var_name, "alpha") == 0)
			{
				alpha = var_value;
				input_system_received[14] = true;
			}
			else if (strcmp(var_name, "eta") == 0)
			{
				eta = var_value;
				input_system_received[15] = true;
			}
			else if (strcmp(var_name, "elements") == 0)
			{
				elements = (int) var_value;
				input_system_received[16] = true;
			}
		}
		fclose(in1);

		/* system input verification */
		for (int i = 0; i < number_system_inputs; i++)
		{
			if(input_system_received[i] == false)
			{
				fprintf(stderr, "Error: there is at least one missing input ");
				fprintf(stderr, "from %s.\n", argv[2]);
				exit(13);
			}
		}
		if (fabs(alpha) < 1e-14 || fabs(eta) < 1e-14)
		{
			fprintf(stderr, "Error: nor alpha nor eta should be zero.\n");
			exit(13);
		}

		/* setting up Voigt elements */
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

			FILE *in1_elements = fopen(argv[2], "r");
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

		/* position and velocity at periapsis given by Murray */
		tilde_x[0] 		= a * (1.0 - e);
		tilde_x[1] 		= 0.0;
		tilde_x[2] 		= 0.0;	
		tilde_x_dot[0]	= 0.0;
		tilde_x_dot[1] 	= ((2.0 * PI) / kepler_period(m1, m2, G, a)) 
							* a * sqrt((1.0 + e)/(1.0 - e));
		tilde_x_dot[2] 	= 0.0;

		/* for testing */
		// print_vector(tilde_x);
		// print_vector(tilde_x_dot);
		// exit(99);

	}
	else if (atoi(argv[1]) == 2) /* convert input if necessary */
	{
		convert_input(&m1, &m2, &I0, &R, &kf,
			omega, &alpha, &eta,
			tilde_x, tilde_x_dot,
			G,
			argv[2]);
		// exit(99);
	}
	else
	{
	   	fprintf(stderr, "Type of system file must be 1 or 2.\n");
	   	exit(5);
	}

	/* for testing */
	// print_vector(tilde_x);
	// print_vector(tilde_x_dot);
	// print_vector(omega);
	// print_vector(b0_diag);
	// printf("%1.5e\n", m1);
	// printf("%1.5e\n", m2);
	// printf("%1.5e\n", I0);
	// printf("%1.5e\n", R);
	// printf("%1.5e\n", kf);
	// printf("%1.5e\n", alpha);
	// printf("%1.5e\n", eta);
	// exit(99);

	/* integrator parameters */
	double	t_init = 0.0, t_trans = 0.0, t_final = 0.0;
	double 	t_step = 0.0;
	double	eps_abs = 0.0, eps_rel = 0.0;
	int		data_skip = 0;

	/* verification variables for integrator input */
	int 	number_integrator_inputs = 7;
	bool	input_integrator_received[number_integrator_inputs];
	for (int i = 0; i < number_integrator_inputs; i++)
	{
		input_integrator_received[i] = false;
	}

	/* reading integrator specs from user */
	FILE *in2 = fopen(argv[3], "r");
	if	(in2 == NULL)
	{
		fprintf(stderr, "Warning: could not read integrator input file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	while(fscanf(in2, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "t_init") == 0)
		{
			t_init = var_value;
			input_integrator_received[0] = true;
		}
		else if (strcmp(var_name, "t_trans") == 0)
		{
			t_trans = var_value;
			input_integrator_received[1] = true;
		}
		else if (strcmp(var_name, "t_final") == 0)
		{
			t_final = var_value;
			input_integrator_received[2] = true;
		}
		else if (strcmp(var_name, "t_step") == 0)
		{
			t_step = var_value;
			input_integrator_received[3] = true;
		}
		else if (strcmp(var_name, "eps_abs") == 0)
		{
			eps_abs = var_value;
			input_integrator_received[4] = true;
		}
		else if (strcmp(var_name, "eps_rel") == 0)
		{
			eps_rel = var_value;
			input_integrator_received[5] = true;
		}
		else if (strcmp(var_name, "data_skip") == 0)
		{
			data_skip = (int) var_value;
			input_integrator_received[6] = true;
		}
	}
	fclose(in2);

	/* parameter input verification */
	for (int i = 0; i < number_integrator_inputs; i++)
	{
		if(input_integrator_received[i] == false)
		{
			fprintf(stderr, "Error: there is at least one missing input ");
			fprintf(stderr, "from %s.\n", argv[3]);
			exit(13);
		}
	}

	/* make a copy of the input files */
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}
	FILE *in1_to_copy = fopen(argv[2], "r");
	FILE *in1_copy = fopen("output/input_pars_system_copy.txt" , "w");
	char ch = fgetc(in1_to_copy);
    while(ch != EOF)
    {
        fputc(ch, in1_copy);
        ch = fgetc(in1_to_copy);
    }
	fclose(in1_to_copy);
	fclose(in1_copy);
	FILE *in2_to_copy = fopen(argv[3], "r");
	FILE *in2_copy = fopen("output/input_pars_integrator_copy.txt" , "w");
	char ch2 = fgetc(in2_to_copy);
    while(ch2 != EOF)
    {
        fputc(ch2, in2_copy);
        ch2 = fgetc(in2_to_copy);
    }
	fclose(in2_to_copy);
	fclose(in2_copy);

	/* complete b0_me */
	double b0_me[5];
	b0_me[0] = b0_diag[0];
	b0_me[1] = 0.0;
	b0_me[2] = 0.0;
	b0_me[3] = b0_diag[1];
	b0_me[4] = 0.0;

	/* variables not given by user */
	double gamma = parameter_gamma(G, I0, R, kf);
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

	/* calculate I, J2, and C22 */
	// I do not use these as of yet
	// double I[9];
	// calculate_inertia_tensor(I, I0, b);
	// double J2;
	// J2 = calculate_J2(m1, R, I);
	// double C22;
	// C22 = calculate_C22(m1, R, I);

	/* for testing */
	// printf("I = \n");
	// print_square_matrix(I);
	// printf("J2 = %e\n", J2);
	// printf("C22 = %e\n", C22);
	// exit(43);

	/* for testing */
	// printf("tilde_x = \n");
	// print_vector(tilde_x);
	// printf("\ntilde_x_dot = \n");
	// print_vector(tilde_x_dot);
	// printf("\nl = \n");
	// print_vector(l);
	// printf("\nb0_diag = \n");
	// print_vector(b0_diag);
	// printf("\nomega = \n");
	// print_vector(omega);
	// printf("\nG = %e\n", G);
	// printf("\nm1 = %e\n", m1);
	// printf("\nm2 = %e\n", m2);
	// printf("\nI0 = %e\n", I0);
	// printf("\ngamma = %e\n", gamma);
	// printf("\nalpha = %e\n", alpha);
	// printf("\neta = %e\n", eta);
	// printf("\nalpha_0 = %e\n", alpha_0);
	// for (int i  = 0; i < elements; i++)
	// {
	// 	printf("\nalpha_%d = %f\n", i+1, alpha_elements[i]);
	// 	printf("\neta%d = %f\n", i+1, eta_elements[i]);
	// }
	// exit(42);

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
	double 	t = t_init;
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
	while (t < t_final)
	// while (counter < 1) // for testing
	{
		/* for testing */
		// printf("omega 1 = \n");
		// print_vector(omega);

		if (fabs(t_final-t) < t_step)
		{
			t_step = t_final - t; //smaller last step
		}

	  	gsl_odeiv2_system sys = {field_1EB1PM, NULL, dim, params};
	
		int status = 
			gsl_odeiv2_evolve_apply_fixed_step (ode_evolve, 
				ode_control, ode_step, &sys, &t, t_step, y);

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

		/* calculate b */
		calculate_b(b, G, m2, gamma, alpha_0, alpha,
			tilde_x, omega, b0_me, u_me, elements, bk_me);

		/* calculate total angular momentum */
		double l_total[3];
		total_angular_momentum(l_total, m1, m2, 
			tilde_x, tilde_x_dot, l);

		/* writes output */
		if (t > t_trans)
		{
			if (counter % data_skip == 0)
			{
				/* time, position, and velocity */
				printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e", 
					t, tilde_x[0], tilde_x[1], tilde_x[2],
					tilde_x_dot[0], tilde_x_dot[1], tilde_x_dot[2]);
				/* angular momentum and total angular momentum */
				printf (" %.5e %.5e %.5e %.5e %.5e %.5e", 
					l[0], l[1], l[2],
					l_total[0], l_total[1], l_total[2]);
				/* deformation matrix */
				printf (" %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e",
					b[0], b[1], b[2],
					b[3], b[4], b[5],
					b[6], b[7], b[8]);
				printf("\n");
			}
		}

		/* for testing */
		// printf("omega 2 = \n");
		// print_vector(omega);

		/* for testing */
		// printf("tilde_x = \n");
		// print_vector(tilde_x);
		// printf("tilde_x_dot = \n");
		// print_vector(tilde_x_dot);
		// printf("l = \n");
		// print_vector(l);
		// printf("b0_me = \n");
		// for (int i  = 0; i < 5; i++) printf("%e ", b0_me[i]);
		// printf("\n");
		// printf("u_me = \n");
		// for (int i  = 0; i < 5; i++) printf("%e ", u_me[i]);
		// printf("\n");
		// printf("bk_me = \n");
		// for (int i = 0; i < (elements * 5); i++) printf("%e ", bk_me[i]);
		// printf("\n");
		// printf("omega = \n");
		// print_vector(omega);
		// printf("G = %e\n", G);
		// printf("m1 = %e\n", m1);
		// printf("m2 = %e\n", m2);
		// printf("I0 = %e\n", I0);
		// printf("gamma = %e\n", gamma);
		// printf("alpha = %e\n", alpha);
		// printf("eta = %e\n", eta);
		// printf("alpha_0 = %e\n", alpha_0);
		// for (int i  = 0; i < elements; i++)
		// {
		// 	printf("alpha_%d = %f\n", i+1, alpha_elements[i]);
		// 	printf("eta%d = %f\n", i+1, eta_elements[i]);
		// }
		// printf("b = \n");
		// print_square_matrix(b);

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
