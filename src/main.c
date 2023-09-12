#define _GNU_SOURCE // getline
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h> // ssize_t

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "celmec.h"
#include "dynamical_system.h"
#include "parsing.h"

#define t(n) printf("Here %d\n", n)

int
main(int argc, char *argv[]) 
{
	/* variables for parsing input file */
    char 	*line = NULL;
    size_t 	len = 0;
    ssize_t	read;
	char 	first_col[100];
	char 	second_col[100];
	char	sim_name[100];			// simulation name
	char	input_folder[100];
	char	output_folder[100];
	char	system_specs[100];
	char	integrator_specs[100];
	char	dev_specs[100];
	int		system_file_type = 1;
	int		number_of_bodies = 0;
	bool	system_file_type_check = false;
	bool	number_of_bodies_check = false;
	FILE	*in = fopen(argv[1], "r");

	/* parse input file */
   	while ((read = getline(&line, &len, in)) != -1)
	{
		sscanf(line, "%s %s", first_col, second_col);
		if (strcmp(first_col, "name") == 0)
		{
			strcpy(sim_name, second_col);
		}
		else if (strcmp(first_col, "input_folder") == 0)
		{
			strcpy(input_folder, second_col);
		}
		else if (strcmp(first_col, "output_folder") == 0)
		{
			strcpy(output_folder, second_col);
			// create output folder if it does not exist
			struct stat st = {0};
			if (stat(output_folder, &st) == -1) {
				mkdir(output_folder, 0700);
			}
		}
		else if (strcmp(first_col, "system_specs") == 0)
		{
			strcpy(system_specs, input_folder);
			strcat(system_specs, second_col);
		}
		else if (strcmp(first_col, "integrator_specs") == 0)
		{
			strcpy(integrator_specs, input_folder);
			strcat(integrator_specs, second_col);
		}
		else if (strcmp(first_col, "dev_specs") == 0)
		{
			strcpy(dev_specs, second_col);
		}
		else if (strcmp(first_col, "system_file_type") == 0)
		{
			system_file_type = atoi(second_col);
			system_file_type_check = true;
		}
		else if (strcmp(first_col, "number_of_bodies") == 0)
		{
			number_of_bodies = atoi(second_col);
			number_of_bodies_check = true;
		}
	}

	/* checking if file type and number of bodies were received */
	if (system_file_type_check == false)
	{
		fprintf(stderr, "Error: Please provide type of system file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(15);
	}
	else if (system_file_type != 1 && system_file_type != 2)
	{
		fprintf(stderr, "Error: Type of system file should be 1 or 2.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(15);
	}
	if (number_of_bodies_check == false)
	{
		fprintf(stderr, "Error: Please provide number of bodies.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(15);
	}
	else if (number_of_bodies < 1)
	{
		fprintf(stderr, "Error: Number of bodies should be greater than 0.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(15);
	}

	/* gravitational constant */
	double  G = 4.0 * M_PI * M_PI; // AU Msun year

	/* dev variables */
	char	units[100] = "AU_Msun_year";

	/* read dev file */
	FILE *in3 = fopen(dev_specs, "r");
	if	(in3 != NULL)
	{
		char 	*line = NULL;
    	size_t 	len = 0;
	    ssize_t read;
		char 	first_col[100];
		char 	second_col[100];
		while ((read = getline(&line, &len, in3)) != -1)
		{
			sscanf(line, "%s %s",
				first_col, second_col);
			if (strcmp(first_col, "Units") == 0)
			{
				if (strcmp(second_col, "SI") == 0)
				{
					G = 6.6743e-11; // SI
					strcpy(units, second_col);
				}
				else if (strcmp(second_col, "G_unity") == 0)
				{
					G = 1.0; // non-dimensional
					strcpy(units, second_col);
				}
			}
		}
		fclose(in3);
	}

	/* array of celestial bodies */
	cltbdy	*bodies;

	/* orbital parameters given by user */
	double 	e = 0.0, a = 0.0;
	/* state variables given by user*/
	double 	b0_diag[] = {0.0, 0.0, 0.0};
	/* non-state variables given by user */
	double 	omega[3];
	omega[0] = omega[0];
	omega[1] = omega[1];
	omega[2] = omega[2];
	/* system parameters given by user */
	int		elements = 0; // number of Voigt elements
	double  m1 = 0.0, m2 = 0.0;
	double	I0 = I0, R = R, kf = kf;
	double	alpha = 0.0, eta = 0.0, alpha_0 = 0.0;
	/* Voigt elements for Maxwell generalized rheology */
	double	*alpha_elements = *(&alpha_elements);
	double	*eta_elements = *(&eta_elements);
	/* position and velocity */
	double	tilde_x[3];
	tilde_x[0] = tilde_x[0];
	tilde_x[1] = tilde_x[1];
	tilde_x[2] = tilde_x[2];
	double 	tilde_x_dot[3];
	tilde_x_dot[0] = tilde_x_dot[0];
	tilde_x_dot[1] = tilde_x_dot[1];
	tilde_x_dot[2] = tilde_x_dot[2];
	/* deformation settings */
	// bool	centrifugal = false;
	// bool	tidal = false;

	/* auxiliary variables for fscanf */
	char 	var_name[100];
	double 	var_value;

	if (system_file_type == 1)
	{
		// Warning for dev
		fprintf(stderr, "Deformable settings not implemented");
		fprintf(stderr, " yet for this type of file!\n");
		exit(29);

		/* verification variables for system input */
		int 	number_system_inputs = 17;
		bool	input_system_received[number_system_inputs];
		for (int i = 0; i < number_system_inputs; i++)
		{
			input_system_received[i] = false;
		}

		/* reading system specs from user */
		FILE *in1 = fopen(system_specs, "r");
		if	(in1 == NULL)
		{
			fprintf(stderr, "Warning: could not read system specs file.\n");
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
				fprintf(stderr, "from %s.\n", system_specs);
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

			FILE *in1_elements = fopen(system_specs, "r");
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
		tilde_x_dot[1] 	= ((2.0 * M_PI) / kepler_period(m1, m2, G, a)) 
							* a * sqrt((1.0 + e)/(1.0 - e));
		tilde_x_dot[2] 	= 0.0;

		/* for testing */
		// print_vector(tilde_x);
		// print_vector(tilde_x_dot);
		// e = calculate_eccentricity(G, m1, m2, tilde_x, tilde_x_dot);
		// a = calculate_semi_major_axis(G, m1, m2, tilde_x, tilde_x_dot);
		// printf("e = %e a = %e\n", e, a);
		// exit(99);

	}
	else if (system_file_type == 2) /* convert input if necessary */
	{
		convert_input(&bodies,
			number_of_bodies, 
			G,
			system_specs,
			units);

		/* for testing */
		// printf("eta = %e alpha = %e tau = %e\n", eta, alpha, eta/alpha);
		// for (int i = 0; i < number_of_bodies; i++)
		// {
		// 	double relative_x[3];
		// 	linear_combination_vector(relative_x,
		// 		1.0, bodies[i].x,
		// 		-1.0, bodies[0].x);
		// 	double relative_x_dot[3];
		// 	linear_combination_vector(relative_x_dot,
		// 		1.0, bodies[i].x_dot,
		// 		-1.0, bodies[0].x_dot);
		// 	e = calculate_eccentricity(G, bodies[0].mass, 
		// 		bodies[i].mass, relative_x, relative_x_dot);
		// 	a = calculate_semi_major_axis(G, bodies[0].mass, 
		// 		bodies[i].mass, relative_x, relative_x_dot);
		// 	printf("bodies = %s e = %e a = %e\n", bodies[i].name, e, a);
		// }
		// exit(99);
	}
	else
	{
	   	fprintf(stderr, "Type of system file must be 1 or 2.\n");
	   	exit(5);
	}

	/* for testing */
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%s ", bodies[i].name);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].mass);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].lod);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].obl);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].psi);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].R);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].rg);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].J2);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].C22);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].lib);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].kf);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].Dt);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].tau);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].a);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].e);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].I);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].M);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].w);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].Omega);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", bodies[i].tidal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", bodies[i].centrifugal);
	// printf("\n");
	// exit(99);

	// /* Calibration and plots of k2 */
	// double rate = 2.54014310646e-13; // 3.8 cm/yr in AU/yr
	// double a0 = calculate_semi_major_axis(G, m1, m2, tilde_x, tilde_x_dot);
	// double Imk2;
	// Imk2 = calibrate_Imk2(rate, a0, m1, m2, I0, R, omega[2], G);
	// printf("b = %1.5e\n", Imk2);
	// double nu = 0.680; // tau_v = nu * tau
	// double tau_v_pair[2];
	// double tau_pair[2];
	// calculate_tau_v_and_tau(tau_v_pair, tau_pair, 
	// 	nu, Imk2, a0, m1, m2, kf, omega[2], G);
	// printf("tau_v_minus = %1.5e = %1.5e s tau_minus = %1.5e\n", tau_v_pair[0],
	// 	tau_v_pair[0] * (365.25 * 24.0 * 60.0 * 60.0), tau_pair[0]);
	// printf("tau_v_plus = %1.5e = %1.5e s tau_plus = %1.5e\n", tau_v_pair[1],
	// 	tau_v_pair[0] * (365.25 * 24.0 * 60.0 * 60.0), tau_pair[1]);
	// /* Love number as a function of frequency */
	// for (double sigma_loop = 1e-10; sigma_loop < 1e20; sigma_loop *= 1.01)
	// {
	// 	double real_loop, imag_loop;
	// 	/* measure k2 for lower tau */
	// 	calculate_k2(&real_loop, &imag_loop, sigma_loop, kf, 
	// 		tau_v_pair[0], tau_pair[0]);
	// 	printf("%1.15e %1.15e %1.15e %1.15e ", 
	// 		sigma_loop, real_loop, -1.0 * imag_loop,
	// 		sqrt(real_loop * real_loop + imag_loop * imag_loop));
	// 	/* measure k2 for higher tau */
	// 	calculate_k2(&real_loop, &imag_loop, sigma_loop, kf, 
	// 		tau_v_pair[1], tau_pair[1]);
	// 	printf("%1.15e %1.15e %1.15e\n", 
	// 		real_loop, -1.0 * imag_loop,
	// 		sqrt(real_loop * real_loop + imag_loop * imag_loop));
	// }
	// /* Chandler wobble */
	// double sigma_chandler = (2.0 * M_PI) / (433.0 / 365.25); // rad/yr
	// printf("Chandler wobble:\n");
	// printf("Chandler wobble frequency Cwf = %1.5e\n", sigma_chandler);
	// double real, imag;
	// calculate_k2(&real, &imag, sigma_chandler, kf, 
	// 	tau_v_pair[0], tau_pair[0]);
	// printf("tau = %1.5e\n", tau_pair[0]);
	// printf("Re(k2) = %1.5e Im(k2) = %1.5e |k_2| = %1.5e\n", 
	// 	real, imag, sqrt(real*real+imag*imag));
	// calculate_k2(&real, &imag, sigma_chandler, kf, 
	// 	tau_v_pair[1], tau_pair[1]);
	// printf("tau = %1.5e\n", tau_pair[1]);
	// printf("Re(k2) = %1.5e Im(k2) = %1.5e |k_2| = %1.5e\n", 
	// 	real, imag, sqrt(real*real+imag*imag));
	// exit(182);

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
	FILE *in2 = fopen(integrator_specs, "r");
	if	(in2 == NULL)
	{
		fprintf(stderr, "Warning: could not read integrator specs file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
   	while ((read = getline(&line, &len, in2)) != -1)
	{
		sscanf(line, "%s %lf", var_name, &var_value);
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
			fprintf(stderr, "from %s.\n", integrator_specs);
			exit(13);
		}
	}

	/* make a copy of the input files */
	// struct stat st = {0};
	// if (stat("output", &st) == -1) {
	// 	mkdir("output", 0700);
	// }
	// FILE *in1_to_copy = fopen(system_specs, "r");
	// FILE *in1_copy = fopen("output/input_pars_system_copy.txt" , "w");
	// char ch = fgetc(in1_to_copy);
    // while(ch != EOF)
    // {
    //     fputc(ch, in1_copy);
    //     ch = fgetc(in1_to_copy);
    // }
	// fclose(in1_to_copy);
	// fclose(in1_copy);
	// FILE *in2_to_copy = fopen(integrator_specs, "r");
	// FILE *in2_copy = fopen("output/input_pars_integrator_copy.txt" , "w");
	// char ch2 = fgetc(in2_to_copy);
    // while(ch2 != EOF)
    // {
    //     fputc(ch2, in2_copy);
    //     ch2 = fgetc(in2_to_copy);
    // }
	// fclose(in2_to_copy);
	// fclose(in2_copy);

	/* complete b0_me */
	// I am using same b0_diag and alpha_0 for every body for now!
	// double b0_me[5];
	for (int i = 0; i < number_of_bodies; i++)
	{
		bodies[i].b0_me[0] = b0_diag[0];
		bodies[i].b0_me[1] = 0.0;
		bodies[i].b0_me[2] = 0.0;
		bodies[i].b0_me[3] = b0_diag[1];
		bodies[i].b0_me[4] = 0.0;
		bodies[i].alpha_0 = alpha_0;
	}

	/* variables not given by user */
	// double u_me[5], *bk_me = *(&bk_me);
	for (int i  = 0; i < number_of_bodies; i++)
	{
		for (int j = 0; j < 5; j++) bodies[i].u_me[i] = 0.0;
		if (bodies[i].elements > 0)
		{
			bodies[i].bk_me = (double *) calloc(bodies[i].elements * 5, sizeof(double));
		}
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
	// double 	b[9], l[3];
	for (int i = 0; i < number_of_bodies; i++)
	{
		calculate_b(i, bodies, number_of_bodies, G);
		calculate_l(i, bodies, number_of_bodies, G);
		// printf("omega of body %d before\n", i+1);
		// print_vector(bodies[i].omega);
		// calculate_omega(i, bodies, number_of_bodies, G);
		// printf("omega of body %d after\n", i+1);
		// print_vector(bodies[i].omega);
	}

	// exit(98);

	/* for testing */
	// double omega_seed_test[3];
	// copy_vector(omega_seed_test, omega);
	// calculate_omega(omega, omega_seed_test, G, m2, I0, gamma, alpha_0, 
	// 	alpha, tilde_x, l, b0_me, u_me, elements, bk_me);
	// print_vector(l);
	// print_vector(omega);
	// exit(99);

	/* for testing */
	// printf("b = \n");
	// print_square_matrix(b);
	// printf("\nomega = \n");
	// print_vector(omega);
	// exit(42);
	// null_matrix(b);

	/* calculate I, J2, and C22 */
	// I do not use these as of yet
	// double Inert[9];
	// calculate_inertia_tensor(Inert, I0, b);
	// double J2;
	// J2 = calculate_J2(m1, R, Inert);
	// double C22;
	// C22 = calculate_C22(m1, R, Inert);

	/* for testing */
	// printf("I = \n");
	// print_square_matrix(Inert);
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

	/* variables and parameters passed as field arguments */

	/* total number of Voigt elements */
	int	elements_total = 0;
	for (int i = 0; i < number_of_bodies; i++)
	{
		elements_total += bodies[i].elements;
	}
	
	/* state variables */
	int		dim_state_per_body_without_elements = 19;
	int		dim_state = (dim_state_per_body_without_elements * number_of_bodies) + (5 * elements_total);
	double 	y[dim_state];
	int		elements_counter = 0;
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_state_skip = i * dim_state_per_body_without_elements + 5 * elements_counter;
		for (int j = 0; j < 3; j++)
		{
			y[0 + dim_state_skip + j] 	= bodies[i].x[j];
			y[3 + dim_state_skip + j] 	= bodies[i].x_dot[j];
			y[6 + dim_state_skip + j] 	= bodies[i].l[j];
		}
		for (int j = 0; j < 5; j++)
		{
			y[9 + dim_state_skip + j] 	= bodies[i].b0_me[j];
			y[14 + dim_state_skip + j] 	= bodies[i].u_me[j];
		}
		for (int j = 0; j < 5 * bodies[i].elements; j++)
		{
			y[19 + dim_state_skip + j] 	= bodies[i].bk_me[j];
		}
		elements_counter += bodies[i].elements;

		/* for testing */
		// printf("Body %d\n", i+1);
		// printf("x = \n");
		// print_vector(bodies[i].x);
		// printf("x_dot = \n");
		// print_vector(bodies[i].x_dot);
		// printf("l = \n");
		// print_vector(bodies[i].l);
		// double b0_print[9], u_print[9];
		// construct_traceless_symmetric_matrix(b0_print, bodies[i].b0_me);
		// construct_traceless_symmetric_matrix(u_print, bodies[i].u_me);
		// printf("b0 = \n");
		// print_square_matrix(b0_print);
		// printf("u = \n");
		// print_square_matrix(u_print);
		// exit(99);
	
	}

	// for (int i = 0; i < dim_state; i++)
	// {
	// 	printf("y[%d] = %1.10e\n", i, y[i]);
	// }
	// exit(98);

	/* parameters */
	int		dim_params_per_body_without_elements = 12;
	int		dim_params = 2 + (dim_params_per_body_without_elements * number_of_bodies) + (2 * elements_total);
	double	params[dim_params];
	elements_counter = 0; 
	params[0] = G;
	params[1] = (double) number_of_bodies;
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		params[2 + 0 + dim_params_skip] = (double) bodies[i].centrifugal;
		params[2 + 1 + dim_params_skip] = (double) bodies[i].tidal;
		for (int j = 0; j < 3; j++)
		{
			params[2 + 2 + dim_params_skip + j] = bodies[i].omega[j];
		}
		params[2 + 5 + dim_params_skip] = bodies[i].mass;
		params[2 + 6 + dim_params_skip] = bodies[i].I0;
		params[2 + 7 + dim_params_skip] = bodies[i].gamma;
		params[2 + 8 + dim_params_skip] = bodies[i].alpha;
		params[2 + 9 + dim_params_skip] = bodies[i].eta;
		params[2 + 10 + dim_params_skip] = bodies[i].alpha_0;
		params[2 + 11 + dim_params_skip] = (double) bodies[i].elements;
		for (int j = 0; j < bodies[i].elements; j++)
		{
			params[2 + 12 + dim_params_skip + j] = bodies[i].alpha_elements[j];
		}
		for (int j = 0; j < bodies[i].elements; j++)
		{
			params[2 + 13 + dim_params_skip + j + bodies[i].elements - 1] = bodies[i].eta_elements[j];
		}
		elements_counter += bodies[i].elements;		
	}

	// for (int i = 0; i < dim_state; i++)
	// {
	// 	printf("y[%d] = %1.10e\n", i, y[i]);
	// }
	// for (int i = 0; i < dim_params; i++)
	// {
	// 	printf("%f\n", params[i]);
	// }
	// exit(98);

	/* GSL variables */
	const gsl_odeiv2_step_type * ode_type
		= gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step * ode_step
    	= gsl_odeiv2_step_alloc (ode_type, dim_state);
  	gsl_odeiv2_control * ode_control
    	= gsl_odeiv2_control_y_new (eps_abs, eps_rel);
  	gsl_odeiv2_evolve * ode_evolve
    	= gsl_odeiv2_evolve_alloc (dim_state);

	// double e_init = calculate_eccentricity(G, m1, m2, tilde_x, tilde_x_dot);
	// double a_init = calculate_semi_major_axis(G, m1, m2, tilde_x, tilde_x_dot);
	// double norm_omega_init = norm_vector(omega);
	// double norm_l_init = norm_vector(l);
	// double l_total_init[3];
	// total_angular_momentum(l_total_init, m1, m2, tilde_x, tilde_x_dot, l);
	// double norm_l_total_init = norm_vector(l_total_init);
	// double I[9];
	// calculate_inertia_tensor(I, I0, b);
	// double J2_init = calculate_J2(m1, R, I);
	// double C22_init = calculate_C22(m1, R, I);

	// print_vector(omega);
	// printf("%.5e\n", norm_vector(omega));
	// exit(42);

	FILE *out[number_of_bodies + 1];
	char filename[150];
	for (int i = 0; i < number_of_bodies; i++)
	{
		/* names and opens output files for each body */
		strcpy(filename, output_folder);
		strcat(filename, "results_");
		strcat(filename, sim_name);
		strcat(filename, "_");
		strcat(filename, bodies[i].name);
		strcat(filename, ".dat");
		out[i] = fopen(filename, "w");
		/* writes output headers */
		fprintf(out[i], "time(yr)");
		fprintf(out[i], " x(AU) y(AU) z(AU) vx(Au/yr) vy(Au/yr) vz(Au/yr)");
		fprintf(out[i], " a(AU)");
		fprintf(out[i], " |omega|");
		fprintf(out[i], " b3");
		fprintf(out[i], " |b|");
		fprintf(out[i], "\n");
	}
	/* names and opens general output file */
	strcpy(filename, output_folder);
	strcat(filename, "results_");
	strcat(filename, sim_name);
	strcat(filename, "_");
	strcat(filename, "general");
	strcat(filename, ".dat");
	out[number_of_bodies] = fopen(filename, "w");
	fprintf(out[number_of_bodies], "time(yr)");
	fprintf(out[number_of_bodies], " |Total_angular_momentum|");
	fprintf(out[number_of_bodies], "\n");

	/* integration loop */
	int		counter = 0;	
	double 	t = t_init;
	while (t < t_final)
	// while (counter < 1) // for testing
	{
		/* for testing */
		// printf("omega 1 = \n");
		// print_vector(omega);
		// printf("%1.5e\n", t);

		if (fabs(t_final-t) < t_step)
		{
			t_step = t_final - t; //smaller last step
		}

	  	gsl_odeiv2_system sys = {field_GV, NULL, dim_state, params};
	
		int status = 
			gsl_odeiv2_evolve_apply_fixed_step (ode_evolve, 
				ode_control, ode_step, &sys, &t, t_step, y);

		if (status != GSL_SUCCESS)
		{
			fprintf(stderr, "Error: GSL odeiv2 status = %d\n", status);
			break;
		}

		/* update variables */
		elements_counter = 0; 
		for (int i = 0; i < number_of_bodies; i++)
		{
			int	dim_state_skip = i * dim_state_per_body_without_elements + 5 * elements_counter;
			
			for (int j = 0; j < 3; j++)
			{
				bodies[i].x[j] 		= y[0 + dim_state_skip + j];
				bodies[i].x_dot[j] 	= y[3 + dim_state_skip + j];
				bodies[i].l[j] 		= y[6 + dim_state_skip + j];
			}
			for (int j = 0; j < 5; j++)
			{
				bodies[i].b0_me[j] 	= y[9 + dim_state_skip + j];
				bodies[i].u_me[j] 	= y[14 + dim_state_skip + j];
			}
			if (bodies[i].elements > 0)
			{
				bodies[i].bk_me = (double *) malloc(5 * bodies[i].elements * sizeof(double));
				for (int j = 0; j < 5 * bodies[i].elements; j++)
				{
					bodies[i].bk_me[j] = y[19 + dim_state_skip + j];
				}
			}

			elements_counter += bodies[i].elements;
		}

		/* for testing */
		// print_vector(l);
		// printf("%.5e\n", norm_vector(l));
		// exit(42);

		/* update omega */
		elements_counter = 0;
		for (int i = 0; i < number_of_bodies; i++)
		{
			int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		
			calculate_omega(i, bodies, number_of_bodies, G);
		
			for (int j = 0; j < 3; j++)
			{
				params[2 + 2 + dim_params_skip + j] = bodies[i].omega[j];
			}

			elements_counter += bodies[i].elements;
		}


		/* for testing */
		// print_vector(omega);
		// printf("%.5e\n", norm_vector(omega));
		// exit(42);

		/* calculate b */
		for (int i = 0; i < number_of_bodies; i++)
		{
			calculate_b(i, bodies, number_of_bodies, G);
		}

		/* calculate total angular momentum */
		// double l_total[3];
		// total_angular_momentum(l_total, m1, m2, 
		// 	tilde_x, tilde_x_dot, l);

		/* calculate e and a */
		// e = calculate_eccentricity(G, m1, m2, tilde_x, tilde_x_dot);
		// a = calculate_semi_major_axis(G, m1, m2, tilde_x, tilde_x_dot);

		/* calculate differences */
		// double e_dif = e - e_init;
		// double a_dif = a - a_init;
		// double omega_dif = norm_vector(omega) - norm_omega_init;
		// double l_dif = norm_vector(l) - norm_l_init;
		// double l_total_dif = norm_vector(l_total) - norm_l_total_init;

		/* construct b0 and u for printing */
		// double b0[9], u[9];
		// construct_traceless_symmetric_matrix(b0, b0_me);
		// construct_traceless_symmetric_matrix(u, u_me);

		/* calculate J2 and C22 */
		// calculate_inertia_tensor(I, I0, b);
		// double J2 = calculate_J2(m1, R, I);
		// double C22 = calculate_C22(m1, R, I);
		// double J2_dif = J2 - J2_init;
		// double C22_dif = C22 - C22_init; 

		/* writes output */
		if (t > t_trans)
		{
			if (counter % data_skip == 0)
			{
				for (int i = 0; i < number_of_bodies; i++)
				{
					/* time */

					fprintf (out[i], "%.15e", t);

					/* position, and velocity */

					// transform to body 1-centered system
					double relative_x[3];
					linear_combination_vector(relative_x,
						1.0, bodies[i].x,
						-1.0, bodies[0].x);
					double relative_x_dot[3];
					linear_combination_vector(relative_x_dot,
						1.0, bodies[i].x_dot,
						-1.0, bodies[0].x_dot);

					fprintf (out[i], " %.15e %.15e %.15e %.15e %.15e %.15e", 
						relative_x[0], relative_x[1], relative_x[2],
						relative_x_dot[0], relative_x_dot[1], relative_x_dot[2]);

					/* orbital eccentricity and semimajor axis */

					// printf (" %.15e %.15e", e, e_dif);
					// printf (" %.15e %.15e", a, a_dif);

					double a_body =
						calculate_semi_major_axis(G, bodies[0].mass, bodies[i].mass, 
							relative_x, relative_x_dot);
					fprintf (out[i], " %.15e", a_body);

					/* angular velocity */

					fprintf (out[i], " %.15e", norm_vector(bodies[i].omega));
					// fprintf (out[i], " %.15e", omega_dif);
					// fprintf (out[i], " %.15e %.15e %.15e", 
					// 	omega[0], omega[1], omega[2]);

					/* angular momentum and total angular momentum */

					// printf (" %.15e %.15e", norm_vector(l), l_dif);
					// printf (" %.15e %.15e", norm_vector(l_total), l_total_dif);
					// printf (" %.15e %.15e %.15e %.15e %.15e %.15e", 
					// 	l[0], l[1], l[2],
					// 	l_total[0], l_total[1], l_total[2]);

					/* deformation matrix */

					fprintf (out[i], " %.15e", bodies[i].b[8]);
					fprintf (out[i], " %.15e", norm_square_matrix(bodies[i].b));
					// printf (" %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e",
					// 	b[0], b[1], b[2],
					// 	b[3], b[4], b[5],
					// 	b[6], b[7], b[8]);

					/* prestress */

					// printf (" %.15e", norm_square_matrix(b0));
					// printf (" %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e",
					// 	b0[0], b0[1], b0[2],
					// 	b0[3], b0[4], b0[5],
					// 	b0[6], b0[7], b0[8]);

					/* rheology */

					// printf (" %.15e", norm_square_matrix(u));
					// printf (" %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e",
					// 	u[0], u[1], u[2],
					// 	u[3], u[4], u[5],
					// 	u[6], u[7], u[8]);

					/* classical gravitational potential terms */

					// printf (" %.5e %.5e", J2, J2_dif);
					// printf (" %.5e %.5e", C22, C22_dif);

					fprintf(out[i], "\n");

				} // end loop over bodies

				/* print on general file */

				/* time */

				fprintf (out[number_of_bodies], "%.15e", t);

				/* print total angular momentum on general file */

				double l_total[3];
				total_angular_momentum(l_total, bodies, number_of_bodies, G);
				
				fprintf (out[number_of_bodies], " %.15e", norm_vector(l_total));

				fprintf(out[number_of_bodies], "\n");

			} // end if (counter % data_skip == 0)

			counter++;

		} // end if (t > t_trans)

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

	}

	/* free GSL variables */
	gsl_odeiv2_evolve_free (ode_evolve);
	gsl_odeiv2_control_free (ode_control);
	gsl_odeiv2_step_free (ode_step);

	for (int i = 0; i < number_of_bodies + 1; i++)
	{
		fclose(out[i]);
	}

	/* free Voigt elements */
	for (int i = 0; i < number_of_bodies; i++)
	{
		if (bodies[i].elements > 0)
		{
			free(bodies[i].alpha_elements);
			free(bodies[i].eta_elements);	
			free(bodies[i].bk_me);
		}
	}

	/* free array of celestial bodies */
	free(bodies);

	return 0;
}
