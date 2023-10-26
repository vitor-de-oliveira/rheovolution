#include "data_processing.h"

int
count_columns(const char *s)
{
    int columns = 0;
     
    while (isspace(*s))	// skip leading whitespace
    {
		s++;	// go to next character in string
	}
    while (*s != '\0')	// while not at end of string
	{
        columns++;
        while (*s != '\0' && !isspace(*s))	// while not spaces
		{    
            s++;	// go to next character in string
        }
        while (isspace(*s))	// skip inter-column whitespace
		{
            s++;
		}
	}
     
    return columns;
}

int
parse_input(siminf *simulation,
			const char input_file[])
{
	/* auxiliary variables */
    char 	*line = NULL;
    size_t 	len = 0;
    ssize_t	read;
	char 	first_col[100];
	char 	second_col[100];
	char	hold[100];
	double	second_col_double;
	bool	system_file_type_received = false;
	bool	number_of_bodies_received = false;
	bool	dev_specs_file_received = false;

	/* parse input file */
	FILE	*in = fopen(input_file, "r");
	if (in == NULL)
	{
		fprintf(stderr, "Warning: could not read input file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	else
	{
		strcpy((*simulation).main_input, input_file);
	}
   	while ((read = getline(&line, &len, in)) != -1)
	{
		sscanf(line, "%s %s", first_col, second_col);
		if (strcmp(first_col, "name") == 0)
		{
			strcpy((*simulation).name, second_col);
		}
		else if (strcmp(first_col, "input_folder") == 0)
		{
			strcpy((*simulation).input_folder, second_col);
		}
		else if (strcmp(first_col, "output_folder") == 0)
		{
			strcpy((*simulation).output_folder, second_col);
			strcat((*simulation).output_folder, "output_");
			strcpy(hold, (*simulation).name);
			strcat((*simulation).output_folder, hold);
			strcat((*simulation).output_folder, "/");
			// create output folder if it does not exist
			struct stat st = {0};
			if (stat((*simulation).output_folder, &st) == -1) {
				mkdir((*simulation).output_folder, 0700);
			}
		}
		else if (strcmp(first_col, "system_specs") == 0)
		{
			strcpy((*simulation).system_specs, (*simulation).input_folder);
			strcat((*simulation).system_specs, second_col);
		}
		else if (strcmp(first_col, "integrator_specs") == 0)
		{
			strcpy((*simulation).integrator_specs, (*simulation).input_folder);
			strcat((*simulation).integrator_specs, second_col);
		}
		else if (strcmp(first_col, "dev_specs") == 0)
		{
			strcpy((*simulation).dev_specs, (*simulation).input_folder);
			strcat((*simulation).dev_specs, second_col);
			dev_specs_file_received = true;
		}
		else if (strcmp(first_col, "system_file_type") == 0)
		{
			(*simulation).system_file_type = atoi(second_col);
			system_file_type_received = true;
		}
		else if (strcmp(first_col, "number_of_bodies") == 0)
		{
			(*simulation).number_of_bodies = atoi(second_col);
			number_of_bodies_received = true;
		}
	}
	fclose(in);

	/* verify if any of the given files does not exist */
	FILE *in_system = fopen((*simulation).system_specs, "r");
	if (in_system == NULL)
	{
		fprintf(stderr, "Warning: could not read system specs.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	fclose(in_system);
	FILE *in_integrator = fopen((*simulation).integrator_specs, "r");
	if (in_integrator == NULL)
	{
		fprintf(stderr, "Warning: could not read integrator specs.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	fclose(in_integrator);
	if (dev_specs_file_received == true)
	{
		FILE *in_dev = fopen((*simulation).dev_specs, "r");
		if (in_integrator == NULL)
		{
			fprintf(stderr, "Warning: could not read integrator specs.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}
		fclose(in_dev);	
	}

	/* check and set file type */
	if (system_file_type_received == true)
	{
		if ((*simulation).system_file_type != 1 && (*simulation).system_file_type != 2)
		{
			fprintf(stderr, "Error: Type of system file should be either 1 or 2.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(15);
		}
	}
	else
	{
		(*simulation).system_file_type = 2;
	}

	/* check and set number of bodies */
	FILE *in_col = fopen((*simulation).system_specs, "r");
	read = getline(&line, &len, in_col);
	int col_num = count_columns(line);
	fclose(in_col);
	if (number_of_bodies_received == true)
	{
		if (col_num < (*simulation).number_of_bodies)
		{
			fprintf(stderr, "Warning: number of bodies cannot\n");
			fprintf(stderr, "exceed number of columns in system file.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}
		else if ((*simulation).number_of_bodies < 1)
		{
			fprintf(stderr, "Warning: Number of bodies should be greater than 0.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(15);
		}
	}
	else
	{
		(*simulation).number_of_bodies = col_num - 1;
	}

	/* gravitational constant */
	(*simulation).G = 4.0 * M_PI * M_PI; // AU Msun year

	/* simulation units */
	strcpy((*simulation).units, "AU_Msun_year");

	/* read dev file */
	FILE *in3 = fopen((*simulation).dev_specs, "r");
	if	(in3 != NULL)
	{
		while ((read = getline(&line, &len, in3)) != -1)
		{
			sscanf(line, "%s %s",
				first_col, second_col);
			if (strcmp(first_col, "Units") == 0)
			{
				if (strcmp(second_col, "SI") == 0)
				{
					(*simulation).G = 6.6743e-11; // SI
					strcpy((*simulation).units, second_col);
				}
				else if (strcmp(second_col, "G_unity") == 0)
				{
					(*simulation).G = 1.0; // non-dimensional
					strcpy((*simulation).units, second_col);
				}
			}
		}
		fclose(in3);
	}

	/* verification variables for integrator input */
	int 	number_integrator_inputs = 7;
	bool	input_integrator_received[number_integrator_inputs];
	for (int i = 0; i < number_integrator_inputs; i++)
	{
		input_integrator_received[i] = false;
	}

	/* reading integrator specs from user */
	FILE *in2 = fopen((*simulation).integrator_specs, "r");
	if	(in2 == NULL)
	{
		fprintf(stderr, "Warning: could not read integrator specs file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
   	while ((read = getline(&line, &len, in2)) != -1)
	{
		sscanf(line, "%s %lf", first_col, &second_col_double);
		if (strcmp(first_col, "t_init(yr)") == 0)
		{
			(*simulation).t_init = second_col_double;
			input_integrator_received[0] = true;
		}
		else if (strcmp(first_col, "t_trans(yr)") == 0)
		{
			(*simulation).t_trans = second_col_double;
			input_integrator_received[1] = true;
		}
		else if (strcmp(first_col, "t_final(yr)") == 0)
		{
			(*simulation).t_final = second_col_double;
			input_integrator_received[2] = true;
		}
		else if (strcmp(first_col, "t_step(yr)") == 0)
		{
			(*simulation).t_step = second_col_double;
			input_integrator_received[3] = true;
		}
		else if (strcmp(first_col, "eps_abs") == 0)
		{
			(*simulation).eps_abs = second_col_double;
			input_integrator_received[4] = true;
		}
		else if (strcmp(first_col, "eps_rel") == 0)
		{
			(*simulation).eps_rel = second_col_double;
			input_integrator_received[5] = true;
		}
		else if (strcmp(first_col, "data_skip") == 0)
		{
			(*simulation).data_skip = (int) second_col_double;
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
			fprintf(stderr, "from %s.\n", (*simulation).integrator_specs);
			exit(13);
		}
	}

	return 0;
}

int
fill_in_bodies_data	(cltbdy	**bodies,
				 	 const siminf simulation)
{
	if (simulation.system_file_type == 1)
	{
		fprintf(stderr, "Error: simulation.system_file_type == 1 ");
		fprintf(stderr, "not implemented yet.\n");
		exit(22);

		// /* orbital parameters given by user */
		// double 	e = 0.0, a = 0.0;
		// /* state variables given by user*/
		// double 	b0_diag[] = {0.0, 0.0, 0.0};
		// /* non-state variables given by user */
		// double 	omega[3];
		// omega[0] = omega[0];
		// omega[1] = omega[1];
		// omega[2] = omega[2];
		// /* system parameters given by user */
		// int		elements = 0; // number of Voigt elements
		// double  m1 = 0.0, m2 = 0.0;
		// double	I0 = I0, R = R, kf = kf;
		// double	alpha = 0.0, eta = 0.0, alpha_0 = 0.0;
		// /* Voigt elements for Maxwell generalized rheology */
		// double	*alpha_elements = *(&alpha_elements);
		// double	*eta_elements = *(&eta_elements);
		// /* position and velocity */
		// double	tilde_x[3];
		// tilde_x[0] = tilde_x[0];
		// tilde_x[1] = tilde_x[1];
		// tilde_x[2] = tilde_x[2];
		// double 	tilde_x_dot[3];
		// tilde_x_dot[0] = tilde_x_dot[0];
		// tilde_x_dot[1] = tilde_x_dot[1];
		// tilde_x_dot[2] = tilde_x_dot[2];
		// /* deformation settings */
		// // bool	centrifugal = false;
		// // bool	tidal = false;

		// /* auxiliary variables for fscanf */
		// char 	var_name[100];
		// double 	var_value;

		// // Warning for dev
		// fprintf(stderr, "Deformable settings not implemented");
		// fprintf(stderr, " yet for this type of file!\n");
		// exit(29);

		// /* verification variables for system input */
		// int 	number_system_inputs = 17;
		// bool	input_system_received[number_system_inputs];
		// for (int i = 0; i < number_system_inputs; i++)
		// {
		// 	input_system_received[i] = false;
		// }

		// /* reading system specs from user */
		// FILE *in1 = fopen(simulation.system_specs, "r");
		// if	(in1 == NULL)
		// {
		// 	fprintf(stderr, "Warning: could not read system specs file.\n");
		// 	fprintf(stderr, "Exiting the program now.\n");
		// 	exit(13);
		// }
		// while(fscanf(in1, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
		// {
		// 	if (strcmp(var_name, "e") == 0)
		// 	{
		// 		e = var_value;
		// 		input_system_received[0] = true;
		// 	}					
		// 	else if (strcmp(var_name, "a") == 0)
		// 	{
		// 		a = var_value;
		// 		input_system_received[1] = true;
		// 	}
		// 	else if (strcmp(var_name, "m1") == 0)
		// 	{
		// 		m1 = var_value;
		// 		input_system_received[2] = true;
		// 	}
		// 	else if (strcmp(var_name, "m2") == 0)
		// 	{
		// 		m2 = var_value;
		// 		input_system_received[3] = true;
		// 	}
		// 	else if (strcmp(var_name, "I0") == 0)
		// 	{
		// 		I0 = var_value;
		// 		input_system_received[4] = true;
		// 	}
		// 	else if (strcmp(var_name, "R") == 0)
		// 	{
		// 		R = var_value;
		// 		input_system_received[5] = true;
		// 	}
		// 	else if (strcmp(var_name, "kf") == 0)
		// 	{
		// 		kf = var_value;
		// 		input_system_received[6] = true;
		// 	}
		// 	else if (strcmp(var_name, "b0_x") == 0)
		// 	{
		// 		b0_diag[0] = var_value;
		// 		input_system_received[7] = true;
		// 	}
		// 	else if (strcmp(var_name, "b0_y") == 0)
		// 	{
		// 		b0_diag[1] = var_value;
		// 		input_system_received[8] = true;
		// 	}
		// 	else if (strcmp(var_name, "b0_z") == 0)
		// 	{
		// 		b0_diag[2] = var_value;
		// 		input_system_received[9] = true;
		// 	}
		// 	else if (strcmp(var_name, "omega_x") == 0)
		// 	{
		// 		omega[0] = var_value;
		// 		input_system_received[10] = true;
		// 	}
		// 	else if (strcmp(var_name, "omega_y") == 0)
		// 	{
		// 		omega[1] = var_value;
		// 		input_system_received[11] = true;
		// 	}
		// 	else if (strcmp(var_name, "omega_z") == 0)
		// 	{
		// 		omega[2] = var_value;
		// 		input_system_received[12] = true;
		// 	}
		// 	else if (strcmp(var_name, "alpha_0") == 0)
		// 	{
		// 		alpha_0 = var_value;
		// 		input_system_received[13] = true;
		// 	}
		// 	else if (strcmp(var_name, "alpha") == 0)
		// 	{
		// 		alpha = var_value;
		// 		input_system_received[14] = true;
		// 	}
		// 	else if (strcmp(var_name, "eta") == 0)
		// 	{
		// 		eta = var_value;
		// 		input_system_received[15] = true;
		// 	}
		// 	else if (strcmp(var_name, "elements") == 0)
		// 	{
		// 		elements = (int) var_value;
		// 		input_system_received[16] = true;
		// 	}
		// }
		// fclose(in1);

		// /* system input verification */
		// for (int i = 0; i < number_system_inputs; i++)
		// {
		// 	if(input_system_received[i] == false)
		// 	{
		// 		fprintf(stderr, "Error: there is at least one missing input ");
		// 		fprintf(stderr, "from %s.\n", simulation.system_specs);
		// 		exit(13);
		// 	}
		// }
		// if (fabs(alpha) < 1e-14 || fabs(eta) < 1e-14)
		// {
		// 	fprintf(stderr, "Error: nor alpha nor eta should be zero.\n");
		// 	exit(13);
		// }

		// /* setting up Voigt elements */
		// if (elements > 0)
		// {
		// 	alpha_elements 	= (double *) malloc(elements * sizeof(double));
		// 	eta_elements 	= (double *) malloc(elements * sizeof(double));

		// 	char name_element_alpha[20], name_element_eta[20];
		// 	bool element_alpha_found[elements];
		// 	bool element_eta_found[elements];
		// 	for (int i = 0; i < elements; i++)
		// 	{
		// 		element_alpha_found[i] = false;
		// 		element_eta_found[i] = false;
		// 	}

		// 	FILE *in1_elements = fopen(simulation.system_specs, "r");
		// 	while(fscanf(in1_elements, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
		// 	{
		// 		for (int i = 0; i < elements; i++)
		// 		{
		// 			sprintf(name_element_alpha, "alpha_%d", i+1);
		// 			sprintf(name_element_eta, "eta_%d", i+1);
		// 			if (strcmp(var_name, name_element_alpha) == 0)
		// 			{
		// 				alpha_elements[i] = var_value;
		// 				element_alpha_found[i] = true;
		// 			}
		// 			else if (strcmp(var_name, name_element_eta) == 0)
		// 			{
		// 				eta_elements[i] = var_value;
		// 				element_eta_found[i] = true;
		// 			}
		// 		}
		// 	}
		// 	fclose(in1_elements);
		// 	for (int i = 0; i < elements; i++)
		// 	{
		// 		if (element_alpha_found[i] == false || 
		// 			element_eta_found[i] == false)
		// 		{
		// 			fprintf(stderr, "Error: parameters missing for Voigt elements.\n");
		// 			exit(10);
		// 		}
		// 		if (alpha_elements[i] < 1e-13 || eta_elements[i] < 1e-13)
		// 		{
		// 			fprintf(stderr, "Warning: nor alpha nor eta should be zero.\n");
		// 			exit(13);
		// 		}
		// 	}
		// }

		// /* position and velocity at periapsis given by Murray */
		// tilde_x[0] 		= a * (1.0 - e);
		// tilde_x[1] 		= 0.0;
		// tilde_x[2] 		= 0.0;	
		// tilde_x_dot[0]	= 0.0;
		// tilde_x_dot[1] 	= ((2.0 * M_PI) / kepler_period(m1, m2, simulation.G, a)) 
		// 					* a * sqrt((1.0 + e)/(1.0 - e));
		// tilde_x_dot[2] 	= 0.0;

		// /* for testing */
		// // print_vector(tilde_x);
		// // print_vector(tilde_x_dot);
		// // e = calculate_eccentricity(G, m1, m2, tilde_x, tilde_x_dot);
		// // a = calculate_semi_major_axis(G, m1, m2, tilde_x, tilde_x_dot);
		// // printf("e = %e a = %e\n", e, a);
		// // exit(99);

	} // end if (simulation.system_file_type == 1)
	else if (simulation.system_file_type == 2) /* convert input if necessary */
	{
		/* allocate memory for bodies */
		*bodies = (cltbdy *) malloc(simulation.number_of_bodies * sizeof(cltbdy));

		/* auxiliary variables for reading input */
		char 	*line = NULL;
		size_t 	len = 0;
		ssize_t read;

		/* verification variables for input */
		int 	number_par_inputs = 18;
		bool	input_par_received[number_par_inputs];
		for (int i = 0; i < number_par_inputs; i++)
		{
			input_par_received[i] = false;
		}
		/* verification variables for names, orbit and deformable settings */
		bool	input_name_received = false;
		bool	input_keplerian_received = false;
		bool	input_orbit_2body_received = false;
		int 	number_deformable_inputs = 3; 
		bool	input_deformable_received[number_deformable_inputs];
		for (int i = 0; i < number_deformable_inputs; i++)
		{
			input_deformable_received[i] = false;
		}

		/* reading input parameters */
		FILE 	*in1 = fopen(simulation.system_specs, "r");
		while ((read = getline(&line, &len, in1)) != -1)
		{
			const char tok_del[6] = " \t\n";	// token delimiter
			char *token = strtok(line, tok_del);
			if(token == NULL) break;	// in case there is a newline
			if (strcmp(token, "Name") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					strcpy((*bodies)[i].name, token);
				}
				input_name_received = true;
			}
			else if (strcmp(token, "mass(Msun)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].mass = atof(token);
				}
				input_par_received[0] = true;
			}
			else if (strcmp(token, "lod(day)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].lod = atof(token);
				}
				input_par_received[1] = true;
			}
			else if (strcmp(token, "obl(deg)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].obl = atof(token);
				}
				input_par_received[2] = true;
			}
			else if (strcmp(token, "psi(deg)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].psi = atof(token);
				}
				input_par_received[3] = true;
			}
			else if (strcmp(token, "R(km)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].R = atof(token);
				}
				input_par_received[4] = true;
			}
			else if (strcmp(token, "rg") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].rg = atof(token);
				}
				input_par_received[5] = true;
			}
			else if (strcmp(token, "J2") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].J2 = atof(token);
				}
				input_par_received[6] = true;
			}
			else if (strcmp(token, "C22") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].C22 = atof(token);
				}
				input_par_received[7] = true;
			}
			else if (strcmp(token, "lib(deg)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].lib = atof(token);
				}
				input_par_received[8] = true;
			}
			else if (strcmp(token, "kf") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].kf = atof(token);
				}
				input_par_received[9] = true;
			}
			else if (strcmp(token, "Dt(s)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].Dt = atof(token);
				}
				input_par_received[10] = true;
			}
			else if (strcmp(token, "tau(yr)") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].tau = atof(token);
				}
				input_par_received[11] = true;
			}
			else if (strcmp(token, "a(AU)") == 0)
			{
				(*bodies)[0].a = NAN;
				for (int i = 1; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].a = atof(token);
				}
				input_par_received[12] = true;
			}
			else if (strcmp(token, "e") == 0)
			{
				(*bodies)[0].e = NAN;
				for (int i = 1; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].e = atof(token);
				}
				input_par_received[13] = true;
			}
			else if (strcmp(token, "I(deg)") == 0)
			{
				(*bodies)[0].I = NAN;
				for (int i = 1; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].I = atof(token);
				}
				input_par_received[14] = true;
			}
			else if (strcmp(token, "M(deg)") == 0)
			{
				(*bodies)[0].M = NAN;
				for (int i = 1; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].M = atof(token);
				}
				input_par_received[15] = true;
			}
			else if (strcmp(token, "w(deg)") == 0)
			{
				(*bodies)[0].w = NAN;
				for (int i = 1; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].w = atof(token);
				}
				input_par_received[16] = true;
			}
			else if (strcmp(token, "OMEGA(deg)") == 0)
			{
				(*bodies)[0].Omega = NAN;
				for (int i = 1; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					(*bodies)[i].Omega = atof(token);
				}
				input_par_received[17] = true;
			}
			else if (strcmp(token, "keplerian") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					if (strcmp(token, "yes") == 0)
					{
						(*bodies)[i].keplerian = true;
					}
					else if (strcmp(token, "no") == 0)
					{
						(*bodies)[i].keplerian = false;
					}
					else
					{
						printf("%s\n", token);
						fprintf(stderr, "Please provide yes or no ");
						fprintf(stderr, "for keplerian variable\n");
						exit(14);
					}
				}
				input_keplerian_received = true;
			}
			else if (strcmp(token, "orbit_2body") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					if (strcmp(token, "yes") == 0)
					{
						(*bodies)[i].orbit_2body = true;
					}
					else if (strcmp(token, "no") == 0)
					{
						(*bodies)[i].orbit_2body = false;
					}
					else
					{
						printf("%s\n", token);
						fprintf(stderr, "Please provide yes or no ");
						fprintf(stderr, "for orbit_2body variable\n");
						exit(14);
					}
				}
				input_orbit_2body_received = true;
			}
			else if (strcmp(token, "centrifugal") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					if (strcmp(token, "yes") == 0)
					{
						(*bodies)[i].centrifugal = true;
					}
					else if (strcmp(token, "no") == 0)
					{
						(*bodies)[i].centrifugal = false;
					}
					else
					{
						printf("%s\n", token);
						fprintf(stderr, "Please provide yes or no ");
						fprintf(stderr, "for centrifugal variable\n");
						exit(14);
					}
				}
				input_deformable_received[0] = true;
			}
			else if (strcmp(token, "tidal") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					if (strcmp(token, "yes") == 0)
					{
						(*bodies)[i].tidal = true;
					}
					else if (strcmp(token, "no") == 0)
					{
						(*bodies)[i].tidal = false;
					}
					else
					{
						fprintf(stderr, "Please provide yes or no ");
						fprintf(stderr, "for tidal variable\n");
						exit(14);
					}
				}
				input_deformable_received[1] = true;
			}
			else if (strcmp(token, "point_mass") == 0)
			{
				for (int i = 0; i < simulation.number_of_bodies; i++)
				{
					token = strtok(NULL, tok_del);
					if (strcmp(token, "yes") == 0)
					{
						(*bodies)[i].point_mass = true;
					}
					else if (strcmp(token, "no") == 0)
					{
						(*bodies)[i].point_mass = false;
					}
					else
					{
						fprintf(stderr, "Please provide yes or no ");
						fprintf(stderr, "for point mass variable\n");
						exit(14);
					}
				}
				input_deformable_received[2] = true;
			}
		}
		fclose(in1);

		/* parameter input verification */
		for (int i = 0; i < number_par_inputs; i++)
		{
			if(input_par_received[i] == false)
			{
				fprintf(stderr, "Error: there is at least one missing input ");
				fprintf(stderr, "from %s.\n", simulation.system_specs);
				fprintf(stderr, "Exiting the program now.\n");
				// fprintf(stderr, "Missing input number %d\n", i); // for testing
				exit(14);
			}
		}
		/* name  input verification */
		if (input_name_received == false)
		{
			fprintf(stderr, "Error: could not read body names ");
			fprintf(stderr, "from %s.\n", simulation.system_specs);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
		/* orbital variable input verification */
		if (input_keplerian_received == false)
		{
			fprintf(stderr, "Error: missing keplerian status ");
			fprintf(stderr, "from %s.\n", simulation.system_specs);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
		if (input_orbit_2body_received == false)
		{
			fprintf(stderr, "Error: missing orbit 2 body status ");
			fprintf(stderr, "from %s.\n", simulation.system_specs);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
		/* deformable variables input verification */
		for (int i = 0; i < number_deformable_inputs; i++)
		{
			if(input_deformable_received[i] == false)
			{
				fprintf(stderr, "Error: there is at least one missing input ");
				fprintf(stderr, "for the deformation variables ");
				fprintf(stderr, "from %s.\n", simulation.system_specs);
				fprintf(stderr, "Exiting the program now.\n");
				exit(14);
			}
		}

		/* converting units */

		double deg_to_rad = M_PI / 180.0;
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].obl *= deg_to_rad;
			(*bodies)[i].psi *= deg_to_rad;
			(*bodies)[i].lib *= deg_to_rad;
			(*bodies)[i].I *= deg_to_rad;
			(*bodies)[i].M *= deg_to_rad;
			(*bodies)[i].w *= deg_to_rad;
			(*bodies)[i].Omega *= deg_to_rad;
		}

		if (strcmp(simulation.units, "SI") == 0)
		{
			/* conversion units to SI */
			double Msun_to_kg = 1988500.0e24;
			double day_to_s = 24.0 * 60.0 * 60.0;
			double km_to_m = 1e3;
			double year_to_s = 365.25 * day_to_s;
			double AU_to_m = 1.495978707e11;

			/* variables in SI*/
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				(*bodies)[i].mass *= Msun_to_kg;
				(*bodies)[i].R *= km_to_m;
				(*bodies)[i].lod *= day_to_s;
				(*bodies)[i].tau *= year_to_s;
				(*bodies)[i].a *= AU_to_m;
			}
		}
		else
		{
			/* conversion units to AU Msun year */
			double km_to_AU = 1.0 / 1.495978707e8;
			double day_to_year = 1.0 / 365.25;
			double s_to_year = day_to_year / (24.0 * 60.0 * 60.0);

			/* variables in AU Msun year*/
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				(*bodies)[i].R *= km_to_AU;
				(*bodies)[i].lod *= day_to_year;
				(*bodies)[i].Dt *= s_to_year;
			}
		}
		
		/* calculating state variables */
		
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			double	m = (*bodies)[i].mass;
			double	R = (*bodies)[i].R;
			double  Td = (*bodies)[i].lod;

			double	rg = (*bodies)[i].rg;
			double	J2 = (*bodies)[i].J2;

			double  theta = (*bodies)[i].obl;
			double	psi = (*bodies)[i].psi;
			double	phi = (*bodies)[i].lib;

			double	a = (*bodies)[i].a;
			double	e = (*bodies)[i].e;
			double	I = (*bodies)[i].I;
			double	M = (*bodies)[i].M;
			double	w = (*bodies)[i].w;
			double	Omega = (*bodies)[i].Omega;

			double	kf = (*bodies)[i].kf;
			double	Dt = (*bodies)[i].Dt;
			double	tau = (*bodies)[i].tau;
			
			/* 1st set of variables - x and x_dot */
			double qr_1_I[4];
			double qr_3_Omega[4];

			if (i == 0)
			{
				null_vector((*bodies)[i].x);
				null_vector((*bodies)[i].x_dot);
				identity_quaternion(qr_1_I);		// for 2nd set of variables
				identity_quaternion(qr_3_Omega);	// for 2nd set of variables
			}
			else
			{
				double T = 
				kepler_period((*bodies)[0].mass, m, simulation.G, a);
				double n = (2.0 * M_PI) / T;

				double E = kepler_equation(e, M);
				double r = a * (1.0 - e * cos(E));
				double f = atan2(sqrt(1.0 - e * e) * sin(E), cos(E) - e);

				double position_in_plane[] 
					= {r * cos(f), r * sin(f), 0.0};
				double velocity_in_plane[] 
					= {-1.0 * n * a / sqrt(1.0 - e * e) * sin(f), 
						n * a / sqrt(1.0 - e * e) * (e + cos(f)), 
						0.0};

				double qr_3_w[9];
				rotation_quaternion_z(qr_3_w, w);
				rotation_quaternion_x(qr_1_I, I);
				rotation_quaternion_z(qr_3_Omega, Omega);

				double full_rotation_orbit_quaternion[4];
				quaternion_times_quaternion(full_rotation_orbit_quaternion,
					qr_3_Omega, qr_1_I);
				quaternion_times_quaternion(full_rotation_orbit_quaternion,
					full_rotation_orbit_quaternion, qr_3_w);

				rotate_vector_with_quaternion((*bodies)[i].x, 
					full_rotation_orbit_quaternion, position_in_plane);
				rotate_vector_with_quaternion((*bodies)[i].x_dot, 
					full_rotation_orbit_quaternion, velocity_in_plane);
			}

			/* 2nd set of variables - omega, body frame variables and b0_diag */
			// b0_diag not implemented yet. 
			// at the moment, it is being dealt with by the main
			double qr_3_psi[4];
			rotation_quaternion_z(qr_3_psi, psi);
			double qr_1_theta[4];
			rotation_quaternion_x(qr_1_theta, theta);
			double qr_3_phi[4];
			rotation_quaternion_z(qr_3_phi, phi);

			double full_rotation_body_quaternion[4];
			quaternion_times_quaternion(full_rotation_body_quaternion,
				qr_3_Omega, qr_1_I);
			quaternion_times_quaternion(full_rotation_body_quaternion,
				full_rotation_body_quaternion, qr_3_psi);
			quaternion_times_quaternion(full_rotation_body_quaternion,
				full_rotation_body_quaternion, qr_1_theta);
			quaternion_times_quaternion(full_rotation_body_quaternion,
				full_rotation_body_quaternion, qr_3_phi);

			copy_quaternion((*bodies)[i].q, full_rotation_body_quaternion);

			double omega_on_body[] = {0.0, 0.0, 0.0};
			double omega_direction_on_body[] = {0.0, 0.0, 1.0}; // strong assumption
			scale_vector(omega_on_body, 2.0 * M_PI / Td, omega_direction_on_body);
			rotate_vector_with_quaternion((*bodies)[i].omega,
				full_rotation_body_quaternion, omega_on_body);

			rotation_matrix_from_quaternion((*bodies)[i].Y, 
				full_rotation_body_quaternion);

			/* 3rd set of variables - I0 */
			(*bodies)[i].I0 = (3.0 * rg - 2.0 * J2) * m * R * R / 3.0;

			/* 4th set of variables - gamma, alpha and eta */
			(*bodies)[i].gamma = parameter_gamma(simulation.G, (*bodies)[i].I0, R, kf);
			(*bodies)[i].alpha = (*bodies)[i].gamma * Dt / (tau - Dt);
			(*bodies)[i].eta = (*bodies)[i].gamma * Dt;

			/* Kelvin-Voigt elements */
			(*bodies)[i].elements = 0;

			if ((*bodies)[i].elements > 0)
			{
				(*bodies)[i].alpha_elements = (double *) malloc((*bodies)[i].elements * sizeof(double));
				(*bodies)[i].eta_elements = (double *) malloc((*bodies)[i].elements * sizeof(double));
				
				for (int j = 0; j < (*bodies)[i].elements; j++)
				{
					(*bodies)[i].alpha_elements[j] = 1.0;
					(*bodies)[i].eta_elements[j] = 1.0;
				}
			}

		} // end loop over bodies

		/* for testing */
		// for (int i = 0; i < number_of_bodies; i++)
		// 	print_CelestialBody((*bodies)[i]);
		// exit(99);

	} // end else if (simulation.system_file_type == 2)

	return 0;
}

int
write_simulation_overview	(const int time_spent_in_seconds,
							 const siminf simulation)
{
	// creates info file
	char filename[150];
	strcpy(filename, simulation.output_folder);
	strcat(filename, "results_");
	strcat(filename, simulation.name);
	strcat(filename, "_");
	strcat(filename, "sim_info");
	strcat(filename, ".dat");

	// simulation time
	int sec, min, hr, days;
	days = time_spent_in_seconds/(24*3600);
	hr	= (time_spent_in_seconds - 24*3600*days) / 3600;
	min = (time_spent_in_seconds - 24*3600*days - 3600*hr) / 60;
	sec = (time_spent_in_seconds - 24*3600*days - 3600*hr - 60*min) / 1;
	FILE *out_sim_info;
	out_sim_info = fopen(filename, "w");
	fprintf(out_sim_info, "Time spent on simulation:");
	fprintf(out_sim_info, " %d days %d hours %d minutes %d seconds.\n\n",
		days, hr, min, sec);
	fclose(out_sim_info);

	// copy input files info
	FILE *in0_to_copy = fopen(simulation.main_input, "r");
	FILE *in0_copy = fopen(filename, "a");
	char ch0 = fgetc(in0_to_copy);
    while(ch0 != EOF)
    {
        fputc(ch0, in0_copy);
        ch0 = fgetc(in0_to_copy);
    }
	fprintf(in0_copy, "\n\n");
	fclose(in0_copy);
	fclose(in0_to_copy);
	FILE *in1_to_copy = fopen(simulation.system_specs, "r");
	FILE *in1_copy = fopen(filename, "a");
	char ch1 = fgetc(in1_to_copy);
    while(ch1 != EOF)
    {
        fputc(ch1, in1_copy);
        ch1 = fgetc(in1_to_copy);
    }
	fprintf(in1_copy, "\n\n");
	fclose(in1_copy);
	fclose(in1_to_copy);
	FILE *in2_to_copy = fopen(simulation.integrator_specs, "r");
	FILE *in2_copy = fopen(filename, "a");
	char ch2 = fgetc(in2_to_copy);
    while(ch2 != EOF)
    {
        fputc(ch2, in2_copy);
        ch2 = fgetc(in2_to_copy);
    }
	fprintf(in2_copy, "\n\n");
	fclose(in2_copy);
	fclose(in2_to_copy);
	FILE *in3_to_copy = fopen(simulation.dev_specs, "r");
	if (in3_to_copy != NULL)
	{
		FILE *in3_copy = fopen(filename, "a");
		char ch3 = fgetc(in3_to_copy);
		while(ch3 != EOF)
		{
			fputc(ch3, in3_copy);
			ch3 = fgetc(in3_to_copy);
		}
		fprintf(in3_copy, "\n\n");
		fclose(in3_copy);
		fclose(in3_to_copy);
	}

	return 0;
}

int
output_to_orbit(cltbdy *bodies,
			 	const siminf simulation)
{
	/* auxiliary variables */
    char 	*line = NULL;
    size_t 	len = 0;
    ssize_t	read;

	/* input filename of first body */
	char filename_read_first_body[150];
	strcpy(filename_read_first_body, simulation.output_folder);
	strcat(filename_read_first_body, "results_");
	strcat(filename_read_first_body, simulation.name);
	strcat(filename_read_first_body, "_");
	strcat(filename_read_first_body, bodies[0].name);
	strcat(filename_read_first_body, ".dat");

	/* loop over bodies except first one */
	for (int i = 1; i < simulation.number_of_bodies; i++)
	{
		/* open input file for first body */
		FILE *in_state_first_body = fopen(filename_read_first_body, "r");
		if (in_state_first_body == NULL)
		{
			fprintf(stderr, "Warning: could not read state evolution file.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}

		/* auxiliary variables */
		char	*line_first_body = NULL;
		size_t 	len_first_body = 0;
		ssize_t	read_first_body;

		/* input files of remaining bodies */
		char filename_read[150];
		strcpy(filename_read, simulation.output_folder);
		strcat(filename_read, "results_");
		strcat(filename_read, simulation.name);
		strcat(filename_read, "_");
		strcat(filename_read, bodies[i].name);
		strcat(filename_read, ".dat");
		FILE *in_state = fopen(filename_read, "r");
		if (in_state == NULL)
		{
			fprintf(stderr, "Warning: could not read state evolution file.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}

		/* output files for orbital elements */
		char filename_orbital[150];
		FILE *out_orbital;

		/* names and opens output files for each body */
		strcpy(filename_orbital, simulation.output_folder);
		strcat(filename_orbital, "results_");
		strcat(filename_orbital, simulation.name);
		strcat(filename_orbital, "_");
		strcat(filename_orbital, bodies[i].name);
		strcat(filename_orbital, "_orbital_elements");
		strcat(filename_orbital, ".dat");
		out_orbital = fopen(filename_orbital, "w");
		/* writes output headers */
		fprintf(out_orbital, "time(yr)");
		fprintf(out_orbital, " a(AU)");
		fprintf(out_orbital, " e");
		fprintf(out_orbital, " I(°)");
		fprintf(out_orbital, " M(°)");
		fprintf(out_orbital, " w(°)");
		fprintf(out_orbital, " Omega(°)");
		fprintf(out_orbital, "\n");

		// convertion units
		double rad_to_deg = 180.0 / M_PI;

		/* reading input and printing output */
		read = getline(&line, &len, in_state); // discards first line
		read_first_body = 
			getline(&line_first_body, &len_first_body, in_state_first_body);
		while (true)
		{
			/* reading data */
			read = getline(&line, &len, in_state);
			read_first_body = 
				getline(&line_first_body, &len_first_body, in_state_first_body);
			if (read == -1 || read_first_body == -1) break;

			/* tokenizing data with reentrant version of strtok */
			const char tok_del[6] = " \t\n";		// token delimiter
			char *token;
			char *line_track = line;
			char *token_first_body;
			char *line_first_body_track = line_first_body;
			token = strtok_r(line_track, tok_del, &line_track);
			token_first_body = strtok_r(line_first_body_track, tok_del, &line_first_body_track);
			if(token == NULL || token_first_body == NULL) break; // in case there is a newline
					
			/* filling in data */
			double t = atof(token_first_body);
			for (int j = 0; j < 3; j++)
			{
				token_first_body = strtok_r(line_first_body_track, tok_del, &line_first_body_track);
				bodies[0].x[j] = atof(token_first_body);
				token = strtok_r(line_track, tok_del, &line_track);
				bodies[i].x[j] = atof(token);
			}
			for (int j = 0; j < 3; j++)
			{
				token_first_body = strtok_r(line_first_body_track, tok_del, &line_first_body_track);
				bodies[0].x_dot[j] = atof(token_first_body);
				token = strtok_r(line_track, tok_del, &line_track);
				bodies[i].x_dot[j] = atof(token);
			}
				
			/* calculate orbital elements from state vector */
			calculate_orbital_elements(&bodies[i], bodies[0], simulation.G);

			/* time */
			fprintf (out_orbital, "%.15e", t);

			/* semimajor axis */
			fprintf (out_orbital, " %.15e", bodies[i].a);

			/* orbital eccentricity */
			fprintf (out_orbital, " %.15e", bodies[i].e);

			/* inclination */
			fprintf (out_orbital, " %.15e", bodies[i].I * rad_to_deg);

			/* mean anomaly */
			fprintf (out_orbital, " %.15e", bodies[i].M * rad_to_deg);

			/* argument of periapsis */
			fprintf (out_orbital, " %.15e", bodies[i].w * rad_to_deg);

			/* longitude of the ascending node */
			fprintf (out_orbital, " %.15e", bodies[i].Omega * rad_to_deg);

			/* next line */
			fprintf(out_orbital, "\n");		
		}

		/* close all files */				
		fclose(in_state);
		fclose(out_orbital);
		fclose(in_state_first_body);

	} // end loop over bodies except first one

	return 0;
}

int
output_to_spin	(cltbdy *bodies,
				 const siminf simulation)
{
	/* auxiliary variables */
    char 	*line = NULL;
    size_t 	len = 0;
    ssize_t	read;

	/* input filename of first body */
	char filename_read_first_body[150];
	strcpy(filename_read_first_body, simulation.output_folder);
	strcat(filename_read_first_body, "results_");
	strcat(filename_read_first_body, simulation.name);
	strcat(filename_read_first_body, "_");
	strcat(filename_read_first_body, bodies[0].name);
	strcat(filename_read_first_body, ".dat");

	/* loop over all bodies */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if(bodies[i].point_mass == false)
		{
			// obs.: it reads the first body file twice for i = 0

			/* open input file for first body */
			FILE *in_state_first_body = fopen(filename_read_first_body, "r");
			if (in_state_first_body == NULL)
			{
				fprintf(stderr, "Warning: could not read state evolution file.\n");
				fprintf(stderr, "Exiting the program now.\n");
				exit(13);
			}

			/* auxiliary variables */
			char	*line_first_body = NULL;
			size_t 	len_first_body = 0;
			ssize_t	read_first_body;

			/* input files of remaining bodies */
			char filename_read[150];
			strcpy(filename_read, simulation.output_folder);
			strcat(filename_read, "results_");
			strcat(filename_read, simulation.name);
			strcat(filename_read, "_");
			strcat(filename_read, bodies[i].name);
			strcat(filename_read, ".dat");
			FILE *in_state = fopen(filename_read, "r");
			if (in_state == NULL)
			{
				fprintf(stderr, "Warning: could not read state evolution file.\n");
				fprintf(stderr, "Exiting the program now.\n");
				exit(13);
			}

			/* output files for orientation variables */
			char filename_orientation[150];
			FILE *out_orientation;

			/* names and opens output files for each body */
			strcpy(filename_orientation, simulation.output_folder);
			strcat(filename_orientation, "results_");
			strcat(filename_orientation, simulation.name);
			strcat(filename_orientation, "_");
			strcat(filename_orientation, bodies[i].name);
			strcat(filename_orientation, "_orientation");
			strcat(filename_orientation, ".dat");
			out_orientation = fopen(filename_orientation, "w");
			/* writes output headers */
			fprintf(out_orientation, "time(yr)");
			fprintf(out_orientation, " |omega|");
			fprintf(out_orientation, " |l|");
			fprintf(out_orientation, " |b|");
			fprintf(out_orientation, " obl(°)");
			fprintf(out_orientation, "\n");

			// convertion units
			double rad_to_deg = 180.0 / M_PI;

			/* reading input and printing output */
			read = getline(&line, &len, in_state); // discards first line
			read_first_body = 
				getline(&line_first_body, &len_first_body, in_state_first_body);
			while (true)
			{
				/* reading data */
				read = getline(&line, &len, in_state);
				read_first_body = 
					getline(&line_first_body, &len_first_body, in_state_first_body);
				if (read == -1 || read_first_body == -1) break;

				/* tokenizing data with reentrant version of strtok */
				const char tok_del[6] = " \t\n";		// token delimiter
				char *token;
				char *line_track = line;
				char *token_first_body;
				char *line_first_body_track = line_first_body;
				token = strtok_r(line_track, tok_del, &line_track);
				token_first_body = strtok_r(line_first_body_track, tok_del, &line_first_body_track);
				if(token == NULL || token_first_body == NULL) break; // in case there is a newline
						
				/* filling in data */
				double t = atof(token_first_body);
				for (int j = 0; j < 3; j++)
				{
					token_first_body = strtok_r(line_first_body_track, tok_del, &line_first_body_track);
					bodies[0].x[j] = atof(token_first_body);
					token = strtok_r(line_track, tok_del, &line_track);
					bodies[i].x[j] = atof(token);
				}
				for (int j = 0; j < 3; j++)
				{
					token_first_body = strtok_r(line_first_body_track, tok_del, &line_first_body_track);
					bodies[0].x_dot[j] = atof(token_first_body);
					token = strtok_r(line_track, tok_del, &line_track);
					bodies[i].x_dot[j] = atof(token);
				}
				for (int j = 0; j < 3; j++)
				{
					token = strtok_r(line_track, tok_del, &line_track);
					bodies[i].omega[j] = atof(token);
				}
				for (int j = 0; j < 3; j++)
				{
					token = strtok_r(line_track, tok_del, &line_track);
					bodies[i].l[j] = atof(token);
				}
				token = strtok_r(line_track, tok_del, &line_track);
				double b_norm = atof(token);
					
				double relative_x[3];
				linear_combination_vector(relative_x,
					1.0, bodies[i].x,
					-1.0, bodies[0].x);
				double relative_x_dot[3];
				linear_combination_vector(relative_x_dot,
					1.0, bodies[i].x_dot,
					-1.0, bodies[0].x_dot);

				/* time */
				fprintf (out_orientation, "%.15e", t);

				/* norm of angular velocity vector */
				fprintf (out_orientation, " %.15e", norm_vector(bodies[i].omega));

				/* norm of angular momentum vector */
				fprintf (out_orientation, " %.15e", norm_vector(bodies[i].l));

				/* norm of deformation matrix */
				fprintf (out_orientation, " %.15e", b_norm);

				/* obliquity */
				double obliquity;
				if (i == 0)
				{
					obliquity = calculate_obliquity_free_body(bodies[i].omega);
				}
				else
				{
					obliquity = calculate_obliquity_on_orbit(relative_x, relative_x_dot, bodies[i].omega);
				}	
				fprintf (out_orientation, " %.15e", obliquity * rad_to_deg);

				/* next line */
				fprintf(out_orientation, "\n");		
			}

			/* close all files */				
			fclose(in_state);
			fclose(out_orientation);
			fclose(in_state_first_body);
		}

	} // end loop over all bodies

	return 0;
}

int
plot_output_comma_orbit_and_spin(const cltbdy *bodies,
								 const siminf simulation)
{
	/* suppresses warning messages by moving
		* stderr to /dev/null
		* should be used with care
	*/
	if(freopen("/dev/null", "w", stderr) == NULL)
	{
		fprintf(stderr, "Error: stderr could not be moved to\n");
		fprintf(stderr, "/dev/null in plot command.\n");
		exit(22);
	}
	/* create figures dir */
	char plot_folder[150];
	strcpy(plot_folder, simulation.output_folder);
	strcat(plot_folder, "/figures/");
	struct stat st_plot = {0};
	if (stat(plot_folder, &st_plot) == -1) {
		mkdir(plot_folder, 0700);
	}
	char filename_get[150];
	char filename_plot[150];
	char filename_get_first_body[150];
	/* get output filenames of first body */
	strcpy(filename_get_first_body, simulation.output_folder);
	strcat(filename_get_first_body, "results_");
	strcat(filename_get_first_body, simulation.name);
	strcat(filename_get_first_body, "_");
	strcat(filename_get_first_body, bodies[0].name);
	strcat(filename_get_first_body, ".dat");
	FILE *gnuplotPipe;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		/* get output filenames of each body */
		strcpy(filename_get, simulation.output_folder);
		strcat(filename_get, "results_");
		strcat(filename_get, simulation.name);
		strcat(filename_get, "_");
		strcat(filename_get, bodies[i].name);
		strcat(filename_get, ".dat");
		/* plot outputs */
		/* orbit */
		if (i > 0)
		{
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_orbit");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"x(AU)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"y(AU)\"\n");
			fprintf(gnuplotPipe, "set zlabel \"z(AU)\" rotate parallel\n");
			fprintf(gnuplotPipe, "set title \"%s's orbit\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "splot \'%s\' u 2:3:4 w d, \'%s\' u 2:3:4 w p pt 7 ps 3", 
				filename_get, filename_get_first_body);
			pclose(gnuplotPipe);
		}
		/* reads output files of orbital elements of each body */
		strcpy(filename_get, simulation.output_folder);
		strcat(filename_get, "results_");
		strcat(filename_get, simulation.name);
		strcat(filename_get, "_");
		strcat(filename_get, bodies[i].name);
		strcat(filename_get, "_orbital_elements");
		strcat(filename_get, ".dat");
		/* plot outputs */
		if (i > 0)
		{
			/* semi-major axis */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_semi_major_axis");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"a(AU)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's semi-major axis\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:2 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* eccentricity */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_eccentricity");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"e\"\n");
			fprintf(gnuplotPipe, "set title \"%s's eccentricity\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:3 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* inclination */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_inclination");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"I(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's inclination\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:4 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* mean anomaly */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_mean_anomaly");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"M(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's mean anomaly\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:5 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* argument of periapsis */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_argument_of_periapsis");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"w(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's argument of periapsis\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:6 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* longitude of the asceding node */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_longitude_of_the_asceding_node");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"{/Symbol W}(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's longitude of the ascending node\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:7 w l lw 3", filename_get);
			pclose(gnuplotPipe);
		} // end if i > 0
		/* reads output files of orientation variables of each body */
		if (bodies[i].point_mass == false)
		{
			strcpy(filename_get, simulation.output_folder);
			strcat(filename_get, "results_");
			strcat(filename_get, simulation.name);
			strcat(filename_get, "_");
			strcat(filename_get, bodies[i].name);
			strcat(filename_get, "_orientation");
			strcat(filename_get, ".dat");
			/* plot outputs */
			/* angular velocity */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_angular_velocity");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"|angular velocity|(rad/yr)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's angular velocity\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:2 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* angular momentum */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_angular_momentum");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"|l|\"\n");
			fprintf(gnuplotPipe, "set title \"%s's angular momentum\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:3 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* deformation */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_deformation");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"|b|\"\n");
			fprintf(gnuplotPipe, "set title \"%s's deformation\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:4 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* obliquity */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_obliquity");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"Obliquity(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's obliquity\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:5 w l lw 3", filename_get);
			pclose(gnuplotPipe);
		}
	} // end loop over bodies
	/* reads output file of full system */
	strcpy(filename_get, simulation.output_folder);
	strcat(filename_get, "results_");
	strcat(filename_get, simulation.name);
	strcat(filename_get, "_");
	strcat(filename_get, "full_system");
	strcat(filename_get, ".dat");
	/* plot outputs */
	/* total angular momentum */
	strcpy(filename_plot, plot_folder);
	strcat(filename_plot, "figure_");
	strcat(filename_plot, simulation.name);
	strcat(filename_plot, "_");
	strcat(filename_plot, "full_system");
	strcat(filename_plot, "_total_angular_momentum");
	strcat(filename_plot, ".png");
	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
	fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
	fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
	fprintf(gnuplotPipe, "set ylabel \"|Total angular momentum|\"\n");
	fprintf(gnuplotPipe, "set title \"System's total angular momentum\"\n");
	fprintf(gnuplotPipe, "plot \'%s\' u 1:5 w l lw 3", filename_get);
	pclose(gnuplotPipe);
	
	return 0;
}
