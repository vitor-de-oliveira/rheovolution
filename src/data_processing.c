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
print_SimulationInfo(siminf simulation)
{
	printf("name = %s\n", simulation.name);

	printf("main input = %s\n", simulation.main_input);

	printf("input folder = %s\n", simulation.input_folder);
	printf("output folder = %s\n", simulation.output_folder);
	printf("system specs = %s\n", simulation.system_specs);
	printf("integration specs = %s\n", simulation.integration_specs);
	printf("dev specs = %s\n", simulation.dev_specs);

	printf("G = %1.10e\n", simulation.G);
	printf("units = %s\n", simulation.units);
	printf("rheology model = %s\n", simulation.rheology_model);
	printf("number of bodies = %d\n", simulation.number_of_bodies);
	printf("keplerian motion = %d\n", simulation.keplerian_motion);
	printf("two bodies approx = %d\n", simulation.two_bodies_aprox);

	printf("t init = %1.10e\n", simulation.t_init);
	printf("t trans = %1.10e\n", simulation.t_trans);
	printf("t final = %1.10e\n", simulation.t_final);
	printf("t step = %1.10e\n", simulation.t_step);
	printf("t step received = %d\n", simulation.t_step_received);

	printf("t step init = %1.10e\n", simulation.t_step_init);
	printf("t step min = %1.10e\n", simulation.t_step_min);
	printf("eps abs = %1.10e\n", simulation.error_abs);
	printf("eps rel = %1.10e\n", simulation.error_rel);

	printf("output size = %1.5e\n", simulation.output_size);
	printf("data skip = %d\n", simulation.data_skip);

	return 0;
}

int
parse_input(siminf *simulation,
			const char input_file[])
{
	/* auxiliary variables */
    char 	*line = NULL;
    size_t 	len = 0;
    ssize_t	read;
	char 	first_col[300];
	char 	second_col[300];
	char	hold[300];
	bool	rheology_model_received = false;
	bool	dev_specs_file_received = false;
	bool	number_of_bodies_received = false;
	bool	keplerian_motion_received = false;
	bool	two_bodies_aprox_received = false;

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
		strcpy(simulation->main_input, input_file);
	}
   	while ((read = getline(&line, &len, in)) != -1)
	{
		sscanf(line, "%s %s", first_col, second_col);
		if (strcmp(first_col, "name") == 0)
		{
			strcpy(simulation->name, second_col);
		}
		else if (strcmp(first_col, "input_folder") == 0)
		{
			strcpy(simulation->input_folder, second_col);
		}
		else if (strcmp(first_col, "output_folder") == 0)
		{
			strcpy(simulation->output_folder, second_col);
			strcat(simulation->output_folder, "output_");
			strcpy(hold, simulation->name);
			strcat(simulation->output_folder, hold);
			strcat(simulation->output_folder, "/");
			// create output folder if it does not exist
			struct stat st = {0};
			if (stat(simulation->output_folder, &st) == -1) {
				mkdir(simulation->output_folder, 0700);
			}
		}
		else if (strcmp(first_col, "system_specs") == 0)
		{
			strcpy(simulation->system_specs, simulation->input_folder);
			strcat(simulation->system_specs, second_col);
		}
		else if (strcmp(first_col, "integration_specs") == 0)
		{
			strcpy(simulation->integration_specs, simulation->input_folder);
			strcat(simulation->integration_specs, second_col);
		}
		else if (strcmp(first_col, "dev_specs") == 0)
		{
			strcpy(simulation->dev_specs, simulation->input_folder);
			strcat(simulation->dev_specs, second_col);
			dev_specs_file_received = true;
		}
		else if (strcmp(first_col, "rheology_model") == 0)
		{
			if (strcmp(second_col, "Maxwell") == 0 ||
				strcmp(second_col, "gen_Voigt") == 0)
			{
				strcpy(simulation->rheology_model, second_col);
				rheology_model_received = true;
			}
			else
			{
				printf("%s\n", second_col);
				fprintf(stderr, "Please provide Maxwell or gen_Voigt ");
				fprintf(stderr, "for the rheology model\n");
				exit(14);
			}
			rheology_model_received = true;
		}
		else if (strcmp(first_col, "number_of_bodies") == 0)
		{
			simulation->number_of_bodies = atoi(second_col);
			number_of_bodies_received = true;
		}
		else if (strcmp(first_col, "keplerian_motion") == 0)
		{
			if (strcmp(second_col, "yes") == 0)
			{
				simulation->keplerian_motion = true;
			}
			else if (strcmp(second_col, "no") == 0)
			{
				simulation->keplerian_motion = false;
			}
			else
			{
				printf("%s\n", second_col);
				fprintf(stderr, "Please provide yes or no ");
				fprintf(stderr, "for keplerian motion\n");
				exit(14);
			}
			keplerian_motion_received = true;
		}
		else if (strcmp(first_col, "two_bodies_aprox") == 0)
		{
			if (strcmp(second_col, "yes") == 0)
			{
				simulation->two_bodies_aprox = true;
			}
			else if (strcmp(second_col, "no") == 0)
			{
				simulation->two_bodies_aprox = false;
			}
			else
			{
				printf("%s\n", second_col);
				fprintf(stderr, "Please provide yes or no ");
				fprintf(stderr, "for 2 bodies approximation\n");
				exit(14);
			}
			two_bodies_aprox_received = true;
		}
	}
	fclose(in);

	/* verify if any of the given files does not exist */
	FILE *in_system = fopen(simulation->system_specs, "r");
	if (in_system == NULL)
	{
		fprintf(stderr, "Warning: could not read system specs.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	fclose(in_system);
	FILE *in_integration = fopen(simulation->integration_specs, "r");
	if (in_integration == NULL)
	{
		fprintf(stderr, "Warning: could not read integration specs.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
	fclose(in_integration);
	if (dev_specs_file_received == true)
	{
		FILE *in_dev = fopen(simulation->dev_specs, "r");
		if (in_dev == NULL)
		{
			fprintf(stderr, "Warning: could not read dev specs.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}
		fclose(in_dev);	
	}
	else
	{
		strcpy(simulation->dev_specs, "Not received.");
	}

	/* check and set number of bodies */
	FILE *in_col = fopen(simulation->system_specs, "r");
	read = getline(&line, &len, in_col);
	int col_num = count_columns(line);
	fclose(in_col);
	if (number_of_bodies_received == true)
	{
		if (col_num < simulation->number_of_bodies)
		{
			fprintf(stderr, "Warning: number of bodies cannot\n");
			fprintf(stderr, "exceed number of columns in system file.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}
		else if (simulation->number_of_bodies < 1)
		{
			fprintf(stderr, "Warning: Number of bodies should be greater than 0.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(15);
		}
	}
	else
	{
		simulation->number_of_bodies = col_num - 1;
	}

	/* check rheology model */
	if (rheology_model_received == false)
	{
		strcpy(simulation->rheology_model, "gen_Voigt");
	}

	/* check and set keplerian motion */
	if (keplerian_motion_received == false)
	{
		simulation->keplerian_motion = false;
	}

	/* check and set two bodies approximation */
	if (two_bodies_aprox_received == false)
	{
		simulation->two_bodies_aprox = false;
	}

	/* presetting some values */
	simulation->t_step_received = false;
	simulation->t_step_min = NAN;
	simulation->error_abs = NAN;
	simulation->error_rel = NAN;

	/* verification variables for integration input */
	int 	number_integration_inputs = 3;
	bool	input_integration_received[number_integration_inputs];
	for (int i = 0; i < number_integration_inputs; i++)
	{
		input_integration_received[i] = false;
	}

	/* reading integration specs from user */
	FILE *in2 = fopen(simulation->integration_specs, "r");
   	while ((read = getline(&line, &len, in2)) != -1)
	{
		sscanf(line, "%s %s", first_col, second_col);
		if (strcmp(first_col, "t_init(yr)") == 0)
		{
			simulation->t_init = atof(second_col);
			input_integration_received[0] = true;
		}
		else if (strcmp(first_col, "t_trans(yr)") == 0)
		{
			simulation->t_trans = atof(second_col);
			input_integration_received[1] = true;
		}
		else if (strcmp(first_col, "t_final(yr)") == 0)
		{
			simulation->t_final = atof(second_col);
			input_integration_received[2] = true;
		}
		else if (strcmp(first_col, "t_step(yr)") == 0)
		{
			simulation->t_step = atof(second_col);
			simulation->t_step_received = true;
		}
	}
	fclose(in2);

	/* parameter input verification */
	for (int i = 0; i < number_integration_inputs; i++)
	{
		if(input_integration_received[i] == false)
		{
			fprintf(stderr, "Error: there is at least one missing input ");
			fprintf(stderr, "from %s.\n", simulation->integration_specs);
			exit(13);
		}
	}

	/* gravitational constant */
	simulation->G = 4.0 * M_PI * M_PI; // AU Msun year

	/* simulation units */
	strcpy(simulation->units, "AU_Msun_year");

	/* read dev file */
	FILE *in3 = fopen(simulation->dev_specs, "r");
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
					simulation->G = 6.6743e-11; // SI
					strcpy(simulation->units, second_col);
				}
				else if (strcmp(second_col, "G_unity") == 0)
				{
					simulation->G = 1.0; // non-dimensional
					strcpy(simulation->units, second_col);
				}
			}
		}
		fclose(in3);
	}

	return 0;
}

int
fill_in_bodies_data	(cltbdy	**bodies,
				 	 const siminf simulation)
{
	/* allocate memory for bodies */
	*bodies = (cltbdy *) malloc(simulation.number_of_bodies * sizeof(cltbdy));

	/* auxiliary variables for reading input */
	char 	*line = NULL;
	size_t 	len = 0;
	ssize_t read;

	/* verification variables for input */
	bool	input_mass_received = false;
	bool 	input_R_received = false;
	bool	input_orb_received = false;
	bool	input_rot_received = false;
	bool	input_initial_rot_received = false;
	/* verification variables for names and deformable settings */
	bool	input_name_received = false;
	bool	input_point_mass_received = false;
	bool	input_deformable_received = false;
	bool	input_prestress_received = false;
	bool	input_centrifugal_received = false;
	bool	input_tidal_received = false;
	/* verification variables for rheology */
	int 	number_maxwell_inputs = 3;
	bool	input_maxwell_received[number_maxwell_inputs];
	for (int i = 0; i < number_maxwell_inputs; i++)
	{
		input_maxwell_received[i] = false;
	}
	int 	number_gV_inputs = 3;
	bool	input_gV_received[number_gV_inputs];
	for (int i = 0; i < number_gV_inputs; i++)
	{
		input_gV_received[i] = false;
	}
	bool	input_elements_received = false;
	/* verification variables for orbital parameters */
	bool	input_a_received = false;
	bool	input_e_received = false;
	bool	input_I_received = false;
	bool	input_M_received = false;
	bool	input_w_received = false;
	bool	input_OMEGA_received = false;
	/* verification variables for angular velocity vector */
	bool	input_omega_azi_received = false;
	bool	input_omega_pol_received = false;
	/* verification variables for Stokes coefficients */
	bool	input_rg_received = false;
	bool	input_J2_received = false;
	bool	input_C22_received = false;
	bool	input_S22_received = false;
	bool	input_C21_received = false;
	bool	input_S21_received = false;
	/* verification variables for some orientation parameters */
	bool	input_obl_received = false;
	bool	input_psi_received = false;
	bool	input_lib_received = false;

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
			input_mass_received = true;
		}
		else if (strcmp(token, "orb(day)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].orb = atof(token);
			}
			input_orb_received = true;
		}
		else if (strcmp(token, "rot(day)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].rot = atof(token);
			}
			input_rot_received = true;
		}
		else if (strcmp(token, "rot_ini(day)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].rot_ini = atof(token);
			}
			input_initial_rot_received = true;
		}
		else if (strcmp(token, "R(km)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].R = atof(token);
			}
			input_R_received = true;
		}
		else if (strcmp(token, "rg") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].rg = atof(token);
			}
			input_rg_received = true;
		}
		else if (strcmp(token, "J2") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].J2 = atof(token);
			}
			input_J2_received = true;
		}
		else if (strcmp(token, "C22") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].C22 = atof(token);
			}
			input_C22_received = true;
		}
		else if (strcmp(token, "C21") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].C21 = atof(token);
			}
			input_C21_received = true;
		}
		else if (strcmp(token, "S21") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].S21 = atof(token);
			}
			input_S21_received = true;
		}
		else if (strcmp(token, "S22") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].S22 = atof(token);
			}
			input_S22_received = true;
		}
		else if (strcmp(token, "obl(deg)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].obl = atof(token);
			}
			input_obl_received = true;
		}
		else if (strcmp(token, "psi(deg)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].psi = atof(token);
			}
			input_psi_received = true;
		}
		else if (strcmp(token, "lib(deg)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].lib = atof(token);
			}
			input_lib_received = true;
		}
		else if (strcmp(token, "a(AU)") == 0)
		{
			(*bodies)[0].a = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].a = atof(token);
			}
			input_a_received = true;
		}
		else if (strcmp(token, "e") == 0)
		{
			(*bodies)[0].e = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].e = atof(token);
			}
			input_e_received = true;			
		}
		else if (strcmp(token, "I(deg)") == 0)
		{
			(*bodies)[0].I = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].I = atof(token);
			}
			input_I_received = true;
		}
		else if (strcmp(token, "M(deg)") == 0)
		{
			(*bodies)[0].M = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].M = atof(token);
			}
			input_M_received = true;
		}
		else if (strcmp(token, "w(deg)") == 0)
		{
			(*bodies)[0].w = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].w = atof(token);
			}
			input_w_received = true;
		}
		else if (strcmp(token, "OMEGA(deg)") == 0)
		{
			(*bodies)[0].Omega = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].Omega = atof(token);
			}
			input_OMEGA_received = true;
		}
		else if (strcmp(token, "azi(deg)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].azi = atof(token);
			}
			input_omega_azi_received = true;
		}
		else if (strcmp(token, "pol(deg)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].pol = atof(token);
			}
			input_omega_pol_received = true;
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
			input_point_mass_received = true;
		}
		else if (strcmp(token, "deformable") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				if (strcmp(token, "yes") == 0)
				{
					(*bodies)[i].deformable = true;
				}
				else if (strcmp(token, "no") == 0)
				{
					(*bodies)[i].deformable = false;
				}
				else
				{
					fprintf(stderr, "Please provide yes or no ");
					fprintf(stderr, "for deformable variable\n");
					exit(14);
				}
			}
			input_deformable_received = true;
		}
		else if (strcmp(token, "prestress") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				if (strcmp(token, "yes") == 0)
				{
					(*bodies)[i].prestress = true;
				}
				else if (strcmp(token, "no") == 0)
				{
					(*bodies)[i].prestress = false;
				}
				else
				{
					fprintf(stderr, "Please provide yes or no ");
					fprintf(stderr, "for prestress variable\n");
					exit(14);
				}
			}
			input_prestress_received = true;
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
			input_centrifugal_received = true;
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
			input_tidal_received = true;
		}
		else if (strcmp(token, "k0") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].k0 = atof(token);
			}
			input_maxwell_received[0] = true;
		}
		else if (strcmp(token, "Dt(s)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].Dt = atof(token);
			}
			input_maxwell_received[1] = true;
		}
		else if (strcmp(token, "tau(yr)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].tau = atof(token);
			}
			input_maxwell_received[2] = true;
		}
		else if (strcmp(token, "gamma_0") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].gamma_0 = atof(token);
			}
			input_gV_received[0] = true;
		}
		else if (strcmp(token, "alpha") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].alpha = atof(token);
			}
			input_gV_received[1] = true;
		}
		else if (strcmp(token, "eta") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].eta = atof(token);
			}
			input_gV_received[2] = true;
		}
		else if (strcmp(token, "elements") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].elements = atoi(token);
			}
			input_elements_received = true;
		}
	}
	fclose(in1);

	/* parameter input verification */
	if(input_mass_received == false)
	{
		fprintf(stderr, "Error: mass missing ");
		fprintf(stderr, "from %s.\n", simulation.system_specs);
		fprintf(stderr, "Exiting the program now.\n");
		exit(14);
	}
	if(input_a_received == false && simulation.number_of_bodies != 1)
	{
		fprintf(stderr, "Error: semi-major axis missing ");
		fprintf(stderr, "from %s.\n", simulation.system_specs);
		fprintf(stderr, "Exiting the program now.\n");
		exit(14);
	}
	if(input_initial_rot_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].rot_ini = (*bodies)[i].rot;
		}
	}
	if(input_e_received == false)
	{
		(*bodies)[0].e = NAN;
		for (int i = 1; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].e = 0.0;
		}
	}
	if(input_I_received == false)
	{
		(*bodies)[0].I = NAN;
		for (int i = 1; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].I = 0.0;
		}
	}
	if(input_M_received == false)
	{
		(*bodies)[0].M = NAN;
		for (int i = 1; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].M = 0.0;
		}
	}
	if(input_w_received == false)
	{
		(*bodies)[0].w = NAN;
		for (int i = 1; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].w = 0.0;
		}
	}
	if(input_OMEGA_received == false)
	{
		(*bodies)[0].Omega = NAN;
		for (int i = 1; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].Omega = 0.0;
		}
	}
	if(input_omega_azi_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].azi = 0.0;
		}
	}
	if(input_omega_pol_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].pol = 0.0;
		}
	}
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{	
		if(input_J2_received == false)
		{
			(*bodies)[i].J2 = 0.0;
		}
		else if (fabs((*bodies)[i].J2) <= 1.0e-13)
		{
			(*bodies)[i].J2 = 0.0;
		}
		if(input_C22_received == false)
		{
			(*bodies)[i].C22 = 0.0;
		}
		else if (fabs((*bodies)[i].C22) <= 1.0e-13)
		{
			(*bodies)[i].C22 = 0.0;
		}
		if(input_S22_received == false)
		{
			(*bodies)[i].S22 = 0.0;
		}
		else if (fabs((*bodies)[i].S22) <= 1.0e-13)
		{
			(*bodies)[i].S22 = 0.0;
		}
		if(input_C21_received == false)
		{
			(*bodies)[i].C21 = 0.0;
		}
		else if (fabs((*bodies)[i].C21) <= 1.0e-13)
		{
			(*bodies)[i].C21 = 0.0;
		}
		if(input_S21_received == false)
		{
			(*bodies)[i].S21 = 0.0;
		}
		else if (fabs((*bodies)[i].S21) <= 1.0e-13)
		{
			(*bodies)[i].S21 = 0.0;
		}
	}
	if(input_obl_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].obl = 0.0;
		}
	}
	if(input_psi_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].psi = 0.0;
		}
	}
	if(input_lib_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].lib = 0.0;
		}
	}
	/* elements input verification */
	if(input_elements_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].elements = 0;
		}
	}
	/* name input verification */
	if (input_name_received == false)
	{
		fprintf(stderr, "Error: could not read body names ");
		fprintf(stderr, "from %s.\n", simulation.system_specs);
		fprintf(stderr, "Exiting the program now.\n");
		exit(14);
	}
	/* deformable variables input verification */
	if(input_point_mass_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].point_mass = true;
		}
	}
	if(input_deformable_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].deformable = false;
		}
	}
	if(input_prestress_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].prestress = false;
		}
	}
	if(input_centrifugal_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].centrifugal = false;
		}
	}
	if(input_tidal_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{	
			(*bodies)[i].tidal = false;
		}
	}
	bool	all_point_masses = true;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{	
		if((*bodies)[i].point_mass == false)
		{
			all_point_masses = false;
			break;
		}
	}
	if (all_point_masses == false)
	{
		/* Essential input verification for non-point masses */
		if(input_rot_received == false)
		{
			fprintf(stderr, "Error: rotational period missing ");
			fprintf(stderr, "from %s.\n", simulation.system_specs);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
		if(input_R_received == false)
		{
			fprintf(stderr, "Error: radius missing ");
			fprintf(stderr, "from %s.\n", simulation.system_specs);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
		if(input_rg_received == false)
		{
			fprintf(stderr, "Error: rg missing ");
			fprintf(stderr, "from %s.\n", simulation.system_specs);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
	}
	bool	any_deformable = false;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{	
		if((*bodies)[i].deformable == true)
		{
			any_deformable = true;
			break;
		}
	}
	/* Maxwell rheology input verification */
	if (strcmp(simulation.rheology_model, "Maxwell") == 0)
	{
		if (any_deformable == true)	
		{	
			for (int n = 0; n < number_maxwell_inputs; n++)
			{
				if(input_maxwell_received[n] == false)
				{
					fprintf(stderr, "Error: there is at least one missing input ");
					fprintf(stderr, "for the Maxwell rheology variables ");
					fprintf(stderr, "from %s.\n", simulation.system_specs);
					fprintf(stderr, "Exiting the program now.\n");
					exit(14);
				}
			}
		}
		// zeroing Kelvin-Voigt elements
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].elements = 0;
		}
	}
	/* Generalised Voigt rheology input verification and setting of elements */
	if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
	{
		if (any_deformable == true)	
		{	
			for (int n = 0; n < number_gV_inputs; n++)
			{
				if(input_gV_received[n] == false)
				{
					fprintf(stderr, "Error: there is at least one missing input ");
					fprintf(stderr, "for the generalised Voigt rheology variables ");
					fprintf(stderr, "from %s.\n", simulation.system_specs);
					fprintf(stderr, "Exiting the program now.\n");
					exit(14);
				}
			}
		}

		// allocate memory for Voigt element
		int elements_max = 0;
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			if ((*bodies)[i].point_mass == true || (*bodies)[i].deformable == false)
			{
				(*bodies)[i].elements = 0;
			}
			if ((*bodies)[i].elements > 0)
			{
				(*bodies)[i].alpha_elements 
					= (double *) malloc((*bodies)[i].elements * sizeof(double));
				(*bodies)[i].eta_elements 	
					= (double *) malloc((*bodies)[i].elements * sizeof(double));
			}
			if ((*bodies)[i].elements > elements_max)
			{
				elements_max = (*bodies)[i].elements;
			}
		}

		/* verification for input variables */
		bool **input_element_alpha_received, **input_element_eta_received;
		if (elements_max > 0)
		{
			input_element_alpha_received = malloc(simulation.number_of_bodies * sizeof(bool*));
			input_element_eta_received	 = malloc(simulation.number_of_bodies * sizeof(bool*));
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				if ((*bodies)[i].elements > 0)
				{
					input_element_alpha_received[i] = malloc((*bodies)[i].elements * sizeof(bool));
					input_element_eta_received[i] 	= malloc((*bodies)[i].elements * sizeof(bool));
					for (int j = 0; j < (*bodies)[i].elements; j++)
					{
						input_element_alpha_received[i][j] = false;
						input_element_eta_received[i][j] = false;
					}
				}
			}

			char name_element_alpha[20], name_element_eta[20];
			FILE *in1_elements = fopen(simulation.system_specs, "r");
			while ((read = getline(&line, &len, in1_elements)) != -1)
			{
				const char tok_del[6] = " \t\n";	// token delimiter
				char *token = strtok(line, tok_del);
				if(token == NULL) break;	// in case there is a newline
				for (int j = 0; j < elements_max; j++)
				{
					sprintf(name_element_alpha, "alpha_%d", j+1);
					sprintf(name_element_eta, "eta_%d", j+1);
					if (strcmp(token, name_element_alpha) == 0)
					{
						for (int i = 0; i < simulation.number_of_bodies; i++)
						{
							token = strtok(NULL, tok_del);
							if (j < (*bodies)[i].elements)
							{
								(*bodies)[i].alpha_elements[j] = atof(token);
							}
							if ((*bodies)[i].elements > 0)
							{
								input_element_alpha_received[i][j] = true;
							}
						}
					}
					else if (strcmp(token, name_element_eta) == 0)
					{
						for (int i = 0; i < simulation.number_of_bodies; i++)
						{
							token = strtok(NULL, tok_del);
							if (j < (*bodies)[i].elements)
							{
								(*bodies)[i].eta_elements[j] = atof(token);
							}
							if ((*bodies)[i].elements > 0)
							{
								input_element_eta_received[i][j] = true;
							}
						}
					}
				}
			}
			fclose(in1_elements);

			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				if ((*bodies)[i].elements > 0)
				{
					for (int j = 0; j < (*bodies)[i].elements; j++)
					{
						if (input_element_alpha_received[i][j] == false ||
							input_element_eta_received[i][j] == false)
						{
							fprintf(stderr, "Error: there are missing parameters ");
							fprintf(stderr, "for the generalised Voigt rheology ");
							fprintf(stderr, "from %s.\n", simulation.system_specs);
							fprintf(stderr, "Exiting the program now.\n");
							exit(15);				
						}
					}
					free(input_element_alpha_received[i]);
					free(input_element_eta_received[i]);
				}
			}
			free(input_element_alpha_received);
			free(input_element_eta_received);
		}
	}

	/* hierarchy conditions (just to be safe) */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if ((*bodies)[i].point_mass == true)
		{
			(*bodies)[i].deformable = false;
			(*bodies)[i].prestress = false;
			(*bodies)[i].centrifugal = false;
			(*bodies)[i].tidal = false;
		}
		else if ((*bodies)[i].deformable == false)
		{
			(*bodies)[i].prestress = false;
			(*bodies)[i].centrifugal = false;
			(*bodies)[i].tidal = false;
		}
	}

	/* converting units and setting orbital period */
	/* also transforming some values to nan */

	double deg_to_rad = M_PI / 180.0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if ((*bodies)[i].point_mass == true)
		{
			(*bodies)[i].azi = NAN;
			(*bodies)[i].pol = NAN;
			(*bodies)[i].obl = NAN;
			(*bodies)[i].psi = NAN;
			(*bodies)[i].lib = NAN;
		}
		else
		{
			(*bodies)[i].azi *= deg_to_rad;
			(*bodies)[i].pol *= deg_to_rad;
			(*bodies)[i].obl *= deg_to_rad;
			(*bodies)[i].psi *= deg_to_rad;
			(*bodies)[i].lib *= deg_to_rad;
		}
		if (i > 0)
		{
			(*bodies)[i].I *= deg_to_rad;
			(*bodies)[i].M *= deg_to_rad;
			(*bodies)[i].w *= deg_to_rad;
			(*bodies)[i].Omega *= deg_to_rad;
		}
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
			if ((*bodies)[i].point_mass == true)
			{
				(*bodies)[i].rot = NAN;
				(*bodies)[i].rot_ini = NAN;
			}
			else
			{
				(*bodies)[i].rot *= day_to_s;
				(*bodies)[i].rot_ini *= day_to_s;
			}
			(*bodies)[i].a *= AU_to_m;
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				(*bodies)[i].tau *= year_to_s;
			}
			if (i==0)
			{
				(*bodies)[i].orb = NAN;	
			}
			else
			{
				if(input_orb_received == true)
				{
					(*bodies)[i].orb *= day_to_s;
				}
				else
				{
					(*bodies)[i].orb = kepler_period((*bodies)[0].mass, 
						(*bodies)[i].mass, simulation.G, (*bodies)[i].a);
				}
			}
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
			if ((*bodies)[i].point_mass == true)
			{
				(*bodies)[i].rot = NAN;
				(*bodies)[i].rot_ini = NAN;
			}
			else
			{
				(*bodies)[i].rot *= day_to_year;
				(*bodies)[i].rot_ini *= day_to_year;
			}
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				(*bodies)[i].Dt *= s_to_year;
			}
			if (i==0)
			{
				(*bodies)[i].orb = NAN;	
			}
			else
			{
				if(input_orb_received == true)
				{
					(*bodies)[i].orb *= day_to_year;
				}
				else
				{
					(*bodies)[i].orb = kepler_period((*bodies)[0].mass, 
						(*bodies)[i].mass, simulation.G, (*bodies)[i].a);
				}
			}
		}
	}

	/* verifying the validity of some given input */

	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if ((*bodies)[i].deformable == true)
		{
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				if (fabs((*bodies)[i].k0) < 1e-10)
				{	
					fprintf(stderr, "Warning: k0 cannot be zero.\n");
					fprintf(stderr, "Exiting the program now.\n");
					exit(16);
				}
				if (fabs((*bodies)[i].tau) < 1e-10)
				{	
					fprintf(stderr, "Warning: tau cannot be zero.\n");
					fprintf(stderr, "Exiting the program now.\n");
					exit(16);
				}
				if ((*bodies)[i].tau < (*bodies)[i].Dt)
				{	
					fprintf(stderr, "Warning: tau cannot be lower than Dt.\n");
					fprintf(stderr, "Exiting the program now.\n");
					exit(16);
				}
				if (fabs((*bodies)[i].tau-(*bodies)[i].Dt) < 1e-10)
				{	
					fprintf(stderr, "Warning: tau cannot be equal to Dt.\n");
					fprintf(stderr, "Exiting the program now.\n");
					exit(16);
				}
			}
			if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
			{
				if (fabs((*bodies)[i].gamma_0) < 1e-10)
				{	
					fprintf(stderr, "Warning: gamma cannot be zero.\n");
					fprintf(stderr, "Exiting the program now.\n");
					exit(16);
				}
				if (fabs((*bodies)[i].eta) < 1e-10)
				{	
					fprintf(stderr, "Warning: eta cannot be zero.\n");
					fprintf(stderr, "Exiting the program now.\n");
					exit(16);
				}
				for (int j = 0; j < (*bodies)[i].elements; j++)
				{
					if (fabs((*bodies)[i].eta_elements[j]) < 1e-10)
					{	
						fprintf(stderr, "Warning: all eta elements have ");
						fprintf(stderr, "to be greater than zero.\n");
						fprintf(stderr, "Exiting the program now.\n");
						exit(16);
					}	
				}
			}
		}
	}
	
	/* calculating 4 sets of variables */

	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		double	m = (*bodies)[i].mass;
		double	R = (*bodies)[i].R;
		double	T = (*bodies)[i].orb;

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

		double	k0 = (*bodies)[i].k0;
		double	Dt = (*bodies)[i].Dt;
		double	tau = (*bodies)[i].tau;
		
		/* 1st set of variables - x and x_dot */
		double qr_1_I[4];
		double qr_3_Omega[4];
		if (simulation.number_of_bodies == 1)
		{
			nan_vector((*bodies)[i].x);
			nan_vector((*bodies)[i].x_dot);
			identity_quaternion(qr_1_I);		// for 2nd set of variables
			identity_quaternion(qr_3_Omega);	// for 2nd set of variables
		}
		else
		{
			if (i == 0)
			{
				null_vector((*bodies)[i].x);
				null_vector((*bodies)[i].x_dot);
				identity_quaternion(qr_1_I);		// for 2nd set of variables
				identity_quaternion(qr_3_Omega);	// for 2nd set of variables
			}
			else
			{
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
		}

		/* 2nd set of variables - omega and q */
		if ((*bodies)[i].point_mass == true)
		{
			nan_quaternion((*bodies)[i].q);
			nan_vector((*bodies)[i].omega);
		}
		else
		{
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

			initialize_angular_velocity_on_z_axis(&(*bodies)[i]);
		}

		/* 3rd set of variables - I0 */
		if ((*bodies)[i].point_mass == true)
		{
			(*bodies)[i].I0 = NAN;
		}
		else
		{
			(*bodies)[i].I0 = (3.0 * rg - 2.0 * J2) * m * R * R / 3.0;
		}

		/* 4th set of variables - gamma_0, alpha and eta */
		if ((*bodies)[i].deformable == true)
		{
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				(*bodies)[i].gamma_0 = parameter_gamma_0(simulation.G, (*bodies)[i].I0, R, k0);
				(*bodies)[i].alpha 	 = (*bodies)[i].gamma_0 * Dt / (tau - Dt);
				(*bodies)[i].eta   	 = (*bodies)[i].gamma_0 * Dt;
			}
			else if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
			{
				(*bodies)[i].k0  = NAN;
				(*bodies)[i].tau = NAN;
				(*bodies)[i].Dt  = NAN;
			}
		}
		else
		{
			(*bodies)[i].k0  = NAN;
			(*bodies)[i].tau = NAN;
			(*bodies)[i].Dt  = NAN;

			(*bodies)[i].gamma_0 = NAN;
			(*bodies)[i].alpha 	 = NAN;
			(*bodies)[i].eta   	 = NAN;
		}
	} // end loop over bodies

	/* redefine center of inertial frame to the center of mass */
	if (simulation.number_of_bodies > 1)
	{
		double location_of_center_of_mass[3];
		double velocity_of_center_of_mass[3];
		calculate_center_of_mass(location_of_center_of_mass,
			velocity_of_center_of_mass, *bodies, 
			simulation.number_of_bodies);
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			linear_combination_vector((*bodies)[i].x,
				1.0, (*bodies)[i].x,
				-1.0, location_of_center_of_mass);
			linear_combination_vector((*bodies)[i].x_dot,
				1.0, (*bodies)[i].x_dot,
				-1.0, velocity_of_center_of_mass);
		} // end loop over bodies
	}

	/* calculate Bs_me and P_me, and initialize b_eta and bk */
	/* also update or turn into NAN the Stokes coefficients */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if ((*bodies)[i].point_mass == true)
		{
			double dummy_nan_matrix[9];
			nan_matrix(dummy_nan_matrix);			
			get_main_elements_traceless_symmetric_matrix((*bodies)[i].Bs_me,
				dummy_nan_matrix);
			(*bodies)[i].rg  = NAN;
			(*bodies)[i].J2  = NAN;
			(*bodies)[i].C22 = NAN;
			(*bodies)[i].S22 = NAN;
			(*bodies)[i].C21 = NAN;
			(*bodies)[i].S21 = NAN;
			get_main_elements_traceless_symmetric_matrix((*bodies)[i].P_me,
				dummy_nan_matrix);
			nan_matrix((*bodies)[i].Y);
			nan_matrix((*bodies)[i].Y_trans);
			get_main_elements_traceless_symmetric_matrix((*bodies)[i].p_me,
				dummy_nan_matrix);
			get_main_elements_traceless_symmetric_matrix((*bodies)[i].bs_me,
				dummy_nan_matrix);
			get_main_elements_traceless_symmetric_matrix((*bodies)[i].b_eta_me, 
				dummy_nan_matrix);
		}
		else
		{
			/* real deformation on Tisserand frame */
			double B_stokes_i[9];
			double B_stokes_diag_i[9];
			// initial spherical configuration
			null_matrix(B_stokes_i);
			null_matrix(B_stokes_diag_i);
			// initial non-elipsoidal configuration
			if ((fabs((*bodies)[i].J2)  > 1.0e-13) ||
				(fabs((*bodies)[i].C22) > 1.0e-13) ||
				(fabs((*bodies)[i].S22) > 1.0e-13) ||
				(fabs((*bodies)[i].C21) > 1.0e-13) ||
				(fabs((*bodies)[i].S21) > 1.0e-13))
			{
				body_frame_deformation_from_stokes_coefficients(B_stokes_i, (*bodies)[i]);

				if ((fabs((*bodies)[i].S22) > 1.0e-13) ||
					(fabs((*bodies)[i].C21) > 1.0e-13) ||
					(fabs((*bodies)[i].S21) > 1.0e-13))
				{
					calculate_diagonalized_square_matrix(B_stokes_diag_i, B_stokes_i);
					
					/* update gravity field coefficients to */
					/* the frame of principal inertia momenta */
					double Iner_diag[9];
					calculate_inertia_tensor(Iner_diag, (*bodies)[i].I0, B_stokes_diag_i);
					(*bodies)[i].rg  = calculate_rg ((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
					(*bodies)[i].J2  = calculate_J2 ((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
					(*bodies)[i].C22 = calculate_C22((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
					(*bodies)[i].S22 = calculate_S22((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
					(*bodies)[i].C21 = calculate_C21((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
					(*bodies)[i].S21 = calculate_S21((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
				}
				else
				{
					copy_square_matrix(B_stokes_diag_i, B_stokes_i);
				}
			}
			get_main_elements_traceless_symmetric_matrix((*bodies)[i].Bs_me,
				B_stokes_diag_i);

			/* prestress on Tisserand frame */
			if ((*bodies)[i].prestress == false)
			{
				double dummy_nan_matrix[9];
				nan_matrix(dummy_nan_matrix);
				get_main_elements_traceless_symmetric_matrix((*bodies)[i].P_me,
					dummy_nan_matrix);
			}
			else
			{
				double F_cent_mean_i[9];
				double mean_omega = 2.0 * M_PI / (*bodies)[i].rot;
				calculate_F_cent_mean(F_cent_mean_i, mean_omega);
				double F_tide_mean_i[9];
				null_matrix(F_tide_mean_i); // no permanent tide for now
				double P_i[9];	
				linear_combination_three_square_matrix(P_i,
					(*bodies)[i].gamma_0, B_stokes_diag_i,
					-1.0, F_cent_mean_i,
					-1.0, F_tide_mean_i);
				get_main_elements_traceless_symmetric_matrix((*bodies)[i].P_me,
					P_i);
			}

			/* real deformation and prestress on inertial frame */
			calculate_Y_and_Y_transpose(&(*bodies)[i]);
			calculate_bs_me(&(*bodies)[i]);
			if ((*bodies)[i].prestress == false)
			{
				double dummy_nan_matrix[9];
				nan_matrix(dummy_nan_matrix);
				get_main_elements_traceless_symmetric_matrix((*bodies)[i].p_me,
					dummy_nan_matrix);
			}
			else
			{
				calculate_p_me(&(*bodies)[i]);
			}
			
			/* alloc and initialize bk_me with zeroes */
			if ((*bodies)[i].elements > 0)
			{
				(*bodies)[i].bk_me 
					= (double *) calloc((*bodies)[i].elements * 5, sizeof(double));
			}

			/* initialize b_eta_me */
			if ((*bodies)[i].deformable == false)
			{
				double dummy_nan_matrix[9];
				nan_matrix(dummy_nan_matrix);			
				get_main_elements_traceless_symmetric_matrix((*bodies)[i].b_eta_me, 
					dummy_nan_matrix);
			}
			else
			{
				double f_cent_i[9];
				null_matrix(f_cent_i);
				if ((*bodies)[i].centrifugal == true)
				{
					calculate_f_cent(f_cent_i, (*bodies)[i].omega);
				}
				double f_tide_i[9];
				null_matrix(f_tide_i);
				if ((*bodies)[i].tidal == true)
				{
					calculate_f_tide(f_tide_i, i, (*bodies), 
						simulation.number_of_bodies, simulation.G);
				}
				double f_ps_i[9];
				null_matrix(f_ps_i);
				if ((*bodies)[i].prestress == true)
				{
					calculate_f_ps(f_ps_i, (*bodies)[i]);
				}
				double c_i = calculate_c((*bodies)[i]);
				double b_eta_i[9];
				null_matrix(b_eta_i);
				if ((*bodies)[i].deformable == true)
				{
					double bs_i[9];
					construct_traceless_symmetric_matrix(bs_i,
						(*bodies)[i].bs_me);
					linear_combination_four_square_matrix(b_eta_i,
						c_i/(*bodies)[i].alpha, bs_i,
						-1.0/(*bodies)[i].alpha, f_cent_i,
						-1.0/(*bodies)[i].alpha, f_tide_i,
						-1.0/(*bodies)[i].alpha, f_ps_i);
				}
				get_main_elements_traceless_symmetric_matrix((*bodies)[i].b_eta_me, b_eta_i);
			}
		}
	} // end loop over bodies

	/* calculate b, and initialize l and angular velocity */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if ((*bodies)[i].point_mass == true)
		{
			nan_matrix((*bodies)[i].b);
			nan_vector((*bodies)[i].l);
		}
		else
		{
			initialize_angular_velocity_on_z_axis(&(*bodies)[i]);
			calculate_b(i, (*bodies), simulation.number_of_bodies, simulation.G);
			initialize_angular_velocity(&(*bodies)[i]); // correct alignment
			calculate_l(&(*bodies)[i]);
		}
	} // end loop over bodies

	/* for debugging */
	// for (int i = 0; i < simulation.number_of_bodies; i++)
	// {
	// 	print_CelestialBody((*bodies)[i]);
	// }
	// exit(99);

	return 0;
}

int
calculate_data_skip (siminf *simulation,
					 const cltbdy *bodies)
{
	// checking type of bodies
	bool extended = false;
	bool deformable = false;
	for (int i = 0; i < simulation->number_of_bodies; i++)
	{
		if (bodies[i].point_mass == false)
		{
			extended = true;
			if (bodies[i].deformable == true)
			{
				deformable = true;
				break;
			}
		}
	}

	// calculating data size
	double header_size = 0.0;
	double data_size_for_each_line = 0.0;

	// everything is in bytes
	if (deformable == true)
	{
		header_size = 156.0;
		data_size_for_each_line = 22.0 * 27.0 + 1.0;
	}
	else if (extended == true)
	{
		header_size = 154.0;
		data_size_for_each_line = 22.0 * 26.0 + 1.0;
	}
	else
	{
		header_size = 75.0;
		data_size_for_each_line = 22.0 * 10.0 + 1.0;
	}

	// determine data skip
	double number_of_data_lines_in_file 
		= (simulation->output_size - header_size) / data_size_for_each_line;
	double integration_time = simulation->t_final - simulation->t_trans;
	double number_of_steps = integration_time / simulation->t_step;

	if (number_of_data_lines_in_file > number_of_steps)
	{
		simulation->data_skip = 1;
	}
	else
	{
		double data_skip_double 
			= number_of_steps / number_of_data_lines_in_file;
		if (data_skip_double < (double) INT_MAX)
		{
			simulation->data_skip = (int) data_skip_double;
		}
		else
		{
			fprintf(stderr, "Error: could not calculate a data skip");
			fprintf(stderr, " to guarantee a max data size of %f mb.\n", 
				simulation->output_size / 1e6);
			fprintf(stderr, "Possible reason: final simulation time is");
			fprintf(stderr, " probably much higher than time step.");
			exit(13);
		}
	}

	return 0;
}

int
create_output_files	(const cltbdy *bodies,
			 		 const siminf simulation,
					 FILE *out[])
{
	char filename[300];
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		/* names and opens output files for each body */
		strcpy(filename, simulation.output_folder);
		strcat(filename, "results_");
		strcat(filename, simulation.name);
		strcat(filename, "_");
		strcat(filename, bodies[i].name);
		strcat(filename, ".dat");
		out[i] = fopen(filename, "w");
		/* writes output headers */
		fprintf(out[i], "time(yr)");
		fprintf(out[i], " x(AU) y(AU) z(AU) vx(Au/yr) vy(Au/yr) vz(Au/yr)");
		if (bodies[i].point_mass == false)
		{
			fprintf(out[i], " omega_x omega_y omega_z");
			fprintf(out[i], " l_x l_y l_z");
			fprintf(out[i], " b_11 b_12 b_13");
			fprintf(out[i], " b_21 b_22 b_23");
			fprintf(out[i], " b_31 b_32 b_33");
			fprintf(out[i], " q_1 q_2 q_3 q_4");
			if (bodies[i].deformable == true)
			{
				fprintf(out[i], " D");
			}
		}
		fprintf(out[i], "\n");
	}
	/* names and opens general output file */
	strcpy(filename, simulation.output_folder);
	strcat(filename, "results_");
	strcat(filename, simulation.name);
	strcat(filename, "_");
	strcat(filename, "full_system");
	strcat(filename, ".dat");
	out[simulation.number_of_bodies] = fopen(filename, "w");
	fprintf(out[simulation.number_of_bodies], "time(yr)");
	fprintf(out[simulation.number_of_bodies], " |L|");
	fprintf(out[simulation.number_of_bodies], " E");
	fprintf(out[simulation.number_of_bodies], " E_dot");
	fprintf(out[simulation.number_of_bodies], "\n");

	return 0;
}

int
write_output(const cltbdy *bodies,
			 const siminf simulation,
			 FILE *out[])
{
	double E_dot = 0.0;

	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		/* state variables */

		/* time */
		fprintf (out[i], "%.14e", simulation.t);

		/* position and velocity */
		fprintf (out[i], " %.14e %.14e %.14e %.14e %.14e %.14e", 
			bodies[i].x[0], bodies[i].x[1], bodies[i].x[2],
			bodies[i].x_dot[0], bodies[i].x_dot[1], bodies[i].x_dot[2]);

		if (bodies[i].point_mass == false)
		{
			/* angular velocity vector */
			fprintf (out[i], " %.14e %.14e %.14e", 
				bodies[i].omega[0], bodies[i].omega[1], bodies[i].omega[2]);

			/* angular momentum vector */
			fprintf (out[i], " %.14e %.14e %.14e", 
				bodies[i].l[0], bodies[i].l[1], bodies[i].l[2]);

			/* deformation matrix */
			fprintf (out[i], " %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e",
				bodies[i].b[0], bodies[i].b[1], bodies[i].b[2],
				bodies[i].b[3], bodies[i].b[4], bodies[i].b[5],
				bodies[i].b[6], bodies[i].b[7], bodies[i].b[8]);

			/* quaternion */
			fprintf (out[i], " %.14e %.14e %.14e %.14e", 
				bodies[i].q[0], bodies[i].q[1], bodies[i].q[2], bodies[i].q[3]);

			/* dissipation on body */
			if (bodies[i].deformable == true)
			{
				double D_i = dissipation_function(bodies[i]);
				fprintf (out[i], " %.14e", D_i);
				E_dot -= 2.0 * D_i;
			}
		}

		/* next line */
		fprintf(out[i], "\n");

	} // end loop over bodies

	/* prints on general file */

	/* time */
	fprintf (out[simulation.number_of_bodies], "%.14e", simulation.t);

	/* total angular momentum */
	double l_total[3];
	calculate_total_angular_momentum(l_total, bodies, simulation.number_of_bodies);
	fprintf (out[simulation.number_of_bodies], " %.14e", norm_vector(l_total));

	/* energy without deformation */
	fprintf (out[simulation.number_of_bodies], " %.14e", 
		total_energy_without_deformation(bodies, 
			simulation.number_of_bodies, simulation.G));

	/* total energy dissipation */
	fprintf (out[simulation.number_of_bodies], " %.14e", E_dot);

	/* next line */
	fprintf(out[simulation.number_of_bodies], "\n");

	return 0;
}

int
close_output_files	(const siminf simulation,
					 FILE *out[])
{
	for (int i = 0; i < simulation.number_of_bodies + 1; i++)
	{
		fclose(out[i]);
	}

	return 0;
}

int
write_simulation_overview	(const int time_spent_in_seconds,
							 const siminf simulation)
{
	// creates info file
	char filename[300];
	strcpy(filename, simulation.output_folder);
	strcat(filename, "results_");
	strcat(filename, simulation.name);
	strcat(filename, "_");
	strcat(filename, "sim_info");
	strcat(filename, ".dat");

	// simulation time
	int sec, min, hr, day;
	day = time_spent_in_seconds / (24*3600);
	hr	= (time_spent_in_seconds - 24*3600*day) / 3600;
	min = (time_spent_in_seconds - 24*3600*day - 3600*hr) / 60;
	sec = (time_spent_in_seconds - 24*3600*day - 3600*hr - 60*min) / 1;
	FILE *out_sim_info;
	out_sim_info = fopen(filename, "w");
	fprintf(out_sim_info, "Time spent on simulation:");
	fprintf(out_sim_info, " %d days %d hours %d minutes %d seconds.\n\n",
		day, hr, min, sec);
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
	FILE *in2_to_copy = fopen(simulation.integration_specs, "r");
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

	// parameters calculated by the code
	FILE *out_sim_info_2;
	out_sim_info_2 = fopen(filename, "a");
	fprintf(out_sim_info_2, "t_step(yr) = %1.10e\n", simulation.t_step);
	fprintf(out_sim_info_2, "data_skip = %d\n\n", simulation.data_skip);
	fclose(out_sim_info_2);

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
	char filename_read_first_body[300];
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
		char filename_read[300];
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
		char filename_orbital[300];
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
		fprintf(out_orbital, " x_orb(AU) y_orb(AU) z_orb(AU)");
		fprintf(out_orbital, " |omega_orbit|");
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
			fprintf (out_orbital, "%.14e", t);

			/* semimajor axis */
			fprintf (out_orbital, " %.14e", bodies[i].a);

			/* orbital eccentricity */
			fprintf (out_orbital, " %.14e", bodies[i].e);

			/* inclination */
			fprintf (out_orbital, " %.14e", bodies[i].I * rad_to_deg);

			/* mean anomaly */
			fprintf (out_orbital, " %.14e", bodies[i].M * rad_to_deg);

			/* argument of periapsis */
			fprintf (out_orbital, " %.14e", bodies[i].w * rad_to_deg);

			/* longitude of the ascending node */
			fprintf (out_orbital, " %.14e", bodies[i].Omega * rad_to_deg);

			if (i > 0)
			{
				/* position in 0-centered system */
				double relative_x[3];
				linear_combination_vector(relative_x,
					1.0, bodies[i].x,
					-1.0, bodies[0].x);
				fprintf (out_orbital, " %.14e %.14e %.14e", 
					relative_x[0], relative_x[1], relative_x[2]);

				/* orbital frequency */
				double relative_x_dot[3];
				linear_combination_vector(relative_x_dot,
					1.0, bodies[i].x_dot,
					-1.0, bodies[0].x_dot);	
				double relative_x_cross_relative_x_dot[3];
				cross_product(relative_x_cross_relative_x_dot, 
					relative_x, relative_x_dot);
				double orbital_frequency[3];
				scale_vector(orbital_frequency, 
					1.0 / norm_squared_vector(relative_x),
					relative_x_cross_relative_x_dot);
				fprintf (out_orbital, " %.14e", norm_vector(orbital_frequency));
			}

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
	char filename_read_first_body[300];
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
			char filename_read[300];
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
			char filename_orientation[300];
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
			fprintf(out_orientation, " wI3(°)"); // angle w and I3
			fprintf(out_orientation, " wI3sf(°)"); // angle w and I3 solid frame
			fprintf(out_orientation, " wl(°)"); // angle w and l
			fprintf(out_orientation, " azi(°)"); // nutation angle
			fprintf(out_orientation, " aziPIM(°)"); // nutation angle on Principal
													 // Inertia Momenta (PIM) frame
			fprintf(out_orientation, " J2"); // J2 on PIM frame
			fprintf(out_orientation, " C22"); // C22 on PIM frame
			fprintf(out_orientation, " |q|");
			if (i > 0)
			{
				fprintf(out_orientation, " soa(°)"); // angle relative coordinate 
												 	 // and I1 (spin-orbit angle)
			}
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
				for (int j = 0; j < 9; j++)
				{
					token = strtok_r(line_track, tok_del, &line_track);
					bodies[i].b[j] = atof(token);
				}
				for (int j = 0; j < 4; j++)
				{
					token = strtok_r(line_track, tok_del, &line_track);
					bodies[i].q[j] = atof(token);
				}

				/* time */
				fprintf (out_orientation, "%.14e", t);

				/* norm of angular velocity vector */
				fprintf (out_orientation, " %.14e", norm_vector(bodies[i].omega));

				/* norm of angular momentum vector */
				fprintf (out_orientation, " %.14e", norm_vector(bodies[i].l));

				/* norm of deformation matrix */
				fprintf (out_orientation, " %.14e", norm_square_matrix(bodies[i].b));

				/* obliquity */
				if (i == 0)
				{
					calculate_obliquity_free_body_from_angular_velocity(&bodies[i]);
					// calculate_obliquity_free_body_from_figure_axis_of_solid_frame(&bodies[i]);
				}
				else
				{
					calculate_obliquity_on_orbit_from_angular_velocity(&bodies[i], bodies[0]);
					// calculate_obliquity_on_orbit_from_figure_axis_of_solid_frame(&bodies[i], bodies[0]);
				}	
				fprintf (out_orientation, " %.14e", bodies[i].obl * rad_to_deg);

				/* angle between spin axis and figure axis */
				double angle_wI3;
				angle_wI3 = angle_between_spin_axis_and_figure_axis(bodies[i]);
				fprintf (out_orientation, " %.14e", angle_wI3 * rad_to_deg);

				/* angle between spin axis and figure axis of solid frame */
				double angle_wI3sf;
				angle_wI3sf = angle_between_spin_axis_and_figure_axis_of_solid_frame(bodies[i]);
				fprintf (out_orientation, " %.14e", angle_wI3sf * rad_to_deg);

				/* angle between spin axis and angular momentum */
				double angle_wl;
				angle_wl = angle_between_spin_axis_and_angular_momentum(bodies[i]);
				fprintf (out_orientation, " %.14e", angle_wl * rad_to_deg);

				/* nutation frequency */
				double Y_i[9], Y_i_trans[9];
				rotation_matrix_from_quaternion(Y_i, bodies[i].q);
				transpose_square_matrix(Y_i_trans, Y_i);		
				double ang_vel_body_frame[3];
				square_matrix_times_vector(ang_vel_body_frame,
					Y_i_trans, bodies[i].omega);
				double ang_vel_body_frame_spherical[3];
				cartesian_to_spherical_coordinates(ang_vel_body_frame_spherical,
					ang_vel_body_frame);
				bodies[i].azi = ang_vel_body_frame_spherical[1];
				fprintf (out_orientation, " %.14e", bodies[i].azi * rad_to_deg);

				/* nutation frequency on Principal Inertia Momenta (PIM) frame */
				double P_i[9], P_i_trans[9];
				calculate_eigenvectors_matrix(P_i, bodies[i].b);
				transpose_square_matrix(P_i_trans, P_i);	
				square_matrix_times_vector(ang_vel_body_frame,
					P_i_trans, bodies[i].omega);
				cartesian_to_spherical_coordinates(ang_vel_body_frame_spherical,
					ang_vel_body_frame);
				fprintf (out_orientation, " %.14e", 
					ang_vel_body_frame_spherical[1] * rad_to_deg);

				/* Stokes coefficients on Principal Inertia Momenta (PIM) frame */
				double B[9];
				square_matrix_times_square_matrix(B, bodies[i].b, Y_i);
				square_matrix_times_square_matrix(B, Y_i_trans, B);
				double B_diag_i[9];
				calculate_diagonalized_square_matrix(B_diag_i, B);
				double Iner_diag_i[9];
				calculate_inertia_tensor(Iner_diag_i, bodies[i].I0, B_diag_i);
				bodies[i].J2  = calculate_J2 (bodies[i].mass, bodies[i].R, Iner_diag_i);
				bodies[i].C22 = calculate_C22(bodies[i].mass, bodies[i].R, Iner_diag_i);
				fprintf (out_orientation, " %.14e %.14e", bodies[i].J2, bodies[i].C22);

				/* norm of quaternion */
				fprintf (out_orientation, " %.14e", norm_quaternion(bodies[i].q));

				/* angle between relative coordinate and I1 */
				if (i > 0)
				{
					double soa = angle_between_relative_x_and_I1(bodies[i], bodies[0]);
					fprintf (out_orientation, " %.14e", soa * rad_to_deg);
				}

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
	/*
	 * suppresses warning messages by moving
	 * stderr to /dev/null
	 * should be used with care
	**/
	if(freopen("/dev/null", "w", stderr) == NULL)
	{
		fprintf(stderr, "Error: stderr could not be moved to\n");
		fprintf(stderr, "/dev/null in plot command.\n");
		exit(22);
	}
	/* create figures dir */
	char plot_folder[300];
	strcpy(plot_folder, simulation.output_folder);
	strcat(plot_folder, "/figures/");
	struct stat st_plot = {0};
	if (stat(plot_folder, &st_plot) == -1) {
		mkdir(plot_folder, 0700);
	}
	char filename_get[300];
	char filename_plot[300];
	FILE *gnuplotPipe;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].point_mass == false && bodies[i].deformable == true)
		{
			/* get output files of each body */
			strcpy(filename_get, simulation.output_folder);
			strcat(filename_get, "results_");
			strcat(filename_get, simulation.name);
			strcat(filename_get, "_");
			strcat(filename_get, bodies[i].name);
			strcat(filename_get, ".dat");
			/* plot outputs */
			/* energy dissipation */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_energy_dissipation");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"D\"\n");
			fprintf(gnuplotPipe, "set title \"%s's energy dissipation\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:27 w l lw 3", filename_get);
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"{/Symbol W}(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's longitude of the ascending node\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:7 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* position in 0-centered system */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_orbit");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"x_{orb}(AU)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"y_{orb}(AU)\"\n");
			fprintf(gnuplotPipe, "set zlabel \"z_{orb}(AU)\" rotate parallel\n");
			fprintf(gnuplotPipe, "set title \"%s's orbit\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "splot \'%s\' u 8:9:10 w d", filename_get);
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"Obliquity(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's obliquity\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:5 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* angle between spin axis and figure axis */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_wI3");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"wI3(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's wI3\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:6 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* angle between spin axis and figure axis of solid frame */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_wI3sf");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"wI3sf(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's wI3sf\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:7 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* angle between spin axis and angular momentum */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_wl");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"wl(°)\"\n");
			fprintf(gnuplotPipe, "set title \"%s's wl\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:8 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* azimuthal angle */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_azimuthal_angle");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"Azimuthal angle(°)\"\n");
			fprintf(gnuplotPipe, "set title \"Azimuthal angle of %s's angular velocity\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:9 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* nutation frequency on Principal Inertia Momenta (PIM) frame */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_nutation_PIM");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"Nutation angle on PIM(°)\"\n");
			fprintf(gnuplotPipe, "set title \"Nutation angle on PIM of %s's angular velocity\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:10 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* Stokes coefficients */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_J2");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"J2\"\n");
			fprintf(gnuplotPipe, "set title \"%s's J2\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:11 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_C22");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"C22\"\n");
			fprintf(gnuplotPipe, "set title \"%s's C22\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:12 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			/* quaternion norm */
			strcpy(filename_plot, plot_folder);
			strcat(filename_plot, "figure_");
			strcat(filename_plot, simulation.name);
			strcat(filename_plot, "_");
			strcat(filename_plot, bodies[i].name);
			strcat(filename_plot, "_quaternion_norm");
			strcat(filename_plot, ".png");
			gnuplotPipe = popen("gnuplot -persistent", "w");
			fprintf(gnuplotPipe, "reset\n");
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"|q|\"\n");
			fprintf(gnuplotPipe, "set title \"%s's quaternion norm\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:13 w l lw 3", filename_get);
			pclose(gnuplotPipe);
			if (i > 0)
			{
				/* spin-orbit angle */
				strcpy(filename_plot, plot_folder);
				strcat(filename_plot, "figure_");
				strcat(filename_plot, simulation.name);
				strcat(filename_plot, "_");
				strcat(filename_plot, bodies[i].name);
				strcat(filename_plot, "_spin_orbit_angle");
				strcat(filename_plot, ".png");
				gnuplotPipe = popen("gnuplot -persistent", "w");
				fprintf(gnuplotPipe, "reset\n");
				fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
				fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
				fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
				fprintf(gnuplotPipe, "set border lw 2 \n");
				fprintf(gnuplotPipe, "unset key\n");
				fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
				fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
				fprintf(gnuplotPipe, "set ylabel \"Spin-orbit angle(°)\"\n");
				fprintf(gnuplotPipe, "set title \"%s's spin-orbit angle\"\n", bodies[i].name);
				fprintf(gnuplotPipe, "plot \'%s\' u 1:14 w l lw 3", filename_get);
				pclose(gnuplotPipe);
			} // if (i > 0)
		} // if (bodies[i].point_mass == false)
	} // end loop over bodies
	/* reads output file of full system */
	strcpy(filename_get, simulation.output_folder);
	strcat(filename_get, "results_");
	strcat(filename_get, simulation.name);
	strcat(filename_get, "_");
	strcat(filename_get, "full_system");
	strcat(filename_get, ".dat");
	/* plot outputs */
	if (simulation.number_of_bodies > 1)
	{
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
		fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
		fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
		fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
		fprintf(gnuplotPipe, "set border lw 2 \n");
		fprintf(gnuplotPipe, "unset key\n");
		fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
		fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
		fprintf(gnuplotPipe, "set ylabel \"|Total angular momentum|\"\n");
		fprintf(gnuplotPipe, "set title \"System's total angular momentum\"\n");
		fprintf(gnuplotPipe, "plot \'%s\' u 1:2 w l lw 3", filename_get);
		pclose(gnuplotPipe);
		/* total angular momentum error */
		strcpy(filename_plot, plot_folder);
		strcat(filename_plot, "figure_");
		strcat(filename_plot, simulation.name);
		strcat(filename_plot, "_");
		strcat(filename_plot, "full_system");
		strcat(filename_plot, "_total_angular_momentum_error");
		strcat(filename_plot, ".png");
		gnuplotPipe = popen("gnuplot -persistent", "w");
		fprintf(gnuplotPipe, "reset\n");
		fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
		fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
		fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
		fprintf(gnuplotPipe, "set border lw 2 \n");
		fprintf(gnuplotPipe, "unset key\n");
		fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
		fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
		fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
		fprintf(gnuplotPipe, "set ylabel \"Total angular momentum error\"\n");
		fprintf(gnuplotPipe, "set title \"System's total angular momentum error\"\n");
		fprintf(gnuplotPipe, "col = 2\n");
		fprintf(gnuplotPipe, "row = 1\n");
		fprintf(gnuplotPipe, "stats \'%s\' every ::row::row using col nooutput\n", filename_get);
		fprintf(gnuplotPipe, "value=STATS_min\n");
		fprintf(gnuplotPipe, "plot \'%s\' u 1:(abs($2-value)) w l lw 3", filename_get);
		pclose(gnuplotPipe);
	}
	/* total energy without deformation */
	strcpy(filename_plot, plot_folder);
	strcat(filename_plot, "figure_");
	strcat(filename_plot, simulation.name);
	strcat(filename_plot, "_");
	strcat(filename_plot, "full_system");
	strcat(filename_plot, "_total_energy_without_deformation");
	strcat(filename_plot, ".png");
	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
	fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
	fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
	fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Total energy without deformation\"\n");
	fprintf(gnuplotPipe, "set title \"System's total energy without deformation\"\n");
	fprintf(gnuplotPipe, "plot \'%s\' u 1:3 w l lw 3", filename_get);
	pclose(gnuplotPipe);
	/* total energy dissipation */
	strcpy(filename_plot, plot_folder);
	strcat(filename_plot, "figure_");
	strcat(filename_plot, simulation.name);
	strcat(filename_plot, "_");
	strcat(filename_plot, "full_system");
	strcat(filename_plot, "_total_energy_dissipation");
	strcat(filename_plot, ".png");
	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,25\"\n");
	fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
	fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "set format x %s\n", "\"%.1tx10^{%T}\"");
	fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Total energy dissipation\"\n");
	fprintf(gnuplotPipe, "set title \"System's total energy dissipation\"\n");
	fprintf(gnuplotPipe, "plot \'%s\' u 1:4 w l lw 3", filename_get);
	pclose(gnuplotPipe);

	return 0;
}

double
find_shortest_time_scale(const cltbdy *bodies,
						 const siminf simulation)
{
	// giant number for comparison
	double giant_number = 1.0e20;

	// define return variable
    double shortest_time_scale = giant_number;

	// orbit
	bool	orbital_motion = false;
    double 	shortest_orbital_period = giant_number;
	if (simulation.number_of_bodies > 1)
	{
		shortest_orbital_period = bodies[1].orb;
		orbital_motion = true;
		shortest_time_scale = shortest_orbital_period;
	}

	// spin
	bool	rotational_motion = false;
	double	shortest_rotational_period = giant_number;
    for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].point_mass == false)
		{
			shortest_rotational_period = bodies[i].rot;
			rotational_motion = true;
			shortest_time_scale = shortest_rotational_period;
			break;
		}
	}
	
	// rheology
	bool	deformation = false;
	double	shortest_relaxation_time = giant_number;
    for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].deformable == true)
		{
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				shortest_relaxation_time = bodies[i].tau;
			}
			else if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
			{   
				shortest_relaxation_time = bodies[i].eta / bodies[i].alpha; 
			}
			deformation = true;
			shortest_time_scale = shortest_relaxation_time;
			break;
		}
	}

	// loop over bodies
    for (int i = 0; i < simulation.number_of_bodies; i++)
    {
        // orbit
		if (orbital_motion == true)
		{
			if (i > 0)
			{
				double body_orbital_period = bodies[i].orb;
				if (body_orbital_period < shortest_orbital_period)
				{
					shortest_orbital_period = body_orbital_period;
				}
			}
		}

        // spin
		if (rotational_motion == true)
		{
			if (bodies[i].point_mass == false)
			{
				double body_rotational_period = bodies[i].rot;
				if (bodies[i].rot_ini < bodies[i].rot)
				{
					body_rotational_period = bodies[i].rot_ini;
				}
				if (body_rotational_period < shortest_rotational_period)
				{
					shortest_rotational_period = body_rotational_period;
				}
			}
		}

        // rheology
		if (deformation == true)
		{
			double body_rheology_min_time = 0.0;
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				body_rheology_min_time = bodies[i].tau;
			}
			else if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
			{   
				body_rheology_min_time = bodies[i].eta / bodies[i].alpha;
				for (int j = 0; j < bodies[i].elements; j++)
				{
					double body_rheology_min_time_element 
						= bodies[i].eta_elements[j] / bodies[i].alpha_elements[j];
					if (body_rheology_min_time_element < body_rheology_min_time)
					{
						body_rheology_min_time = body_rheology_min_time_element;
					}
				}           
			}
			if (body_rheology_min_time < shortest_relaxation_time)
			{
				shortest_relaxation_time = body_rheology_min_time;
			}
		}
    } // end loop over bodies

    if (shortest_orbital_period < shortest_time_scale)
    {
        shortest_time_scale = shortest_orbital_period;
    }
    if (shortest_rotational_period < shortest_time_scale)
    {
        shortest_time_scale = shortest_rotational_period;
    }
    if (shortest_relaxation_time < shortest_time_scale)
    {
        shortest_time_scale = shortest_relaxation_time;
    }

    return shortest_time_scale;
}

double
find_largest_time_scale(const cltbdy *bodies,
						const siminf simulation)
{
	// define return variable
    double largest_time_scale = 0.0;

	// orbit
	bool	orbital_motion = false;
    double 	largest_orbital_period = 0.0;
	if (simulation.number_of_bodies > 1)
	{
		largest_orbital_period = bodies[1].orb;
		orbital_motion = true;
		largest_time_scale = largest_orbital_period;
	}

	// spin
	bool	rotational_motion = false;
	double	largest_rotational_period = 0.0;
    for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].point_mass == false)
		{
			largest_rotational_period = bodies[i].rot;
			rotational_motion = true;
			largest_time_scale = largest_rotational_period;
			break;
		}
	}
	
	// rheology
	bool	deformation = false;
	double	largest_relaxation_time = 0.0;
    for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].deformable == true)
		{
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				largest_relaxation_time = bodies[i].tau;
			}
			else if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
			{   
				largest_relaxation_time = bodies[i].eta / bodies[i].alpha; 
			}
			deformation = true;
			largest_time_scale = largest_relaxation_time;
			break;
		}
	}

	// loop over bodies
    for (int i = 0; i < simulation.number_of_bodies; i++)
    {
        // orbit
		if (orbital_motion == true)
		{
			if (i > 0)
			{
				double body_orbital_period = bodies[i].orb;
				if (body_orbital_period > largest_orbital_period)
				{
					largest_orbital_period = body_orbital_period;
				}
			}
		}

        // spin
		if (rotational_motion == true)
		{
			if (bodies[i].point_mass == false)
			{
				double body_rotational_period = bodies[i].rot;
				if (bodies[i].rot_ini > bodies[i].rot)
				{
					body_rotational_period = bodies[i].rot_ini;
				}
				if (body_rotational_period > largest_rotational_period)
				{
					largest_rotational_period = body_rotational_period;
				}
			}
		}

        // rheology
		if (deformation == true)
		{
			double body_rheology_min_time = 0.0;
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				body_rheology_min_time = bodies[i].tau;
			}
			else if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
			{   
				body_rheology_min_time = bodies[i].eta / bodies[i].alpha;
				for (int j = 0; j < bodies[i].elements; j++)
				{
					double body_rheology_min_time_element 
						= bodies[i].eta_elements[j] / bodies[i].alpha_elements[j];
					if (body_rheology_min_time_element > body_rheology_min_time)
					{
						body_rheology_min_time = body_rheology_min_time_element;
					}
				}           
			}
			if (body_rheology_min_time > largest_relaxation_time)
			{
				largest_relaxation_time = body_rheology_min_time;
			}
		}
    } // end loop over bodies

    if (largest_orbital_period > largest_time_scale)
    {
        largest_time_scale = largest_orbital_period;
    }
    if (largest_rotational_period > largest_time_scale)
    {
        largest_time_scale = largest_rotational_period;
    }
    if (largest_relaxation_time > largest_time_scale)
    {
        largest_time_scale = largest_relaxation_time;
    }

    return largest_time_scale;
}
