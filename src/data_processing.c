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
	bool	rheology_model_received = false;
	bool	dev_specs_file_received = false;
	bool	number_of_bodies_received = false;
	bool	omega_correction_received = false;
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
		else if (strcmp(first_col, "omega_correction") == 0)
		{
			if (strcmp(second_col, "yes") == 0)
			{
				simulation->omega_correction = true;
			}
			else if (strcmp(second_col, "no") == 0)
			{
				simulation->omega_correction = false;
			}
			else
			{
				printf("%s\n", second_col);
				fprintf(stderr, "Please provide yes or no ");
				fprintf(stderr, "for omega correction\n");
				exit(14);
			}
			omega_correction_received = true;
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
		fprintf(stderr, "Error: did not receive rheology model ");
		fprintf(stderr, "from %s.\n", simulation->main_input);
		exit(13);
	}

	/* check and set omega correction */
	if (omega_correction_received == false)
	{
		simulation->omega_correction = false;
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

	/* verification variables for integration input */
	int 	number_integration_inputs = 7;
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
			input_integration_received[3] = true;
		}
		else if (strcmp(first_col, "eps_abs") == 0)
		{
			simulation->eps_abs = atof(second_col);
			input_integration_received[4] = true;
		}
		else if (strcmp(first_col, "eps_rel") == 0)
		{
			simulation->eps_rel = atof(second_col);
			input_integration_received[5] = true;
		}
		else if (strcmp(first_col, "data_skip") == 0)
		{
			simulation->data_skip = (int) atof(second_col);
			input_integration_received[6] = true;
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

	/* always initialize write to file as true */
	simulation->write_to_file = true;

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

	/* pre-initializing some parameters */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{	
		(*bodies)[i].S22 = 0.0;
		(*bodies)[i].C21 = 0.0;
		(*bodies)[i].S21 = 0.0;
	}

	/* auxiliary variables for reading input */
	char 	*line = NULL;
	size_t 	len = 0;
	ssize_t read;

	/* verification variables for input */
	int 	number_par_inputs = 17;
	bool	input_par_received[number_par_inputs];
	for (int i = 0; i < number_par_inputs; i++)
	{
		input_par_received[i] = false;
	}
	bool	input_orb_received = false;
	bool	input_initial_rot_received = false;
	/* verification variables for names and deformable settings */
	bool	input_name_received = false;
	int 	number_deformable_inputs = 5; 
	bool	input_deformable_received[number_deformable_inputs];
	for (int i = 0; i < number_deformable_inputs; i++)
	{
		input_deformable_received[i] = false;
	}
	/* verification variables for rheology */
	int 	number_maxwell_inputs = 3;
	bool	input_maxwell_received[number_maxwell_inputs];
	for (int i = 0; i < number_maxwell_inputs; i++)
	{
		input_maxwell_received[i] = false;
	}
	int 	number_gV_inputs = 4;
	bool	input_gV_received[number_gV_inputs];
	for (int i = 0; i < number_gV_inputs; i++)
	{
		input_gV_received[i] = false;
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
			input_par_received[1] = true;
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
		else if (strcmp(token, "C21") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].C21 = atof(token);
			}
		}
		else if (strcmp(token, "S21") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].S21 = atof(token);
			}
		}
		else if (strcmp(token, "S22") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].S22 = atof(token);
			}
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
		else if (strcmp(token, "a(AU)") == 0)
		{
			(*bodies)[0].a = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].a = atof(token);
			}
			input_par_received[9] = true;
		}
		else if (strcmp(token, "e") == 0)
		{
			(*bodies)[0].e = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].e = atof(token);
			}
			input_par_received[10] = true;
		}
		else if (strcmp(token, "I(deg)") == 0)
		{
			(*bodies)[0].I = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].I = atof(token);
			}
			input_par_received[11] = true;
		}
		else if (strcmp(token, "M(deg)") == 0)
		{
			(*bodies)[0].M = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].M = atof(token);
			}
			input_par_received[12] = true;
		}
		else if (strcmp(token, "w(deg)") == 0)
		{
			(*bodies)[0].w = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].w = atof(token);
			}
			input_par_received[13] = true;
		}
		else if (strcmp(token, "OMEGA(deg)") == 0)
		{
			(*bodies)[0].Omega = NAN;
			for (int i = 1; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].Omega = atof(token);
			}
			input_par_received[14] = true;
		}
		else if (strcmp(token, "azi(deg)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].azi = atof(token);
			}
			input_par_received[15] = true;
		}
		else if (strcmp(token, "pol(deg)") == 0)
		{
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].pol = atof(token);
			}
			input_par_received[16] = true;
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
			input_deformable_received[3] = true;
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
			input_deformable_received[4] = true;
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
			input_gV_received[3] = true;
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
			// fprintf(stderr, "Missing input number %d\n", i); // debugging
			exit(14);
		}
	}
	if(input_initial_rot_received == false)
	{
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			(*bodies)[i].rot_ini = (*bodies)[i].rot;
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
	/* Maxwell rheology input verification */
	if (strcmp(simulation.rheology_model, "Maxwell") == 0)
	{
		for (int i = 0; i < number_maxwell_inputs; i++)
		{
			if(input_maxwell_received[i] == false)
			{
				fprintf(stderr, "Error: there is at least one missing input ");
				fprintf(stderr, "for the Maxwell rheology variables ");
				fprintf(stderr, "from %s.\n", simulation.system_specs);
				fprintf(stderr, "Exiting the program now.\n");
				exit(14);
			}
		}
	}
	
	/* Generalised Voigt rheology input verification and setting of elements */
	if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
	{
		for (int i = 0; i < number_gV_inputs; i++)
		{
			if(input_gV_received[i] == false)
			{
				fprintf(stderr, "Error: there is at least one missing input ");
				fprintf(stderr, "for the generalised Voigt rheology variables ");
				fprintf(stderr, "from %s.\n", simulation.system_specs);
				fprintf(stderr, "Exiting the program now.\n");
				exit(14);
			}
		}

		// allocate memory for Voigt element
		int elements_max = 0;
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
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
		if ((*bodies)[i].deformable == false &&
			(*bodies)[i].prestress == false)
		{
			(*bodies)[i].point_mass = true;
		}
		if ((*bodies)[i].point_mass == true)
		{
			(*bodies)[i].deformable = false;
			(*bodies)[i].prestress = false;
			(*bodies)[i].centrifugal = false;
			(*bodies)[i].tidal = false;
		}
		else if ((*bodies)[i].deformable == false)
		{
			(*bodies)[i].centrifugal = false;
			(*bodies)[i].tidal = false;
		}
	}

	/* converting units and setting orbital period */

	double deg_to_rad = M_PI / 180.0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		(*bodies)[i].azi *= deg_to_rad;
		(*bodies)[i].pol *= deg_to_rad;
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
			(*bodies)[i].rot *= day_to_s;
			(*bodies)[i].rot_ini *= day_to_s;
			(*bodies)[i].a *= AU_to_m;
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				(*bodies)[i].tau *= year_to_s;
			}
			if (i==0)
			{
				(*bodies)[i].orb = 0.0;	
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
			(*bodies)[i].rot *= day_to_year;
			(*bodies)[i].rot_ini *= day_to_year;
			if (strcmp(simulation.rheology_model, "Maxwell") == 0)
			{
				(*bodies)[i].Dt *= s_to_year;
			}
			if (i==0)
			{
				(*bodies)[i].orb = 0.0;	
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

		/* 2nd set of variables - omega and q */
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

		/* 3rd set of variables - I0 */
		(*bodies)[i].I0 = (3.0 * rg - 2.0 * J2) * m * R * R / 3.0;

		/* Kelvin-Voigt elements */
		if (strcmp(simulation.rheology_model, "Maxwell") == 0)
		{
			(*bodies)[i].elements = 0;
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

			(*bodies)[i].gamma_0 = 1.0;
			(*bodies)[i].alpha 	 = 0.0;
			(*bodies)[i].eta   	 = 1.0;

			if (strcmp(simulation.rheology_model, "gen_Voigt") == 0)
			{
				for (int j = 0; j < (*bodies)[i].elements; j++)
				{
					(*bodies)[i].alpha_elements[j] = 0.0;
					(*bodies)[i].eta_elements[j]   = 1.0;
				}
			}
		}
	} // end loop over bodies

	/* set bs, calculate p, and initialize b_eta and bk */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		/* real deformation on Tisserand frame */
		double B_stokes_i[9];
		body_frame_deformation_from_stokes_coefficients(B_stokes_i, (*bodies)[i]);
		double B_stokes_diag_i[9];
		calculate_diagonalized_square_matrix(B_stokes_diag_i, B_stokes_i);

		/* update gravity field coefficients to */
		/* the frame of principal inertia moments */
		double Iner_diag[9];
		calculate_inertia_tensor(Iner_diag, (*bodies)[i].I0, B_stokes_diag_i);
		(*bodies)[i].rg  = calculate_rg ((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
		(*bodies)[i].J2  = calculate_J2 ((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
		(*bodies)[i].C22 = calculate_C22((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
		(*bodies)[i].S22 = calculate_S22((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
		(*bodies)[i].C21 = calculate_C21((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);
		(*bodies)[i].S21 = calculate_S21((*bodies)[i].mass, (*bodies)[i].R, Iner_diag);

		/* prestress on Tisserand frame */
		double P_i[9];	
		null_matrix(P_i);
		if ((*bodies)[i].prestress == true &&
			(*bodies)[i].deformable == true)
		{
			double F_cent_mean_i[9];
			double mean_omega = 2.0 * M_PI / (*bodies)[i].rot;
			calculate_F_cent_mean(F_cent_mean_i, mean_omega);
			double F_tide_mean_i[9];
			null_matrix(F_tide_mean_i); // no permanent tide for now
			linear_combination_three_square_matrix(P_i,
				(*bodies)[i].gamma_0, B_stokes_diag_i,
				-1.0, F_cent_mean_i,
				-1.0, F_tide_mean_i);
		}

		/* real deformation and prestress on inertial frame */
		double Y_i[9], Y_i_trans[9];
		rotation_matrix_from_quaternion(Y_i, (*bodies)[i].q);
		transpose_square_matrix(Y_i_trans, Y_i);
		double b_stokes_i[9];
		square_matrix_times_square_matrix(b_stokes_i, Y_i, B_stokes_diag_i);
		square_matrix_times_square_matrix(b_stokes_i, b_stokes_i, Y_i_trans);
		get_main_elements_traceless_symmetric_matrix((*bodies)[i].bs_me, b_stokes_i);
		double p_i[9];
		square_matrix_times_square_matrix(p_i, Y_i, P_i);
		square_matrix_times_square_matrix(p_i, p_i, Y_i_trans);			
		get_main_elements_traceless_symmetric_matrix((*bodies)[i].p_me, p_i);

		/* alloc and initialize bk at equilibrium */
		if ((*bodies)[i].elements > 0)
		{
			(*bodies)[i].bk_me 
				= (double *) calloc((*bodies)[i].elements * 5, sizeof(double));
		}

		/* initialize b_eta */
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
			linear_combination_four_square_matrix(b_eta_i,
				 c_i/(*bodies)[i].alpha, b_stokes_i,
				-1.0/(*bodies)[i].alpha, f_cent_i,
				-1.0/(*bodies)[i].alpha, f_tide_i,
				-1.0/(*bodies)[i].alpha, f_ps_i);
		}
		get_main_elements_traceless_symmetric_matrix((*bodies)[i].b_eta_me, b_eta_i);
		
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
create_output_files	(const cltbdy *bodies,
			 		 const siminf simulation,
					 FILE *out[])
{
	char filename[150];
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
	fprintf(out[simulation.number_of_bodies], " x_com(AU)");
	fprintf(out[simulation.number_of_bodies], " y_com(AU)");
	fprintf(out[simulation.number_of_bodies], " z_com(AU)");
	fprintf(out[simulation.number_of_bodies], " |L|");
	fprintf(out[simulation.number_of_bodies], "\n");

	return 0;
}

int
write_output(const cltbdy *bodies,
			 const siminf simulation,
			 FILE *out[])
{
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		/* state variables */

		/* time */
		fprintf (out[i], "%.15e", simulation.t);

		/* position, and velocity */
		if (simulation.keplerian_motion == true)
		{
			fprintf (out[i], " %.15e %.15e %.15e %.15e %.15e %.15e", 
				bodies[i].x[0], bodies[i].x[1], bodies[i].x[2],
				bodies[i].x_dot[0], bodies[i].x_dot[1], bodies[i].x_dot[2]);
		}
		else
		{
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
		}

		if (bodies[i].point_mass == false)
		{
			/* angular velocity vector */
			fprintf (out[i], " %.15e %.15e %.15e", 
				bodies[i].omega[0], bodies[i].omega[1], bodies[i].omega[2]);

			/* angular momentum vector */
			fprintf (out[i], " %.15e %.15e %.15e", 
				bodies[i].l[0], bodies[i].l[1], bodies[i].l[2]);

			/* deformation matrix */
			fprintf (out[i], " %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e",
				bodies[i].b[0], bodies[i].b[1], bodies[i].b[2],
				bodies[i].b[3], bodies[i].b[4], bodies[i].b[5],
				bodies[i].b[6], bodies[i].b[7], bodies[i].b[8]);

			/* quaternion */
			fprintf (out[i], " %.15e %.15e %.15e %.15e", 
				bodies[i].q[0], bodies[i].q[1], bodies[i].q[2], bodies[i].q[3]);

		}

		/* next line */
		fprintf(out[i], "\n");

	} // end loop over bodies

	/* prints on general file */

	/* time */
	fprintf (out[simulation.number_of_bodies], "%.15e", simulation.t);

	/* location of center of mass */
	double center_of_mass[3];
	calculate_center_of_mass(center_of_mass, bodies, simulation.number_of_bodies, simulation.G);
	double relative_center_of_mass[3];
	linear_combination_vector(relative_center_of_mass,
		1.0, center_of_mass,
		-1.0, bodies[0].x);
	fprintf (out[simulation.number_of_bodies], " %.15e", relative_center_of_mass[0]);
	fprintf (out[simulation.number_of_bodies], " %.15e", relative_center_of_mass[1]);
	fprintf (out[simulation.number_of_bodies], " %.15e", relative_center_of_mass[2]);

	/* total angular momentum */
	double l_total[3];
	calculate_total_angular_momentum(l_total, bodies, simulation.number_of_bodies, simulation.G);
	fprintf (out[simulation.number_of_bodies], " %.15e", norm_vector(l_total));

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
			fprintf(out_orientation, " wI3(°)"); // angle w and I3
			fprintf(out_orientation, " wI3sf(°)"); // angle w and I3 solid frame
			fprintf(out_orientation, " wl(°)"); // angle w and l
			fprintf(out_orientation, " azi(°)"); // nutation angle
			fprintf(out_orientation, " aziPIM(°)"); // nutation angle on Principal
													 // Inertia Momenta (PIM) frame
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
				fprintf (out_orientation, "%.15e", t);

				/* norm of angular velocity vector */
				fprintf (out_orientation, " %.15e", norm_vector(bodies[i].omega));

				/* norm of angular momentum vector */
				fprintf (out_orientation, " %.15e", norm_vector(bodies[i].l));

				/* norm of deformation matrix */
				fprintf (out_orientation, " %.15e", norm_square_matrix(bodies[i].b));

				/* obliquity */
				if (i == 0)
				{
					// calculate_obliquity_free_body_from_angular_velocity(&bodies[i]);
					calculate_obliquity_free_body_from_figure_axis_of_solid_frame(&bodies[i]);
				}
				else
				{
					linear_combination_vector(bodies[i].relative_x,
						1.0, bodies[i].x,
						-1.0, bodies[0].x);
					linear_combination_vector(bodies[i].relative_x_dot,
						1.0, bodies[i].x_dot,
						-1.0, bodies[0].x_dot);
					// calculate_obliquity_on_orbit_from_angular_velocity(&bodies[i]);
					calculate_obliquity_on_orbit_from_figure_axis_of_solid_frame(&bodies[i]);
				}	
				fprintf (out_orientation, " %.15e", bodies[i].obl * rad_to_deg);

				/* angle between spin axis and figure axis */
				double angle_wI3;
				angle_wI3 = angle_between_spin_axis_and_figure_axis(bodies[i]);
				fprintf (out_orientation, " %.15e", angle_wI3 * rad_to_deg);

				/* angle between spin axis and figure axis of solid frame */
				double angle_wI3sf;
				angle_wI3sf = angle_between_spin_axis_and_figure_axis_of_solid_frame(bodies[i]);
				fprintf (out_orientation, " %.15e", angle_wI3sf * rad_to_deg);

				/* angle between spin axis and angular momentum */
				double angle_wl;
				angle_wl = angle_between_spin_axis_and_angular_momentum(bodies[i]);
				fprintf (out_orientation, " %.15e", angle_wl * rad_to_deg);

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
				fprintf (out_orientation, " %.15e", bodies[i].azi * rad_to_deg);

				/* nutation frequency on Principal Inertia Momenta (PIM) frame */
				double P_i[9], P_i_trans[9];
				calculate_eigenvectors_matrix(P_i, bodies[i].b);
				transpose_square_matrix(P_i_trans, P_i);	
				square_matrix_times_vector(ang_vel_body_frame,
					P_i_trans, bodies[i].omega);
				cartesian_to_spherical_coordinates(ang_vel_body_frame_spherical,
					ang_vel_body_frame);
				fprintf (out_orientation, " %.15e", 
					ang_vel_body_frame_spherical[1] * rad_to_deg);

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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
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
			fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
			fprintf(gnuplotPipe, "set loadpath \"%s\"\n", simulation.output_folder);
			fprintf(gnuplotPipe, "set output \"%s\"\n", filename_plot);
			fprintf(gnuplotPipe, "set border lw 2 \n");
			fprintf(gnuplotPipe, "unset key\n");
			fprintf(gnuplotPipe, "set xlabel \"Time(yr)\"\n");
			fprintf(gnuplotPipe, "set ylabel \"Azimuthal angle(°)\"\n");
			fprintf(gnuplotPipe, "set title \"Azimuthal angle of %s's angular velocity\"\n", bodies[i].name);
			fprintf(gnuplotPipe, "plot \'%s\' u 1:9 w l lw 3", filename_get);
			pclose(gnuplotPipe);
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
