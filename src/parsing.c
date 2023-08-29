#include "parsing.h"

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
convert_input	(cltbdy	**body,
				 const int number_of_bodies,
				 const double G,
				 const char file[],
				 const char units[])
{
	/* allocate memory for body */
	*body = malloc(number_of_bodies * sizeof(cltbdy));

	/* auxiliary variables for reading input */
    char 	*line = NULL;
    size_t 	len = 0;
    ssize_t read;
	int 	col_num;

	/* verification variables for input */
	int 	number_par_inputs = 18;
	bool	input_par_received[number_par_inputs];
	for (int i = 0; i < number_par_inputs; i++)
	{
		input_par_received[i] = false;
	}
	/* verification variables for names and deformable settings */
	bool	input_name_received = false;
	int 	number_deformable_inputs = 3; 
	bool	input_deformable_received[number_deformable_inputs];
	for (int i = 0; i < number_deformable_inputs; i++)
	{
		input_deformable_received[i] = false;
	}

	/* reading input parameters */
	FILE 	*in1 = fopen(file, "r");
	if	(in1 == NULL)
	{
		fprintf(stderr, "Warning: could not read input file.\n");
		fprintf(stderr, "Exiting the program now.\n");
		exit(13);
	}
   	while ((read = getline(&line, &len, in1)) != -1)
	{
		col_num = count_columns(line);
		if (col_num < number_of_bodies)
		{
			fprintf(stderr, "Warning: number of bodies cannot\n");
			fprintf(stderr, "exceed number of columns in system file.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}
		const char tok_del[3] = " \t\n";		// token delimiter
		char *token = strtok(line, tok_del);
		if (strcmp(token, "Name") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				strcpy((*body)[i].name, token);
			}
			input_name_received = true;
		}
		else if (strcmp(token, "mass(Msun)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].mass = atof(token);
			}
			input_par_received[0] = true;
		}
		else if (strcmp(token, "lod(day)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].lod = atof(token);
			}
			input_par_received[1] = true;
		}
		else if (strcmp(token, "obl(deg)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].obl = atof(token);
			}
			input_par_received[2] = true;
		}
		else if (strcmp(token, "psi(deg)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].psi = atof(token);
			}
			input_par_received[3] = true;
		}
		else if (strcmp(token, "R(km)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].R = atof(token);
			}
			input_par_received[4] = true;
		}
		else if (strcmp(token, "rg") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].rg = atof(token);
			}
			input_par_received[5] = true;
		}
		else if (strcmp(token, "J2") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].J2 = atof(token);
			}
			input_par_received[6] = true;
		}
		else if (strcmp(token, "C22") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].C22 = atof(token);
			}
			input_par_received[7] = true;
		}
		else if (strcmp(token, "lib(deg)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].lib = atof(token);
			}
			input_par_received[8] = true;
		}
		else if (strcmp(token, "kf") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].kf = atof(token);
			}
			input_par_received[9] = true;
		}
		else if (strcmp(token, "Dt(s)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].Dt = atof(token);
			}
			input_par_received[10] = true;
		}
		else if (strcmp(token, "tau(yr)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].tau = atof(token);
			}
			input_par_received[11] = true;
		}
		else if (strcmp(token, "a(AU)") == 0)
		{
			(*body)[0].a = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].a = atof(token);
			}
			input_par_received[12] = true;
		}
		else if (strcmp(token, "e") == 0)
		{
			(*body)[0].e = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].e = atof(token);
			}
			input_par_received[13] = true;
		}
		else if (strcmp(token, "I(deg)") == 0)
		{
			(*body)[0].I = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].I = atof(token);
			}
			input_par_received[14] = true;
		}
		else if (strcmp(token, "M(deg)") == 0)
		{
			(*body)[0].M = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].M = atof(token);
			}
			input_par_received[15] = true;
		}
		else if (strcmp(token, "w(deg)") == 0)
		{
			(*body)[0].w = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].w = atof(token);
			}
			input_par_received[16] = true;
		}
		else if (strcmp(token, "OMEGA(deg)") == 0)
		{
			(*body)[0].Omega = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*body)[i].Omega = atof(token);
			}
			input_par_received[17] = true;
		}
		else if (strcmp(token, "centrifugal") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				if (strcmp(token, "yes") == 0)
				{
					(*body)[i].centrifugal = true;
				}
				else if (strcmp(token, "no") == 0)
				{
					(*body)[i].centrifugal = false;
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
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				if (strcmp(token, "yes") == 0)
				{
					(*body)[i].tidal = true;
				}
				else if (strcmp(token, "no") == 0)
				{
					(*body)[i].tidal = false;
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
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				if (strcmp(token, "yes") == 0)
				{
					(*body)[i].point_mass = true;
				}
				else if (strcmp(token, "no") == 0)
				{
					(*body)[i].point_mass = false;
				}
				else
				{
					printf("%s\n", token);
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
			fprintf(stderr, "from %s.\n", file);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
	}
	/* name  input verification */
	if (input_name_received == false)
	{
		fprintf(stderr, "Error: could not read body name ");
		fprintf(stderr, "from %s.\n", file);
		fprintf(stderr, "Exiting the program now.\n");
		exit(14);
	}
	/* deformable variables input verification */
	for (int i = 0; i < number_deformable_inputs; i++)
	{
		if(input_deformable_received[i] == false)
		{
			fprintf(stderr, "Error: there is at least one missing input ");
			fprintf(stderr, "from %s.\n", file);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
	}

	/* converting units */

	double deg_to_rad = M_PI / 180.0;
	for (int i = 0; i < number_of_bodies; i++)
	{
		(*body)[i].obl *= deg_to_rad;
		(*body)[i].psi *= deg_to_rad;
		(*body)[i].lib *= deg_to_rad;
		(*body)[i].I *= deg_to_rad;
		(*body)[i].M *= deg_to_rad;
		(*body)[i].w *= deg_to_rad;
		(*body)[i].Omega *= deg_to_rad;
	}

	if (strcmp(units, "SI") == 0)
	{
		/* conversion units to SI */
		double Msun_to_kg = 1988500.0e24;
		double day_to_s = 24.0 * 60.0 * 60.0;
		double km_to_m = 1e3;
		double year_to_s = 365.25 * day_to_s;
		double AU_to_m = 1.495978707e11;

		/* variables in SI*/
		for (int i = 0; i < number_of_bodies; i++)
		{
			(*body)[i].mass *= Msun_to_kg;
			(*body)[i].R *= km_to_m;
			(*body)[i].lod *= day_to_s;
			(*body)[i].tau *= year_to_s;
			(*body)[i].a *= AU_to_m;
		}
	}
	else
	{
		/* conversion units to AU Msun year */
		double km_to_AU = 1.0 / 1.495978707e8;
		double day_to_year = 1.0 / 365.25;
		double s_to_year = day_to_year / (24.0 * 60.0 * 60.0);

		/* variables in AU Msun year*/
		for (int i = 0; i < number_of_bodies; i++)
		{
			(*body)[i].R *= km_to_AU;
			(*body)[i].lod *= day_to_year;
			(*body)[i].Dt *= s_to_year;
		}
	}
	
	/* calculating state variables */
	
	for (int i = 0; i < number_of_bodies; i++)
	{
		double	m = (*body)[i].mass;
		double	R = (*body)[i].R;
		double  Td = (*body)[i].lod;

		double	rg = (*body)[i].rg;
		double	J2 = (*body)[i].J2;

		double  theta = (*body)[i].obl;
		double	psi = (*body)[i].psi;
		double	phi = (*body)[i].lib;

		double	a = (*body)[i].a;
		double	e = (*body)[i].e;
		double	I = (*body)[i].I;
		double	M = (*body)[i].M;
		double	w = (*body)[i].w;
		double	Omega = (*body)[i].Omega;

		double	kf = (*body)[i].kf;
		double	Dt = (*body)[i].Dt;
		double	tau = (*body)[i].tau;
		
		/* 1st set of variables - x and x_dot */
		double R_1_I[9];
		double R_3_Omega[9];

		if (i == 0)
		{
			null_vector((*body)[i].x);
			null_vector((*body)[i].x_dot);
			identity_matrix(R_1_I);			// for 2nd set of variables
			identity_matrix(R_3_Omega);		// for 2nd set of variables
		}
		else
		{
			double T = 
			kepler_period((*body)[0].mass, m, G, a);
			double n = (2.0 * M_PI) / T;

			double E = kepler_equation(e, M);
			double r = a * (1.0 - e * cos(E));
			// double f = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) * tan(0.5 * E));
			double f =
				atan2(sqrt(1.0 - e * e) * sin(E), cos(E) - e);

			double position_in_plane[] 
				= {r * cos(f), r * sin(f), 0.0};
			double velocity_in_plane[] 
				= {-1.0 * n * a / sqrt(1.0 - e * e) * sin(f), 
					n * a / sqrt(1.0 - e * e) * (e + cos(f)), 
					0.0};

			double R_3_w[9];
			rotation_matrix_3d_z(R_3_w, w);
			rotation_matrix_3d_x(R_1_I, I);
			rotation_matrix_3d_z(R_3_Omega, Omega);

			double full_rotation_orbit[9];
			square_matrix_times_square_matrix(full_rotation_orbit,
				R_3_Omega, R_1_I);
			square_matrix_times_square_matrix(full_rotation_orbit,
				full_rotation_orbit, R_3_w);

			square_matrix_times_vector((*body)[i].x, full_rotation_orbit, position_in_plane);
			square_matrix_times_vector((*body)[i].x_dot, full_rotation_orbit, velocity_in_plane);
		}

		/* 2nd set of variables - omega and b0_diag */
		// b0_diag not implemented yet. 
		// at the moment, it is being dealt with by the main
		double R_3_psi[9];
		rotation_matrix_3d_z(R_3_psi, psi);
		double R_1_theta[9];
		rotation_matrix_3d_x(R_1_theta, theta);
		double R_3_phi[9];
		rotation_matrix_3d_z(R_3_phi, phi);

		double full_rotation_body[9];
		square_matrix_times_square_matrix(full_rotation_body,
			R_3_Omega, R_1_I);
		square_matrix_times_square_matrix(full_rotation_body,
			full_rotation_body, R_3_psi);
		square_matrix_times_square_matrix(full_rotation_body,
			full_rotation_body, R_1_theta);
		square_matrix_times_square_matrix(full_rotation_body,
			full_rotation_body, R_3_phi);
		
		double omega_on_body[] = {0.0, 0.0, 0.0};
		double omega_direction_on_body[] = {0.0, 0.0, 1.0}; // strong assumption
		scale_vector(omega_on_body, 2.0 * M_PI / Td, omega_direction_on_body);
		square_matrix_times_vector((*body)[i].omega, full_rotation_body, omega_on_body);

		/* 3rd set of variables - I0 */
		(*body)[i].I0 = (3.0 * rg - 2.0 * J2) * m * R * R / 3.0;

		/* 4th set of variables - alpha and eta */
		(*body)[i].gamma = parameter_gamma(G, (*body)[i].I0, R, kf);
		(*body)[i].alpha = (*body)[i].gamma * Dt / (tau - Dt);
		(*body)[i].eta = (*body)[i].gamma * Dt;

		/* Kelvin-Voigt elements */
		(*body)[i].elements = 0;

		if ((*body)[i].elements > 0)
		{
			(*body)[i].alpha_elements = malloc((*body)[i].elements * sizeof(double));
			(*body)[i].eta_elements = malloc((*body)[i].elements * sizeof(double));
			
			for (int j = 0; j < (*body)[i].elements; j++)
			{
				(*body)[i].alpha_elements[j] = 1.0;
				(*body)[i].eta_elements[j] = 1.0;
			}
		}

	} // end loop over bodies

	/* for testing */
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%s ", (*body)[i].name);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].mass);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].lod);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].obl);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].psi);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].R);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].rg);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].J2);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].C22);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].lib);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].kf);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].Dt);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].tau);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].a);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].e);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].I);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].M);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].w);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", (*body)[i].Omega);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", (*body)[i].tidal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", (*body)[i].centrifugal);
	// printf("\n");
	// exit(99);

	return 0;
}