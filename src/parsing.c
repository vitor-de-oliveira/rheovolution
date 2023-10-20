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
convert_input	(cltbdy	**bodies,
				 const int number_of_bodies,
				 const double G,
				 const char file[],
				 const char units[])
{
	/* allocate memory for bodies */
	*bodies = (cltbdy *) malloc(number_of_bodies * sizeof(cltbdy));

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
	FILE 	*in1 = fopen(file, "r");
   	while ((read = getline(&line, &len, in1)) != -1)
	{
		const char tok_del[6] = " \t\n";	// token delimiter
		char *token = strtok(line, tok_del);
		if(token == NULL) break;	// in case there is a newline
		if (strcmp(token, "Name") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				strcpy((*bodies)[i].name, token);
			}
			input_name_received = true;
		}
		else if (strcmp(token, "mass(Msun)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].mass = atof(token);
			}
			input_par_received[0] = true;
		}
		else if (strcmp(token, "lod(day)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].lod = atof(token);
			}
			input_par_received[1] = true;
		}
		else if (strcmp(token, "obl(deg)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].obl = atof(token);
			}
			input_par_received[2] = true;
		}
		else if (strcmp(token, "psi(deg)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].psi = atof(token);
			}
			input_par_received[3] = true;
		}
		else if (strcmp(token, "R(km)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].R = atof(token);
			}
			input_par_received[4] = true;
		}
		else if (strcmp(token, "rg") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].rg = atof(token);
			}
			input_par_received[5] = true;
		}
		else if (strcmp(token, "J2") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].J2 = atof(token);
			}
			input_par_received[6] = true;
		}
		else if (strcmp(token, "C22") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].C22 = atof(token);
			}
			input_par_received[7] = true;
		}
		else if (strcmp(token, "lib(deg)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].lib = atof(token);
			}
			input_par_received[8] = true;
		}
		else if (strcmp(token, "kf") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].kf = atof(token);
			}
			input_par_received[9] = true;
		}
		else if (strcmp(token, "Dt(s)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].Dt = atof(token);
			}
			input_par_received[10] = true;
		}
		else if (strcmp(token, "tau(yr)") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].tau = atof(token);
			}
			input_par_received[11] = true;
		}
		else if (strcmp(token, "a(AU)") == 0)
		{
			(*bodies)[0].a = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].a = atof(token);
			}
			input_par_received[12] = true;
		}
		else if (strcmp(token, "e") == 0)
		{
			(*bodies)[0].e = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].e = atof(token);
			}
			input_par_received[13] = true;
		}
		else if (strcmp(token, "I(deg)") == 0)
		{
			(*bodies)[0].I = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].I = atof(token);
			}
			input_par_received[14] = true;
		}
		else if (strcmp(token, "M(deg)") == 0)
		{
			(*bodies)[0].M = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].M = atof(token);
			}
			input_par_received[15] = true;
		}
		else if (strcmp(token, "w(deg)") == 0)
		{
			(*bodies)[0].w = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].w = atof(token);
			}
			input_par_received[16] = true;
		}
		else if (strcmp(token, "OMEGA(deg)") == 0)
		{
			(*bodies)[0].Omega = NAN;
			for (int i = 1; i < number_of_bodies; i++)
			{
				token = strtok(NULL, tok_del);
				(*bodies)[i].Omega = atof(token);
			}
			input_par_received[17] = true;
		}
		else if (strcmp(token, "keplerian") == 0)
		{
			for (int i = 0; i < number_of_bodies; i++)
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
			for (int i = 0; i < number_of_bodies; i++)
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
			for (int i = 0; i < number_of_bodies; i++)
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
			for (int i = 0; i < number_of_bodies; i++)
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
			for (int i = 0; i < number_of_bodies; i++)
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
			fprintf(stderr, "from %s.\n", file);
			fprintf(stderr, "Exiting the program now.\n");
			// fprintf(stderr, "Missing input number %d\n", i); // for testing
			exit(14);
		}
	}
	/* name  input verification */
	if (input_name_received == false)
	{
		fprintf(stderr, "Error: could not read body names ");
		fprintf(stderr, "from %s.\n", file);
		fprintf(stderr, "Exiting the program now.\n");
		exit(14);
	}
	/* orbital variable input verification */
	if (input_keplerian_received == false)
	{
		fprintf(stderr, "Error: missing keplerian status ");
		fprintf(stderr, "from %s.\n", file);
		fprintf(stderr, "Exiting the program now.\n");
		exit(14);
	}
	if (input_orbit_2body_received == false)
	{
		fprintf(stderr, "Error: missing orbit 2 body status ");
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
			fprintf(stderr, "for the deformation variables ");
			fprintf(stderr, "from %s.\n", file);
			fprintf(stderr, "Exiting the program now.\n");
			exit(14);
		}
	}

	/* converting units */

	double deg_to_rad = M_PI / 180.0;
	for (int i = 0; i < number_of_bodies; i++)
	{
		(*bodies)[i].obl *= deg_to_rad;
		(*bodies)[i].psi *= deg_to_rad;
		(*bodies)[i].lib *= deg_to_rad;
		(*bodies)[i].I *= deg_to_rad;
		(*bodies)[i].M *= deg_to_rad;
		(*bodies)[i].w *= deg_to_rad;
		(*bodies)[i].Omega *= deg_to_rad;
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
		for (int i = 0; i < number_of_bodies; i++)
		{
			(*bodies)[i].R *= km_to_AU;
			(*bodies)[i].lod *= day_to_year;
			(*bodies)[i].Dt *= s_to_year;
		}
	}
	
	/* calculating state variables */
	
	for (int i = 0; i < number_of_bodies; i++)
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
			kepler_period((*bodies)[0].mass, m, G, a);
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
		(*bodies)[i].gamma = parameter_gamma(G, (*bodies)[i].I0, R, kf);
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

	return 0;
}
