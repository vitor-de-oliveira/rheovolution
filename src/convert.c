#include "convert.h"

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
read_input(cltbdy **body, int *number_of_bodies, const char file[])
{
	*number_of_bodies = 3; 	// only works for 3 for now
							// because of sscanf
							// if this changes the test
							// area has to change as well

	/* allocate memory for body */
	*body = malloc(*number_of_bodies * sizeof(cltbdy));

	/* auxiliary variables for reading input */
    char 	*line = NULL;
    size_t 	len = 0;
    ssize_t read;
	int 	col_num;
	char 	first_col[100];
	double 	other_col[*number_of_bodies];

	/* verification variables for input */
	int 	number_par_inputs = 18;
	bool	input_par_received[number_par_inputs];
	for (int i = 0; i < number_par_inputs; i++)
	{
		input_par_received[i] = false;
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
		if (col_num == 4)
		{
			sscanf(line, "%s %lf %lf %lf",
				first_col, &other_col[0], &other_col[1], &other_col[2]);
			if (strcmp(first_col, "mass(Msun)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].mass = other_col[i];
				}
				input_par_received[0] = true;
			}
			else if (strcmp(first_col, "lod(day)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].lod = other_col[i];
				}
				input_par_received[1] = true;
			}
			else if (strcmp(first_col, "obl(deg)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].obl = other_col[i];
				}
				input_par_received[2] = true;
			}
			else if (strcmp(first_col, "psi(deg)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].psi = other_col[i];
				}
				input_par_received[3] = true;
			}
			else if (strcmp(first_col, "R(km)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].R = other_col[i];
				}
				input_par_received[4] = true;
			}
			else if (strcmp(first_col, "rg") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].rg = other_col[i];
				}
				input_par_received[5] = true;
			}
			else if (strcmp(first_col, "J2") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].J2 = other_col[i];
				}
				input_par_received[6] = true;
			}
			else if (strcmp(first_col, "C22") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].C22 = other_col[i];
				}
				input_par_received[7] = true;
			}
			else if (strcmp(first_col, "lib(deg)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].lib = other_col[i];
				}
				input_par_received[8] = true;
			}
			else if (strcmp(first_col, "kf") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].kf = other_col[i];
				}
				input_par_received[9] = true;
			}
			else if (strcmp(first_col, "Dt(s)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].Dt = other_col[i];
				}
				input_par_received[10] = true;
			}
			else if (strcmp(first_col, "tau(yr)") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].tau = other_col[i];
				}
				input_par_received[11] = true;
			}
		}
		if (col_num == 3)
		{
			sscanf(line, "%s %lf %lff",
				first_col, &other_col[0], &other_col[1]);
			if (strcmp(first_col, "a(AU)") == 0)
			{
				(*body)[0].a = NAN;
				for (int i = 1; i < *number_of_bodies; i++)
				{
					(*body)[i].a = other_col[i-1];
				}
				input_par_received[12] = true;
			}
			else if (strcmp(first_col, "e") == 0)
			{
				(*body)[0].e = NAN;
				for (int i = 1; i < *number_of_bodies; i++)
				{
					(*body)[i].e = other_col[i-1];
				}
				input_par_received[13] = true;
			}
			else if (strcmp(first_col, "I(deg)") == 0)
			{
				(*body)[0].I = NAN;
				for (int i = 1; i < *number_of_bodies; i++)
				{
					(*body)[i].I = other_col[i-1];
				}
				input_par_received[14] = true;
			}
			else if (strcmp(first_col, "M(deg)") == 0)
			{
				(*body)[0].M = NAN;
				for (int i = 1; i < *number_of_bodies; i++)
				{
					(*body)[i].M = other_col[i-1];
				}
				input_par_received[15] = true;
			}
			else if (strcmp(first_col, "w(deg)") == 0)
			{
				(*body)[0].w = NAN;
				for (int i = 1; i < *number_of_bodies; i++)
				{
					(*body)[i].w = other_col[i-1];
				}
				input_par_received[16] = true;
			}
			else if (strcmp(first_col, "OMEGA(deg)") == 0)
			{
				(*body)[0].Omega = NAN;
				for (int i = 1; i < *number_of_bodies; i++)
				{
					(*body)[i].Omega = other_col[i-1];
				}
				input_par_received[17] = true;
			}
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

	/* auxiliary variables for input names and deformable settings */
	char	other_col_str[*number_of_bodies][100];

	/* verification variables for names and deformable settings */
	bool	input_name_received = false;
	int 	number_deformable_inputs = 2; 
	bool	input_deformable_received[number_deformable_inputs];
	for (int i = 0; i < number_deformable_inputs; i++)
	{
		input_deformable_received[i] = false;
	}

	/* reading input names and deformable settings */
	FILE 	*in1_names = fopen(file, "r");
   	while ((read = getline(&line, &len, in1)) != -1)
	{
		sscanf(line, "%s %s %s %s",
			first_col, other_col_str[0], 
			other_col_str[1], other_col_str[2]);
		if (strcmp(first_col, "Name") == 0)
		{
			for (int i = 0; i < *number_of_bodies; i++)
			{
				strcpy((*body)[i].name, other_col_str[i]);
			}
			input_name_received = true;
		}
		else if (strcmp(first_col, "centrifugal") == 0)
		{
			for (int i = 0; i < *number_of_bodies; i++)
			{
				if (strcmp(other_col_str[i], "yes") == 0)
				{
					(*body)[i].centrifugal = true;
				}
				else if (strcmp(other_col_str[i], "no") == 0)
				{
					(*body)[i].centrifugal = false;
				}
				else
				{
					fprintf(stderr, "Please provide yes or no ");
					fprintf(stderr, "for centrifugal variable\n");
					exit(14);
				}
			}
			input_deformable_received[0] = true;
		}
		else if (strcmp(first_col, "tidal") == 0)
		{
			for (int i = 0; i < *number_of_bodies; i++)
			{
				if (strcmp(other_col_str[i], "yes") == 0)
				{
					(*body)[i].tidal = true;
				}
				else if (strcmp(other_col_str[i], "no") == 0)
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
	}

	/* name  input verification */
	fclose(in1_names);
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

	/* make a copy of the input file */
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}
	FILE *in_to_copy = fopen(file, "r");
	FILE *in_copy = fopen("output/input_sys_copy.txt" , "w");
	char ch = fgetc(in_to_copy);
    while(ch != EOF)
    {
        fputc(ch, in_copy);
        ch = fgetc(in_to_copy);
    }
	fclose(in_to_copy);
	fclose(in_copy);

	return 0;
}

int
convert_input	(double *m1, double *m2, double *I0, double *R,
				 double *kf, double omega[],
				 double *alpha, double *eta,
				 double tilde_x[], double tilde_x_dot[],
				 bool *centrifugal, bool *tidal,
				 const double G,
				 const char file[])
{
	/* array of members of structure CelestialBody */
	int		number_of_bodies;
	cltbdy	*body;

	/* get values from input file */
	read_input(&body, &number_of_bodies, file);

	/* deformable setting variables */
	*centrifugal = body[0].centrifugal;
	*tidal = body[0].tidal;

	/* conversion angles to rad */
	double deg = M_PI / 180.0; // rad

	/* conversion units to SI */
	// double Msun = 1988500.0e24; // kg
	// double day = 24 * 60 * 60; // s
	// double km = 1e3; // m
	// double year = 365.25 * day; // s
	// double AU = 1.495978707e11; // m

	// /* variables in SI*/
	// *m1 = body[0].mass * Msun;
	// *m2 = body[1].mass * Msun;
	// *R	= body[0].R * km;
	// *kf = body[0].kf;

	// double Td = body[0].lod * day;
	// double theta = body[0].obl * deg;
	// double psi = body[0].psi * deg;
	// double rg = body[0].rg;
	// double phi = body[0].lib * deg;
	// double Dt = body[0].Dt;
	// double tau = body[0].tau * year;
	// double a = body[1].a * AU;
	// double e = body[1].e;
	// double I = body[1].I * deg;
	// double M = body[1].M * deg;
	// double w = body[1].w * deg;
	// double Omega = body[1].Omega * deg;

	/* conversion units to AU Msun year */
	double km_AU = 1.0 / 1.495978707e8; // AU
	double day_year = 1.0 / 365.25; // year
	double s_year = 1.0 / (365.25 * 24 * 60 * 60); // year

	/* variables in AU Msun year*/
	*m1 = body[0].mass;
	*m2 = body[1].mass;
	*R	= body[0].R * km_AU;
	*kf = body[0].kf;

	double Td = body[0].lod * day_year;
	double theta = body[0].obl * deg;
	double psi = body[0].psi * deg;
	double rg = body[0].rg;
	double phi = body[0].lib * deg;
	double Dt = body[0].Dt * s_year;
	double tau = body[0].tau;
	double a = body[1].a;
	double e = body[1].e;
	double I = body[1].I * deg;
	double M = body[1].M * deg;
	double w = body[1].w * deg;
	double Omega = body[1].Omega * deg;

	/* auxiliary variables */
	double T = kepler_period(*m1, *m2, G, a);
	// double T = kepler_period_only_m1(*m1, G, a);
	double n = (2.0 * M_PI) / T;

	/* 1st set of variables - tilde_x and tilde_x_dot */
	double E = kepler_equation(e, M);
    double r = a * (1.0 - e * cos(E));
    // double f = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
	// 		* tan(0.5 * E));
	double f = atan2(sqrt(1.0 - e * e) * sin(E), cos(E) - e);

	double position_in_plane[] 
		= {r * cos(f), r * sin(f), 0.0};
	double velocity_in_plane[] 
		= {-1.0 * n * a / sqrt(1.0 - e * e) * sin(f), 
			n * a / sqrt(1.0 - e * e) * (e + cos(f)), 
			0.0};

	double R_3_w[9];
	rotation_matrix_3d_z(R_3_w, w);
	double R_1_I[9];
	rotation_matrix_3d_x(R_1_I, I);
	double R_3_Omega[9];
	rotation_matrix_3d_z(R_3_Omega, Omega);

	double full_rotation_orbit[9];
	square_matrix_times_square_matrix(full_rotation_orbit,
		R_3_Omega, R_1_I);
	square_matrix_times_square_matrix(full_rotation_orbit,
		full_rotation_orbit, R_3_w);

	square_matrix_times_vector(tilde_x, full_rotation_orbit, position_in_plane);
	square_matrix_times_vector(tilde_x_dot, full_rotation_orbit, velocity_in_plane);

	/* for testing */
	// print_vector(tilde_x);
	// print_vector(tilde_x_dot);

	/* 2nd set of variables - omega and b0_diag */
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
	square_matrix_times_vector(omega, full_rotation_body, omega_on_body);

	/* for testing */
	// print_vector(omega_on_body);
	// printf("%.5e\n", norm_vector(omega_on_body));
	// print_vector(omega);
	// printf("%.5e\n", norm_vector(omega));
	// printf("Omega = %e I = %e psi = %e theta = %e phi = %e\n",
	// 	Omega, I, psi, theta, phi);
	// double omega_test[] = {0.0, 0.0, 1.0};
	// copy_vector(omega_test, omega_on_body);
	// print_vector(omega_test);
	// printf("%.5e\n", norm_vector(omega_test));
	// square_matrix_times_vector(omega_test, R_3_phi, omega_test);
	// print_vector(omega_test);
	// printf("%.5e\n", norm_vector(omega_test));
	// print_square_matrix(R_1_theta);
	// double omega_test_copy[] = {0.0, 0.0, 1.0};
	// copy_vector(omega_test_copy, omega_test);
	// square_matrix_times_vector(omega_test, R_1_theta, omega_test);
	// print_vector(omega_test);
	// printf("%.5e\n", norm_vector(omega_test));
	// square_matrix_times_vector(omega_test, R_3_psi, omega_test);
	// print_vector(omega_test);
	// printf("%.5e\n", norm_vector(omega_test));
	// square_matrix_times_vector(omega_test, R_1_I, omega_test);
	// print_vector(omega_test);
	// printf("%.5e\n", norm_vector(omega_test));
	// square_matrix_times_vector(omega_test, R_3_Omega, omega_test);
	// print_vector(omega_test);
	// printf("%.5e\n", norm_vector(omega_test));
	// exit(14);

	/* 3rd set of variables - I0 */
	*I0 = rg * (*m1) * (*R) * (*R);
	// double J2 = body[0].J2;
	// *I0 = (3.0 * rg - 2.0 * J2) * (*m1) * (*R) * (*R) / 3.0;

	/* for testing */
	// printf("%1.5e %1.5e %1.5e\n", rg, *m1, *R);
	// printf("%1.5e\n", *I0);
	// exit(99);

	/* 4th set of variables - alpha and eta */
	double gamma = 3.0 * (*I0) * G / (pow((*R), 5.0) * (*kf));
	*alpha = gamma * Dt / (tau - Dt);
	*eta = gamma * Dt;

	/* for testing */
	// printf("%1.5e\n", *alpha);
	// printf("%1.5e\n", *eta);

	/* for testing */
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%s ", body[i].name);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].mass);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].lod);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].obl);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].psi);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].R);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].rg);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].J2);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].C22);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].lib);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].kf);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].Dt);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].tau);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].a);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].e);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].I);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].M);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].w);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", body[i].Omega);
	// printf("\n");

	free(body);

	return 0;
}