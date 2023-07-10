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
			else if (strcmp(first_col, "k2") == 0)
			{
				for (int i = 0; i < *number_of_bodies; i++)
				{
					(*body)[i].k2 = other_col[i];
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

	/* auxiliary variables for input names */
	char	other_col_names[*number_of_bodies][100];

	/* reading input names */
	bool	input_name_received = false; // verification variable
	FILE 	*in1_names = fopen(file, "r");
   	while ((read = getline(&line, &len, in1)) != -1)
	{
		sscanf(line, "%s %s %s %s",
			first_col, other_col_names[0], 
			other_col_names[1], other_col_names[2]);
		if (strcmp(first_col, "Name") == 0)
		{
			for (int i = 0; i < *number_of_bodies; i++)
			{
				strcpy((*body)[i].name, other_col_names[i]);
			}
			input_name_received = true;
		}
	}
	fclose(in1_names);
	if (input_name_received == false)
	{
		fprintf(stderr, "Error: could not read body name ");
		fprintf(stderr, "from %s.\n", file);
		fprintf(stderr, "Exiting the program now.\n");
		exit(14);
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
				 double *kf, double b0_diag[], double omega[],
				 double *alpha, double *eta,
				 double tilde_x[], double tilde_x_dot[],
				 const double G,
				 const char file[])
{
	/* array of members of structure CelestialBody */
	int		number_of_bodies;
	cltbdy	*body;

	/* get values from input file */
	read_input(&body, &number_of_bodies, file);

	/* implementing */
	// *m1 = body[0].mass;
	// *m2 = body[1].mass;
	// *R	= body[0].R;
	// double e = body[1].e;
	// double a = body[1].a;
	// double T = kepler_period(*m1, *m2, G, a);

	// double E = kepler_equation(e, t);
    // double r = a * (1.0 - e * cos(E));
    // double f = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
	// 		* tan(0.5 * E));

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
	// 	printf("%1.5e ", body[i].k2);
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