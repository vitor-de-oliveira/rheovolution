#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "field.h"

int
main(int argc, char *argv[]) 
{
	// check number of input files
	if( argc < 3 )
	{
	    	printf("Please, provide both the system_pars and the integrator_pars files\n");
	    	exit(2);
   	}
   	else if ( argc > 3 )
   	{
	    	printf("Too many arguments\n");
	    	exit(4);
   	}

	// auxiliary variables for fscanf
	char var_name[100];
	double var_value;

	// system parameters
	double mu = mu;
	double x0 = x0;
	double y0 = y0;
	double z0 = z0;
	double x0_dot = x0_dot;
	double y0_dot = y0_dot;
	double z0_dot = z0_dot;
	FILE *in1 = fopen(argv[1], "r");
	while(fscanf(in1, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "mu") == 0)				mu = var_value;
		else if (strcmp(var_name, "x0") == 0)			x0 = var_value;
		else if (strcmp(var_name, "y0") == 0)			y0 = var_value;
		else if (strcmp(var_name, "z0") == 0)			z0 = var_value;
		else if (strcmp(var_name, "x0_dot") == 0)		x0_dot = var_value;
		else if (strcmp(var_name, "y0_dot") == 0)		y0_dot = var_value;
		else if (strcmp(var_name, "z0_dot") == 0)		z0_dot = var_value;
	}
	fclose(in1);

	// integrator parameters
	int dim = dim;
	double h = h;
	double t0 = t0;
	double t1 = t1;
	double eps_abs = eps_abs;
	double eps_rel = eps_rel;
	int data_step = data_step;
	FILE *in2 = fopen(argv[2], "r");
	while(fscanf(in2, " %99[^' '] = %lf[^\n]", var_name, &var_value) != EOF)
	{
		if (strcmp(var_name, "dim") == 0)				dim = (int) var_value;
		else if (strcmp(var_name, "h") == 0)			h = var_value;
		else if (strcmp(var_name, "t0") == 0)			t0 = var_value;
		else if (strcmp(var_name, "t1") == 0)			t1 = var_value;
		else if (strcmp(var_name, "eps_abs") == 0)		eps_abs = var_value;
		else if (strcmp(var_name, "eps_rel") == 0)		eps_rel = var_value;
		else if (strcmp(var_name, "data_step") == 0)	data_step = (int) var_value;
	}
	fclose(in2);

	// GSL variables
	const gsl_odeiv2_step_type * T
		= gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * s
    	= gsl_odeiv2_step_alloc (T, dim);
  	gsl_odeiv2_control * c
    	= gsl_odeiv2_control_y_new (eps_abs, eps_rel);
  	gsl_odeiv2_evolve * e
    	= gsl_odeiv2_evolve_alloc (dim);
  	gsl_odeiv2_system sys = {func, NULL, dim, &mu};

	// integration loop variables
	double t = t0;
	double y[6];
	y[0] = x0;
	y[1] = y0;
	y[2] = z0;
	y[3] = x0_dot;
	y[4] = y0_dot;
	y[5] = z0_dot;

	// integration loop
	int counter = 0;
	while (t < t1)
	{
		int status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &t, h, y);
		if (status != GSL_SUCCESS) break;
		if (counter % data_step == 0)
		{
			printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
		}
		counter++;
	}

	// free GSL variables
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);
	
	return 0;
}
