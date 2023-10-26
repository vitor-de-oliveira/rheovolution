#define _GNU_SOURCE // getline
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h> // ssize_t
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "celestial_mechanics.h"
#include "dynamical_system.h"
#include "data_processing.h"

#define t(n) printf("Here %d\n", n) // for testing

int
main(int argc, char *argv[]) 
{
	/* start clock */
	clock_t begin_time = clock(), end_time;

	/* define and get simulation info */
	siminf	simulation;
	parse_input(&simulation, argv[1]);

	/* define and fill in array of celestial bodies */
	cltbdy	*bodies;
	fill_in_bodies_data(&bodies, simulation);

	/* other software commands */
	if (argv[2] != NULL)
	{
		if (strcmp(argv[2], "orbit") == 0)
		{
			/* calculate orbital elements from output */
			output_to_orbit(bodies, simulation);
			return 0;
		}
		else if (strcmp(argv[2], "spin") == 0)
		{
			/* calculate spin variables from output */
			output_to_spin(bodies, simulation);
			return 0;		
		}
		else if (strcmp(argv[2], "plot") == 0)
		{
			/* plot everything */
			plot_output_comma_orbit_and_spin(bodies, simulation);
			return 0;
		}
		else
		{
			fprintf(stderr, "Warning: could not recognize 2nd argument of main.\n");
			fprintf(stderr, "Exiting the program now.\n");
			exit(13);
		}
	}

	/* complete b0_me */
	// I am using same b0_diag and alpha_0 for every body for now!
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		bodies[i].b0_me[0] = 0.0;
		bodies[i].b0_me[1] = 0.0;
		bodies[i].b0_me[2] = 0.0;
		bodies[i].b0_me[3] = 0.0;
		bodies[i].b0_me[4] = 0.0;
		bodies[i].alpha_0 = 0.0;
	}

	/* variables not given by user */
	for (int i  = 0; i < simulation.number_of_bodies; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			bodies[i].u_me[j] = 0.01;
		}
		if (bodies[i].elements > 0)
		{
			bodies[i].bk_me = (double *) calloc(bodies[i].elements * 5, sizeof(double));
		}
	}

	/* calculate first l */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		calculate_b(i, bodies, simulation.number_of_bodies, simulation.G);
		calculate_l(&bodies[i], simulation.number_of_bodies, simulation.G);
	}

	/* total number of Voigt elements */
	int	elements_total = 0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		elements_total += bodies[i].elements;
	}
	
	/* state variables */
	int		dim_state_per_body_without_elements = 23;
	int		dim_state = (dim_state_per_body_without_elements * simulation.number_of_bodies) + (5 * elements_total);
	double 	y[dim_state];
	int		elements_counter = 0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
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
		for (int j = 0; j < 4; j++)
		{
			y[19 + 5 * bodies[i].elements + dim_state_skip + j]	= bodies[i].q[j];
		}
		elements_counter += bodies[i].elements;
	}

	/* parameters */
	int		dim_params_per_body_without_elements = 15;
	int		dim_params = 2 + (dim_params_per_body_without_elements * simulation.number_of_bodies) + (2 * elements_total);
	double	params[dim_params];
	elements_counter = 0; 
	params[0] = simulation.G;
	params[1] = (double) simulation.number_of_bodies;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		params[2 + 0 + dim_params_skip] = (double) bodies[i].keplerian;
		params[2 + 1 + dim_params_skip] = (double) bodies[i].orbit_2body;
		params[2 + 2 + dim_params_skip] = (double) bodies[i].point_mass;
		params[2 + 3 + dim_params_skip] = (double) bodies[i].centrifugal;
		params[2 + 4 + dim_params_skip] = (double) bodies[i].tidal;
		for (int j = 0; j < 3; j++)
		{
			params[2 + 5 + dim_params_skip + j] = bodies[i].omega[j];
		}
		params[2 + 8 + dim_params_skip] = bodies[i].mass;
		params[2 + 9 + dim_params_skip] = bodies[i].I0;
		params[2 + 10 + dim_params_skip] = bodies[i].gamma;
		params[2 + 11 + dim_params_skip] = bodies[i].alpha;
		params[2 + 12 + dim_params_skip] = bodies[i].eta;
		params[2 + 13 + dim_params_skip] = bodies[i].alpha_0;
		params[2 + 14 + dim_params_skip] = (double) bodies[i].elements;
		for (int j = 0; j < bodies[i].elements; j++)
		{
			params[2 + 15 + dim_params_skip + j] = bodies[i].alpha_elements[j];
		}
		for (int j = 0; j < bodies[i].elements; j++)
		{
			params[2 + 16 + dim_params_skip + j + bodies[i].elements - 1] = bodies[i].eta_elements[j];
		}
		elements_counter += bodies[i].elements;		
	}

	/* GSL variables */
	const gsl_odeiv2_step_type * ode_type
		= gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step * ode_step
    	= gsl_odeiv2_step_alloc (ode_type, dim_state);
  	gsl_odeiv2_control * ode_control
    	= gsl_odeiv2_control_y_new (simulation.eps_abs, simulation.eps_rel);
  	gsl_odeiv2_evolve * ode_evolve
    	= gsl_odeiv2_evolve_alloc (dim_state);

	/* create output files */
	FILE *out[simulation.number_of_bodies + 1];
	create_output_files(bodies, simulation, out);

	/* integration loop */
	simulation.counter = 0;	
	simulation.t = simulation.t_init;
	while (simulation.t < simulation.t_final)
	{
		if (fabs(simulation.t_final-simulation.t) < simulation.t_step)
		{
			simulation.t_step = simulation.t_final - simulation.t; // smaller last step
		}

	  	gsl_odeiv2_system sys = {field_GV, NULL, dim_state, params};
	
		int status = 
			gsl_odeiv2_evolve_apply_fixed_step (ode_evolve, 
				ode_control, ode_step, &sys, &simulation.t, simulation.t_step, y);

		if (status != GSL_SUCCESS)
		{
			fprintf(stderr, "Error: GSL odeiv2 status = %d\n", status);
			break;
		}

		/* update variables */
		elements_counter = 0; 
		for (int i = 0; i < simulation.number_of_bodies; i++)
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
			for (int j = 0; j < 4; j++)
			{
				bodies[i].q[j] = y[19 + 5 * bodies[i].elements + dim_state_skip + j];
			}

			elements_counter += bodies[i].elements;
		}

		/* update omega */
		elements_counter = 0;
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		
			calculate_omega(i, bodies, simulation.number_of_bodies, simulation.G);
		
			for (int j = 0; j < 3; j++)
			{
				params[2 + 5 + dim_params_skip + j] = bodies[i].omega[j];
			}

			elements_counter += bodies[i].elements;
		}

		/* calculate b */
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			calculate_b(i, bodies, simulation.number_of_bodies, simulation.G);
		}

		/* write output */
		if (simulation.t > simulation.t_trans)
		{
			if (simulation.counter % simulation.data_skip == 0)
			{
				write_output(bodies, simulation, out);
			}
			simulation.counter++;
		}
	}

	/* close output files */
	close_output_files(simulation, out);

	/* free GSL variables */
	gsl_odeiv2_evolve_free (ode_evolve);
	gsl_odeiv2_control_free (ode_control);
	gsl_odeiv2_step_free (ode_step);

	/* free Voigt elements */
	for (int i = 0; i < simulation.number_of_bodies; i++)
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

	/* stop clock */
	end_time = clock();
	int time_spent_in_seconds 
		= (end_time - begin_time) / CLOCKS_PER_SEC;

	/* write overview file */
	write_simulation_overview(time_spent_in_seconds, simulation);
	
	return 0;
}
