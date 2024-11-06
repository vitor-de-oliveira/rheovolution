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
			/* plot everything using Gnuplot */
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

	/* correction for the angular velocity's initial value */
	int		omega_correction_number_of_iterates = 10; // 5
	int		*omega_correction_on_body = *(&omega_correction_on_body);
	double	*omega_correction_step = *(&omega_correction_step);
	double 	*omega_correction_rot = *(&omega_correction_rot);
	double	*omega_correction_rot_after_simulation = *(&omega_correction_rot_after_simulation);
	siminf	simulation_copy = *(&simulation_copy);
	cltbdy	*bodies_copy = *(&bodies_copy);
	if (simulation.omega_correction == true)
	{
		omega_correction_on_body = (int *) malloc(simulation.number_of_bodies * sizeof(int));
		omega_correction_step = (double *) malloc(simulation.number_of_bodies * sizeof(double));
		omega_correction_rot = (double *) malloc(simulation.number_of_bodies * sizeof(double));
		omega_correction_rot_after_simulation = (double *) malloc(simulation.number_of_bodies * sizeof(double));
		bodies_copy = (cltbdy *) malloc(simulation.number_of_bodies * sizeof(cltbdy));
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			bodies_copy[i] = create_and_copy_CelestialBody(bodies[i]);
			if (bodies[i].point_mass == false &&
				bodies[i].deformable == true)
			{
				omega_correction_on_body[i] = 1;
			}
			else
			{
				omega_correction_on_body[i] = 0;
			}
			omega_correction_rot[i] = bodies[i].rot;
		}
		// simulation.omega_correction_t_final = 
		// 	10.0 * find_largest_time_scale(bodies, simulation):
		// printf("%1.5e\n", simulation.omega_correction_t_final);
		simulation.omega_correction_t_final = 100.0; // 50000.0
		simulation.omega_correction_counter = 0;
		simulation_copy = simulation;
		simulation.write_to_file = false;
		simulation.keplerian_motion = true;
	}

	back_omega_correction:;

	if (simulation.omega_correction == true)
	{
		if (simulation.omega_correction_counter > 0) // already run once
		{
			simulation.t_step = simulation_copy.t_step;
			if (simulation.omega_correction_counter == omega_correction_number_of_iterates)
			{
				simulation.omega_correction = false;
				simulation.write_to_file = true;
				simulation.keplerian_motion = simulation_copy.keplerian_motion;
			}
			for (int i = 0; i < simulation.number_of_bodies; i++)
			{
				if (omega_correction_on_body[i] == 1)
				{
					omega_correction_rot_after_simulation[i] 
						= 2.0 * M_PI / norm_vector(bodies[i].omega);

					if (simulation.omega_correction_counter == 1) // defines first step
					{
						omega_correction_step[i] = 
							bodies_copy[i].rot - omega_correction_rot_after_simulation[i];
					}
					else
					{
						double aux = bodies_copy[i].rot - omega_correction_rot_after_simulation[i];
						if (omega_correction_step[i] * aux < 0.0)
						{
							omega_correction_step[i] /= -2.1;
						}
					}
					omega_correction_rot[i] += omega_correction_step[i];

					copy_CelestialBody(&bodies[i], bodies_copy[i]);
					bodies[i].rot_ini = omega_correction_rot[i];
				}
				else
				{
					copy_CelestialBody(&bodies[i], bodies_copy[i]);
				}
			}
		}
	}

	/* initialize l and angular velocity */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].point_mass == true)
		{
			nan_matrix(bodies[i].b);
			nan_vector(bodies[i].l);
		}
		else
		{
			initialize_angular_velocity_on_z_axis(&bodies[i]);
			calculate_b(i, bodies, simulation.number_of_bodies, simulation.G);
			initialize_angular_velocity(&bodies[i]); // correct alignment
			calculate_l(&bodies[i]);
		}
	}

	/* field parameters */
	fldpar params = {simulation, bodies};

	/* state vector */
	size_t	dim_state_vec;
	double  *y;
	mount_state_vector(&y, &dim_state_vec, params);

	/* defining additional simulation parameters */
	if (simulation.t_step_received == false)
	{
		simulation.t_step = find_shortest_time_scale(bodies, simulation);
	}
	simulation.t_step_init = simulation.t_step / 5.0;
	simulation.t_step_min = simulation.t_step / 100.0;
	simulation.error_abs = 1.0e-13;
	simulation.error_rel = 0.0;

	/* set ODE numerical integrator (GSL) */
	gsl_odeiv2_system sys = {field, NULL, dim_state_vec, &params};
	gsl_odeiv2_driver *d = 
		gsl_odeiv2_driver_alloc_y_new(&sys, 
			gsl_odeiv2_step_rk8pd, simulation.t_step_init, 
			simulation.error_abs, simulation.error_rel);
	gsl_odeiv2_driver_set_hmin(d, simulation.t_step_min);

	/* create output files */
	FILE *out[simulation.number_of_bodies + 1];
	if (simulation.write_to_file == true)
	{
		create_output_files(bodies, simulation, out);
	}

	/* determine data skip based on output size */
	simulation.output_size = 100.0e6; // bytes
	calculate_data_skip(&simulation, bodies);

	/* integration loop */
	simulation.counter = 0;	
	simulation.t = simulation.t_init;
	double	final_time;
	int	loop_counter = 0;
	if (simulation.omega_correction == false)
	{
		final_time = simulation.t_final;
	}
	else
	{
		final_time = simulation.omega_correction_t_final;
	}
	// while (simulation.counter < 5) // for testing
	while (simulation.t < final_time)
	{
		if (fabs(final_time-simulation.t) < simulation.t_step)
		{
			simulation.t_step = final_time - simulation.t; // smaller last step
		}

		int status = gsl_odeiv2_driver_apply (d, &simulation.t, 
						simulation.t + simulation.t_step, y);

		if (status != GSL_SUCCESS)
		{
			fprintf(stderr, "Error \"%d\"", status);
			fprintf(stderr, " at time t = %1.5e\n", simulation.t);
			break;
		}

		/* update bodies */
		retrieve_state_vector (&bodies, y, simulation);

		/* update bs and ps */
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			if (bodies[i].point_mass == false)
			{
				if (bodies[i].deformable == false)
				{
					calculate_Y_and_Y_transpose(&bodies[i]);
					calculate_bs_me(&bodies[i]);
				}
				else if (bodies[i].prestress == true)
				{
					calculate_Y_and_Y_transpose(&bodies[i]);
					calculate_p_me(&bodies[i]);
				}
			}
		}

		/* update omega and b */
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{		
			if (bodies[i].point_mass == false)
			{
				calculate_omega(i, bodies, simulation.number_of_bodies, simulation.G);
				calculate_b(i, bodies, simulation.number_of_bodies, simulation.G);
			}
		}

		/* write output */
		if (simulation.write_to_file == true)
		{
			if (simulation.t > simulation.t_trans)
			{
				if (simulation.counter % simulation.data_skip == 0)
				{
					write_output(bodies, simulation, out);
				}
				simulation.counter++;
			}
		}

		loop_counter++;

	} // end while (simulation.t < final_time)

	/* close output files */
	if (simulation.write_to_file == true)
	{
		close_output_files(simulation, out);
	}

	/* free ODE driver (GSL) */
	gsl_odeiv2_driver_free(d);

	/* omega correction loop */
	if (simulation.omega_correction == true)
	{
		simulation.omega_correction_counter++;
		goto back_omega_correction;
	}

	/* free state vector */
	free(y);

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

	if (simulation.omega_correction == true)
	{
		free(omega_correction_on_body);
		free(omega_correction_step);
		free(omega_correction_rot);
		free(omega_correction_rot_after_simulation);
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			if (bodies_copy[i].elements > 0)
			{
				free(bodies_copy[i].alpha_elements);
				free(bodies_copy[i].eta_elements);	
				free(bodies_copy[i].bk_me);
			}
		}
		free(bodies_copy);
	}

	/* stop clock */
	end_time = clock();
	int time_spent_in_seconds 
		= (end_time - begin_time) / CLOCKS_PER_SEC;

	/* write overview file */
	if (simulation.write_to_file == true)
	{
		write_simulation_overview(time_spent_in_seconds, simulation);
	}

	return 0;
}
