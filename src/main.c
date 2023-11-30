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

	/* correction for the angular velocity's initial value */
	int		omega_correction_number_of_iterates = 20;
	bool	omega_correction_on_body[simulation.number_of_bodies];
	double	omega_correction_step[simulation.number_of_bodies];
	double 	omega_correction_lod[simulation.number_of_bodies];
	double	omega_correction_lod_after_simulation[simulation.number_of_bodies];
	siminf	simulation_copy = *(&simulation_copy);
	cltbdy	*bodies_copy = *(&bodies_copy);
	if (simulation.omega_correction == true)
	{
		bodies_copy = (cltbdy *) malloc(simulation.number_of_bodies * sizeof(cltbdy));
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			if (bodies[i].elements > 0)
			{
				bodies_copy[i].alpha_elements = (double *) malloc(bodies[i].elements * sizeof(double));
				bodies_copy[i].eta_elements = (double *) malloc(bodies[i].elements * sizeof(double));
			}
			bodies_copy[i] = bodies[i];
			if (bodies[i].point_mass == false)
			{
				omega_correction_on_body[i] = true;
			}
			else
			{
				omega_correction_on_body[i] = false;
			}
			omega_correction_lod[i] = bodies[i].lod;
		}
		simulation.omega_correction_t_final = 
			10.0 * largest_time_scale(bodies, 
						simulation.number_of_bodies,
						simulation.G);
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
				if (omega_correction_on_body[i] == true)
				{
					omega_correction_lod_after_simulation[i] 
						= 2.0 * M_PI / norm_vector(bodies[i].omega);

					if (simulation.omega_correction_counter == 1) // defines first step
					{
						omega_correction_step[i] = 
							bodies_copy[i].lod - omega_correction_lod_after_simulation[i];
					}
					else
					{
						double aux = bodies_copy[i].lod - omega_correction_lod_after_simulation[i];
						if (omega_correction_step[i] * aux < 0.0)
						{
							omega_correction_step[i] /= -2.1;
						}
					}
					omega_correction_lod[i] += omega_correction_step[i];

					bodies[i] = bodies_copy[i];
					bodies[i].lod = omega_correction_lod[i];

					initialize_angular_velocity_on_figure_axis_of_solid_frame(&bodies[i]);

				}
				else
				{
					bodies[i] = bodies_copy[i]; // reboot to initial state
				}
			}
		}
	}

	/* calculate bk, b0, and u */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		/* bk */
		if (bodies[i].elements > 0)
		{
			bodies[i].bk_me = (double *) calloc(bodies[i].elements * 5, sizeof(double));
		}

		/* transformation matrices */
		double Y_i[9], Y_i_trans[9];
		rotation_matrix_from_quaternion(Y_i, bodies[i].q);
		transpose_square_matrix(Y_i_trans, Y_i);

		/* deformation components - I/II */
		double f_cent_i[9];
		null_matrix(f_cent_i);
		if (bodies[i].centrifugal == true)
		{
			calculate_f_cent(f_cent_i, bodies[i].omega);
		}
		double f_tide_i[9];
		null_matrix(f_tide_i);
		if (bodies[i].tidal == true)
		{
			calculate_f_tide(f_tide_i, i, bodies, 
				simulation.number_of_bodies, simulation.G);
		}

		/* deformation from stokes coefficients */
		double B_stokes_i[9];
		body_frame_deformation_from_stokes_coefficients(B_stokes_i, bodies[i]);
		double b_stokes_i[9];
		square_matrix_times_square_matrix(b_stokes_i, Y_i, B_stokes_i);
		square_matrix_times_square_matrix(b_stokes_i, b_stokes_i, Y_i_trans);

		/* prestress */
		double b0_i[9];	
		null_matrix(b0_i);
		if(bodies[i].prestress == true)
		{
			if (bodies[i].deformable == false)
			{
				copy_square_matrix(b0_i, b_stokes_i);
			}
			else
			{
				linear_combination_three_square_matrix(b0_i,
					(bodies[i].gamma+bodies[i].alpha_0)/bodies[i].alpha_0, b_stokes_i,
					-1.0/bodies[i].alpha_0, f_cent_i,
					-1.0/bodies[i].alpha_0, f_tide_i);
			}
		}
		get_main_elements_traceless_symmetric_matrix(bodies[i].b0_me, b0_i);

		/* deformation components - II/II */
		double f_ps_i[9];
		null_matrix(f_ps_i);
		if (bodies[i].prestress == true)
		{
			calculate_f_ps(f_ps_i, bodies[i]);
		}
		double c_i = calculate_c(bodies[i]);

		/* rheology variables */
		double u_i[9];
		null_matrix(u_i);
		if (bodies[i].deformable == true)
		{
			linear_combination_four_square_matrix(u_i,
				-1.0 * c_i, b_stokes_i,
					   1.0, f_cent_i,
					   1.0, f_tide_i,
					   1.0, f_ps_i);
		}
		get_main_elements_traceless_symmetric_matrix(bodies[i].u_me, u_i);
		
	} // end loop over bodies

	/* calculate first l */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		calculate_b(i, bodies, simulation.number_of_bodies, simulation.G);
		initialize_angular_velocity(&bodies[i]); // correct alignment
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
	int		dim_params_per_body_without_elements = 12;
	int		dim_params = 4 + (dim_params_per_body_without_elements * simulation.number_of_bodies) + (2 * elements_total);
	double	params[dim_params];
	elements_counter = 0; 
	params[0] = simulation.G;
	params[1] = (double) simulation.number_of_bodies;
	params[2] = (double) simulation.keplerian_motion;
	params[3] = (double) simulation.two_bodies_aprox;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		params[4 + 0 + dim_params_skip] = (double) bodies[i].point_mass;
		params[4 + 1 + dim_params_skip] = (double) bodies[i].prestress;
		params[4 + 2 + dim_params_skip] = (double) bodies[i].centrifugal;
		params[4 + 3 + dim_params_skip] = (double) bodies[i].tidal;
		params[4 + 4 + dim_params_skip] = (double) bodies[i].deformable;
		params[4 + 5 + dim_params_skip] = bodies[i].mass;
		params[4 + 6 + dim_params_skip] = bodies[i].I0;
		params[4 + 7 + dim_params_skip] = bodies[i].gamma;
		params[4 + 8 + dim_params_skip] = bodies[i].alpha;
		params[4 + 9 + dim_params_skip] = bodies[i].eta;
		params[4 + 10 + dim_params_skip] = bodies[i].alpha_0;
		params[4 + 11 + dim_params_skip] = (double) bodies[i].elements;
		for (int j = 0; j < bodies[i].elements; j++)
		{
			params[4 + 12 + dim_params_skip + j] = bodies[i].alpha_elements[j];
		}
		for (int j = 0; j < bodies[i].elements; j++)
		{
			params[4 + 13 + dim_params_skip + j + bodies[i].elements - 1] = bodies[i].eta_elements[j];
		}
		elements_counter += bodies[i].elements;		
	}

	/* GSL variables */
	const gsl_odeiv2_step_type * ode_type
		= gsl_odeiv2_step_rk8pd; // gsl_odeiv2_step_rk8pd
	gsl_odeiv2_step * ode_step
    	= gsl_odeiv2_step_alloc (ode_type, dim_state);
  	gsl_odeiv2_control * ode_control
    	= gsl_odeiv2_control_y_new (simulation.eps_abs, simulation.eps_rel);
  	gsl_odeiv2_evolve * ode_evolve
    	= gsl_odeiv2_evolve_alloc (dim_state);

	/* create output files */
	FILE *out[simulation.number_of_bodies + 1];
	if (simulation.write_to_file == true)
	{
		create_output_files(bodies, simulation, out);
	}

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
	while (simulation.t < final_time)
	{
		if (fabs(final_time-simulation.t) < simulation.t_step)
		{
			simulation.t_step = final_time - simulation.t; // smaller last step
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

		/* update omega and b */
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{		
			calculate_omega(i, bodies, simulation.number_of_bodies, simulation.G);
			calculate_b(i, bodies, simulation.number_of_bodies, simulation.G);
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

	/* free GSL variables */
	gsl_odeiv2_evolve_free (ode_evolve);
	gsl_odeiv2_control_free (ode_control);
	gsl_odeiv2_step_free (ode_step);

	if (simulation.omega_correction == true)
	{
		simulation.omega_correction_counter++;
		goto back_omega_correction;
	}

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
