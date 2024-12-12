#include "dynamical_system.h"

int
mount_state_vector (double **state_vector,
					size_t *dim_state_vector,
					const fldpar params)
{
	*state_vector = (double *) malloc (sizeof(double));

	siminf simulation = params.simulation;
	cltbdy *bodies = params.bodies;

	int	counter_mount_state_vec = 0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (simulation.number_of_bodies > 1)
		{
			for (int j = 0; j < 3; j++)
			{
				*state_vector = (double *) realloc (*state_vector, (counter_mount_state_vec + 1) * sizeof(double));
				(*state_vector)[counter_mount_state_vec++] = bodies[i].x[j];
			}
			for (int j = 0; j < 3; j++)
			{
				*state_vector = (double *) realloc (*state_vector, (counter_mount_state_vec + 1) * sizeof(double));
				(*state_vector)[counter_mount_state_vec++] = bodies[i].x_dot[j];
			}
		}
		if (bodies[i].point_mass == false)
		{
			for (int j = 0; j < 3; j++)
			{
				*state_vector = (double *) realloc (*state_vector, (counter_mount_state_vec + 1) * sizeof(double));
				(*state_vector)[counter_mount_state_vec++] = bodies[i].l[j];
			}
			for (int j = 0; j < 4; j++)
			{
				*state_vector = (double *) realloc (*state_vector, (counter_mount_state_vec + 1) * sizeof(double));
				(*state_vector)[counter_mount_state_vec++] = bodies[i].q[j];
			}
			if (bodies[i].deformable == true)
			{
				for (int j = 0; j < 5; j++)
				{
					*state_vector = (double *) realloc (*state_vector, (counter_mount_state_vec + 1) * sizeof(double));
					(*state_vector)[counter_mount_state_vec++] = bodies[i].b_eta_me[j];
				}
				for (int j = 0; j < 5 * bodies[i].elements; j++)
				{
					*state_vector = (double *) realloc (*state_vector, (counter_mount_state_vec + 1) * sizeof(double));
					(*state_vector)[counter_mount_state_vec++] = bodies[i].bk_me[j];
				}
			}
		}
	}
	*dim_state_vector = (size_t) counter_mount_state_vec;

	return 0;
}

int
retrieve_state_vector	(cltbdy **bodies,
						 const double *state_vector,
						 const siminf simulation)
{
	int	counter_retrieve_state_vec = 0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (simulation.number_of_bodies > 1)
		{
			for (int j = 0; j < 3; j++)
			{
				(*bodies)[i].x[j] 
					= state_vector[counter_retrieve_state_vec++];
			}
			for (int j = 0; j < 3; j++)
			{
				(*bodies)[i].x_dot[j] 
					= state_vector[counter_retrieve_state_vec++];
			}
		}
		if ((*bodies)[i].point_mass == false)
		{
			for (int j = 0; j < 3; j++)
			{
				(*bodies)[i].l[j] 
					= state_vector[counter_retrieve_state_vec++];
			}
			for (int j = 0; j < 4; j++)
			{
				(*bodies)[i].q[j] 
					= state_vector[counter_retrieve_state_vec++];
			}
			normalize_quaternion((*bodies)[i].q);
			if ((*bodies)[i].deformable == true)
			{
				for (int j = 0; j < 5; j++)
				{
					(*bodies)[i].b_eta_me[j] 
						= state_vector[counter_retrieve_state_vec++];
				}
				for (int j = 0; j < 5 * (*bodies)[i].elements; j++)
				{
					(*bodies)[i].bk_me[j] 
						= state_vector[counter_retrieve_state_vec++];
				}
			}
		}
	}

	// printf("%d\n", counter_retrieve_state_vec);

	return 0;
}

int
field(double t, 
	  const double y[],
	  double f[],
      void *params)
{
	(void)(t); // reinforce autonomous trait of the system

	/* preparing variables and parameters */
	fldpar *par = (fldpar *) params;
	siminf simulation = par->simulation;
	cltbdy *bodies = par->bodies;

	/* retrieve state vector */
	retrieve_state_vector (&bodies, y, simulation);

	/* calculate bs and ps for every body */
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

	/* calculate omega and b for every body */
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].point_mass == false)
		{
			calculate_omega(i, bodies, simulation.number_of_bodies, simulation.G);
			calculate_b(i, bodies, simulation.number_of_bodies, simulation.G);
		}
	}

	double 	component_x[simulation.number_of_bodies][3];
	double 	component_x_dot[simulation.number_of_bodies][3];
	double 	component_l[simulation.number_of_bodies][3];
	double 	component_b_eta_me[simulation.number_of_bodies][5];
	double 	***component_bk_me = *(&component_bk_me);
	int 	elements_total = 0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		elements_total += bodies[i].elements;
	}
	if (elements_total > 0)
	{
		component_bk_me	
			= (double ***) malloc(simulation.number_of_bodies * sizeof(double **));
		for (int i = 0; i < simulation.number_of_bodies; i++)
		{
			if (bodies[i].elements > 0)
			{
				component_bk_me[i] 
					= (double **) malloc(bodies[i].elements * sizeof(double *));

				for (int j = 0; j < bodies[i].elements; j++)
				{
					component_bk_me[i][j] = (double *) malloc(5 * sizeof(double));
				}
			}
		}
	}
	double component_q[simulation.number_of_bodies][4];

	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		/* calculating components */

		// x component

		if (simulation.number_of_bodies == 1)
		{
			nan_vector(component_x[i]);
		}
		else
		{
			copy_vector (component_x[i], bodies[i].x_dot);
		}

		// x_dot component

		if (simulation.number_of_bodies == 1)
		{
			nan_vector(component_x_dot[i]);
		}
		else
		{
			null_vector(component_x_dot[i]);
			if (simulation.keplerian_motion == true)
			{
				if (i > 0)
				{
					double relative_to_ref_x[3];
					linear_combination_vector(relative_to_ref_x,
						1.0, bodies[i].x,
						-1.0, bodies[0].x);

					double x_relative_to_ref_norm		
						= norm_vector(relative_to_ref_x);
					double x_relative_to_ref_norm_cube	
						= pow(x_relative_to_ref_norm, 3.0);
					
					double minus_G_times_total_mass 
						= -1.0 * simulation.G * (bodies[0].mass + bodies[i].mass);

					scale_vector (component_x_dot[i], 
						minus_G_times_total_mass / x_relative_to_ref_norm_cube,
						relative_to_ref_x);	
				}
			}
			else if (simulation.two_bodies_aprox == true)
			{
				if (i > 0)
				{
					double relative_to_ref_x[3];
					linear_combination_vector(relative_to_ref_x,
						1.0, bodies[i].x,
						-1.0, bodies[0].x);
					
					double x_relative_to_ref_norm			= norm_vector(relative_to_ref_x);
					double x_relative_to_ref_norm_cube		= pow(x_relative_to_ref_norm, 3.0);

					double minus_G_times_total_mass = -1.0 * simulation.G * (bodies[0].mass + bodies[i].mass);

					double component_x_dot_1st_term[] = { 0.0, 0.0, 0.0 };
					scale_vector (component_x_dot_1st_term, 
						1.0 / x_relative_to_ref_norm_cube, relative_to_ref_x);

					if (bodies[i].point_mass == true && bodies[0].point_mass == true)
					{
						scale_vector(component_x_dot[i], 
							minus_G_times_total_mass, component_x_dot_1st_term);
					}
					else
					{
						double scaled_sum_bi_bref[9];
						null_matrix(scaled_sum_bi_bref);
						if (bodies[i].point_mass == false && bodies[0].point_mass == false)
						{
							linear_combination_square_matrix(scaled_sum_bi_bref,
								bodies[i].I0 / bodies[i].mass, bodies[i].b,
								bodies[0].I0 / bodies[0].mass, bodies[0].b);
						}
						else if (bodies[i].point_mass == true)
						{
							scale_square_matrix(scaled_sum_bi_bref, 
								bodies[i].mass * bodies[0].I0, bodies[0].b);
						}
						else
						{
							scale_square_matrix(scaled_sum_bi_bref, 
								bodies[0].mass * bodies[i].I0, bodies[i].b);
						}

						double x_relative_to_ref_norm_seventh = pow(x_relative_to_ref_norm, 7.0);
						double scaled_sum_bi_bref_x[3];
						square_matrix_times_vector(scaled_sum_bi_bref_x, 
							scaled_sum_bi_bref, relative_to_ref_x);
						double scaled_sum_bi_bref_x_dot_x = 
							dot_product(scaled_sum_bi_bref_x, relative_to_ref_x);
						double component_x_dot_2nd_term[] = { 0.0, 0.0, 0.0 };
						scale_vector (component_x_dot_2nd_term, 
							(15.0 * scaled_sum_bi_bref_x_dot_x) / (2.0 * x_relative_to_ref_norm_seventh), 
							relative_to_ref_x);

						double x_relative_to_ref_norm_fifth	= pow(x_relative_to_ref_norm, 5.0);
						double component_x_dot_3rd_term[] = { 0.0, 0.0, 0.0 };
						scale_vector (component_x_dot_3rd_term, 
							-3.0 / x_relative_to_ref_norm_fifth, scaled_sum_bi_bref_x);

						linear_combination_three_vector(component_x_dot[i],
							minus_G_times_total_mass, component_x_dot_1st_term, 
							minus_G_times_total_mass, component_x_dot_2nd_term, 
							minus_G_times_total_mass, component_x_dot_3rd_term);
					}
				}
			}
			else
			{
				for (int j = 0; j < simulation.number_of_bodies; j++)
				{
					if (j != i)
					{
						double relative_x[3];
						linear_combination_vector(relative_x,
							1.0, bodies[i].x,
							-1.0, bodies[j].x);

						double x_relative_norm 		= norm_vector(relative_x);
						double x_relative_norm_cube = pow(x_relative_norm, 3.0);
						double component_x_dot_1st_term[] = { 0.0, 0.0, 0.0 };
						scale_vector (component_x_dot_1st_term, 
							-1.0 * bodies[j].mass / x_relative_norm_cube, relative_x);

						double j_component_x_dot[] = { 0.0, 0.0, 0.0 };
						if (bodies[i].point_mass == true && bodies[j].point_mass == true)
						{
							scale_vector(j_component_x_dot, simulation.G, component_x_dot_1st_term);
						}
						else
						{
							double scaled_sum_bi_bj[9];
							null_matrix(scaled_sum_bi_bj);
							if (bodies[i].point_mass == false && bodies[j].point_mass == false)
							{
								linear_combination_square_matrix(scaled_sum_bi_bj,
									bodies[j].mass * bodies[i].I0, bodies[i].b,
									bodies[i].mass * bodies[j].I0, bodies[j].b);
							}
							else if (bodies[i].point_mass == true)
							{
								scale_square_matrix(scaled_sum_bi_bj, 
									bodies[i].mass * bodies[j].I0, bodies[j].b);
							}
							else
							{
								scale_square_matrix(scaled_sum_bi_bj, 
									bodies[j].mass * bodies[i].I0, bodies[i].b);
							}

							double x_relative_norm_seventh = pow(x_relative_norm, 7.0);
							double scaled_sum_bi_bj_x[3];
							square_matrix_times_vector(scaled_sum_bi_bj_x, 
								scaled_sum_bi_bj, relative_x);
							double scaled_sum_bi_bj_x_dot_x = 
								dot_product(scaled_sum_bi_bj_x, relative_x);
							double component_x_dot_2nd_term[] = { 0.0, 0.0, 0.0 };
							scale_vector (component_x_dot_2nd_term, 
								(-15.0 * scaled_sum_bi_bj_x_dot_x) / (2. * bodies[i].mass * x_relative_norm_seventh), 
								relative_x);

							double x_relative_norm_fifth = pow(x_relative_norm, 5.0);
							double component_x_dot_3rd_term[] = { 0.0, 0.0, 0.0 };
							scale_vector (component_x_dot_3rd_term, 
								3.0 / (bodies[i].mass * x_relative_norm_fifth), scaled_sum_bi_bj_x);
							
							linear_combination_three_vector(j_component_x_dot,
								simulation.G, component_x_dot_1st_term, 
								simulation.G, component_x_dot_2nd_term, 
								simulation.G, component_x_dot_3rd_term);
						}

						linear_combination_vector(component_x_dot[i],
							1.0, component_x_dot[i],
							1.0, j_component_x_dot);

					} // end if (j != i)
				} // end summation on number of bodies
			} // end else for if (bodies[i].keplerian == false)
		}

		// l component

		if (bodies[i].point_mass == true)
		{
			nan_vector(component_l[i]);
		}
		else
		{
			null_vector(component_l[i]);
			for (int j = 0; j < simulation.number_of_bodies; j++)
			{
				if (j != i)
				{
					double relative_x[3];
					linear_combination_vector(relative_x,
						1.0, bodies[i].x,
						-1.0, bodies[j].x);

					double x_relative_norm 		 = norm_vector(relative_x);
					double x_relative_norm_fifth = pow(x_relative_norm, 5.0);

					double bx[3];
					square_matrix_times_vector(bx, bodies[i].b, relative_x);
					double x_cross_bx[3];
					cross_product(x_cross_bx, relative_x, bx);
					double j_component_l[] = { 0.0, 0.0, 0.0 };
					scale_vector (j_component_l, 
						(-3.0 * simulation.G * bodies[j].mass * bodies[i].I0) / x_relative_norm_fifth, x_cross_bx);

					linear_combination_vector(component_l[i],
						1.0, component_l[i],
						1.0, j_component_l);
				} 
			} 
		}

		// q component

		if (bodies[i].point_mass == true)
		{
			nan_quaternion(component_q[i]);
		}
		else
		{
			double half_omega[3];
			scale_vector(half_omega, 0.5, bodies[i].omega);
			double quaternion_half_omega[4];
			quaternion_from_vector(quaternion_half_omega, half_omega);
			quaternion_times_quaternion(component_q[i],
				quaternion_half_omega, bodies[i].q);
		}

		// b_eta_me and bk_me components

		if (bodies[i].point_mass == true || bodies[i].deformable == false)
		{
			double dummy_nan_matrix[9];
			nan_matrix(dummy_nan_matrix);			
			get_main_elements_traceless_symmetric_matrix(component_b_eta_me[i],
				dummy_nan_matrix);
		}
		else
		{
			double omega_hat[9];
			hat_map(omega_hat, bodies[i].omega);

			double b_eta[9];
			construct_traceless_symmetric_matrix(b_eta, bodies[i].b_eta_me);
			double lambda[9];
			linear_combination_square_matrix(lambda, 
				bodies[i].alpha, bodies[i].b, -1.0 * bodies[i].alpha, b_eta);

			if (bodies[i].elements > 0)
			{
				double	bk_me_2d_array[bodies[i].elements][5];
				double	bk[bodies[i].elements][9];
				for (int k = 0; k < bodies[i].elements; k++)
				{
					for (int l = 0; l < 5; l++)
					{
						bk_me_2d_array[k][l] = bodies[i].bk_me[l + (k*5)];

					}
					construct_traceless_symmetric_matrix(bk[k], bk_me_2d_array[k]);

					linear_combination_square_matrix(lambda, 
						1.0, lambda, -1.0 * bodies[i].alpha, bk[k]);
				}

				for (int k = 0; k < bodies[i].elements; k++)
				{
					double omega_hat_comm_bk[9];
					commutator(omega_hat_comm_bk, omega_hat, bk[k]);
					double minus_bk_over_tau_elements[9];
					scale_square_matrix(minus_bk_over_tau_elements, 
						-1.0 * bodies[i].alpha_elements[k] / bodies[i].eta_elements[k],
						bk[k]);
					double lambda_over_eta_elements[9];
					scale_square_matrix(lambda_over_eta_elements,
						1.0 / bodies[i].eta_elements[k], lambda);

					double component_bk[] = { 0.0, 0.0, 0.0,
											  0.0, 0.0, 0.0,
											  0.0, 0.0, 0.0 };
					linear_combination_three_square_matrix(component_bk,
						1.0, omega_hat_comm_bk,
						1.0, minus_bk_over_tau_elements,
						1.0, lambda_over_eta_elements);

					get_main_elements_traceless_symmetric_matrix(component_bk_me[i][k],
						component_bk);
				}
			} // end bodies[i].elements > 0

			double omega_hat_comm_b_eta[9];
			commutator(omega_hat_comm_b_eta, omega_hat, b_eta);
			double lambda_over_eta[9];
			scale_square_matrix(lambda_over_eta,
				1.0 / bodies[i].eta, lambda);
			double component_b_eta[] = { 0.0, 0.0, 0.0,
										 0.0, 0.0, 0.0,
										 0.0, 0.0, 0.0 };
			linear_combination_square_matrix(component_b_eta, 
				1.0, omega_hat_comm_b_eta, 1.0, lambda_over_eta);
			get_main_elements_traceless_symmetric_matrix(component_b_eta_me[i], component_b_eta);

		} // end else < if (bodies[i].point_mass == true || bodies[i].deformable == false) >

	} // end loop over bodies

	/* writing components */

	int counter_mount_field_vec = 0;
	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (simulation.number_of_bodies > 1)
		{
			for (int j = 0; j < 3; j++)
			{
				f[counter_mount_field_vec++] = component_x[i][j];
			}
			for (int j = 0; j < 3; j++)
			{
				f[counter_mount_field_vec++] = component_x_dot[i][j];
			}
		}
		if (bodies[i].point_mass == false)
		{
			for (int j = 0; j < 3; j++)
			{
				f[counter_mount_field_vec++] = component_l[i][j];
			}
			for (int j = 0; j < 4; j++)
			{
				f[counter_mount_field_vec++] = component_q[i][j];
			}
			if (bodies[i].deformable == true)
			{
				for (int j = 0; j < 5; j++)
				{
					f[counter_mount_field_vec++] = component_b_eta_me[i][j];
				}

				for (int k = 0; k < bodies[i].elements; k++)
				{
					for (int l = 0; l < 5; l++)
					{
						f[counter_mount_field_vec++] = component_bk_me[i][k][l];
					}
				}
			}
		}
	}

	for (int i = 0; i < simulation.number_of_bodies; i++)
	{
		if (bodies[i].elements > 0)
		{
			for (int j = 0; j < bodies[i].elements; j++)
			{
				free(component_bk_me[i][j]);
			}
			free(component_bk_me[i]);
		}
	}
	if (elements_total > 0)
	{
		free(component_bk_me);
	}
	
	return GSL_SUCCESS;
}
