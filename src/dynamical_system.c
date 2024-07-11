#include "dynamical_system.h"

int
field_gV(double t, 
		 const double y[],
		 double f[],
       	 void *params)
{
	(void)(t); // reinforce autonomous trait of the system

	/* preparing variables and parameters */

	double 	*par = (double *)params;

	double 	G					= par[0];
	int		number_of_bodies	= (int) par[1];
	bool	keplerian_motion	= (bool) par[2];
	bool	two_bodies_aprox	= (bool) par[3];

	cltbdy 	*bodies;
	bodies = (cltbdy *) malloc (number_of_bodies * sizeof(cltbdy));

	int	dim_params_per_body_without_elements = 21;
	// int	dim_state_per_body_without_elements = 18;
	int	dim_state_per_body_without_elements = 33;
	int elements_total, elements_counter = 0; 
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_params_skip 
			= i * dim_params_per_body_without_elements + 2 * elements_counter;
		
		bodies[i].point_mass  = (bool) par[4 + 0 + dim_params_skip];
		bodies[i].prestress   = (bool) par[4 + 1 + dim_params_skip];
		bodies[i].centrifugal = (bool) par[4 + 2 + dim_params_skip];
		bodies[i].tidal 	  = (bool) par[4 + 3 + dim_params_skip];
		bodies[i].deformable  = (bool) par[4 + 4 + dim_params_skip];
		bodies[i].mass 		  = par[4 + 5 + dim_params_skip];
		bodies[i].I0 		  = par[4 + 6 + dim_params_skip];
		bodies[i].gamma_0 	  = par[4 + 7 + dim_params_skip];
		bodies[i].alpha 	  = par[4 + 8 + dim_params_skip];
		bodies[i].eta		  = par[4 + 9 + dim_params_skip];
		for (int j = 0; j < 5; j++)
		{
			bodies[i].Bs_me[j] = par[4 + 10 + j + dim_params_skip];
		}
		for (int j = 0; j < 5; j++)
		{
			bodies[i].P_me[j] = par[4 + 15 + j + dim_params_skip];
		}
		bodies[i].elements	  = (int) par[4 + 20 + dim_params_skip];
		if (bodies[i].elements > 0)
		{
			bodies[i].alpha_elements 
				= (double *) malloc(bodies[i].elements * sizeof(double));
			for (int j = 0; j < bodies[i].elements; j++)
			{
				bodies[i].alpha_elements[j] 
					= par[4 + 21 + dim_params_skip + j];
			}
			bodies[i].eta_elements 
				= (double *) malloc(bodies[i].elements * sizeof(double));
			for (int j = 0; j < bodies[i].elements; j++)
			{
				bodies[i].eta_elements[j] 
					= par[4 + 22 + dim_params_skip + j + bodies[i].elements - 1];
			}
		}

		int	dim_state_skip = i * dim_state_per_body_without_elements + 5 * elements_counter;
		
		for (int j = 0; j < 3; j++)
		{
			bodies[i].x[j] 		= y[0 + dim_state_skip + j];
			bodies[i].x_dot[j] 	= y[3 + dim_state_skip + j];
			bodies[i].l[j] 		= y[6 + dim_state_skip + j];
		}
		for (int j = 0; j < 5; j++)
		{
			bodies[i].b_eta_me[j] = y[9 + dim_state_skip + j];
		}
		if (bodies[i].elements > 0)
		{
			bodies[i].bk_me 
				= (double *) malloc(5 * bodies[i].elements * sizeof(double));
			for (int j = 0; j < 5 * bodies[i].elements; j++)
			{
				bodies[i].bk_me[j] = y[14 + dim_state_skip + j];
			}
		}
		// for (int j = 0; j < 4; j++)
		// {
		// 	bodies[i].q[j] 
		// 		= y[14 + 5 * bodies[i].elements + dim_state_skip + j];
		// }
		for (int j = 0; j < 9; j++)
		{
			bodies[i].Y[j] = y[14 + 5 * bodies[i].elements + dim_state_skip + j];
		}
		transpose_square_matrix(bodies[i].Y_trans, bodies[i].Y);
		for (int j = 0; j < 5; j++)
		{
			bodies[i].bs_me[j] = y[23 + 5 * bodies[i].elements + dim_state_skip + j];
		}
		for (int j = 0; j < 5; j++)
		{
			bodies[i].p_me[j] = y[28 + 5 * bodies[i].elements + dim_state_skip + j];
		}
		elements_counter += bodies[i].elements;
	}
	elements_total = elements_counter;

	// /* calculate bs and ps for every body */
	// for (int i = 0; i < number_of_bodies; i++)
	// {
	// 	calculate_Y_and_Y_transpose_via_quaternion(&bodies[i]);
	// 	calculate_bs_me(&bodies[i]);
	// 	calculate_p_me(&bodies[i]);
	// }

	/* calculate omega and b for every body */
	for (int i = 0; i < number_of_bodies; i++)
	{
		calculate_omega(i, bodies, number_of_bodies, G);
		calculate_b(i, bodies, number_of_bodies, G);
		// print_CelestialBody(bodies[i]);
	}
	// exit(87);

	double component_x[number_of_bodies][3];
	double component_x_dot[number_of_bodies][3];
	double component_l[number_of_bodies][3];
	double component_b_eta_me[number_of_bodies][5];
	double ***component_bk_me = *(&component_bk_me);
	if (elements_total > 0)
	{
		component_bk_me	
			= (double ***) malloc(number_of_bodies * sizeof(double **));
	}
	for (int i = 0; i < number_of_bodies; i++)
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
	// double component_q[number_of_bodies][4];
	double component_Y[number_of_bodies][9];
	double component_bs_me[number_of_bodies][5];
	double component_p_me[number_of_bodies][5];

	for (int i = 0; i < number_of_bodies; i++)
	{
		double omega_hat[9];
		hat_map(omega_hat, bodies[i].omega);

		/* calculating components */

		// x component

		copy_vector (component_x[i], bodies[i].x_dot);

		// x_dot component

		null_vector(component_x_dot[i]);
		if (keplerian_motion == true)
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
					= -1.0 * G * (bodies[0].mass + bodies[i].mass);

				scale_vector (component_x_dot[i], 
					minus_G_times_total_mass / x_relative_to_ref_norm_cube,
					relative_to_ref_x);	
			}
		}
		else if (two_bodies_aprox == true)
		{
			if (i > 0)
			{
				double relative_to_ref_x[3];
				linear_combination_vector(relative_to_ref_x,
					1.0, bodies[i].x,
					-1.0, bodies[0].x);
				
				double x_relative_to_ref_norm			= norm_vector(relative_to_ref_x);
				double x_relative_to_ref_norm_cube		= pow(x_relative_to_ref_norm, 3.0);
				double x_relative_to_ref_norm_fifth		= pow(x_relative_to_ref_norm, 5.0);
				double x_relative_to_ref_norm_seventh	= pow(x_relative_to_ref_norm, 7.0);

				double minus_G_times_total_mass = -1.0 * G * (bodies[0].mass + bodies[i].mass);

				double component_x_dot_1st_term[] = { 0.0, 0.0, 0.0 };
				scale_vector (component_x_dot_1st_term, 
					1.0 / x_relative_to_ref_norm_cube, relative_to_ref_x);

				double component_x_dot_2nd_term[] = { 0.0, 0.0, 0.0 };
				double scaled_sum_bi_bref[9];
				linear_combination_square_matrix(scaled_sum_bi_bref,
					bodies[i].I0 / bodies[i].mass, bodies[i].b,
					bodies[0].I0 / bodies[0].mass, bodies[0].b);
				double scaled_sum_bi_bref_x[3];
				square_matrix_times_vector(scaled_sum_bi_bref_x, 
					scaled_sum_bi_bref, relative_to_ref_x);
				double scaled_sum_bi_bref_x_dot_x = 
					dot_product(scaled_sum_bi_bref_x, relative_to_ref_x);
				scale_vector (component_x_dot_2nd_term, 
					(15.0 * scaled_sum_bi_bref_x_dot_x) / (2.0 * x_relative_to_ref_norm_seventh), 
					relative_to_ref_x);

				double component_x_dot_3rd_term[] = { 0.0, 0.0, 0.0 };
				scale_vector (component_x_dot_3rd_term, 
					-3.0 / x_relative_to_ref_norm_fifth, scaled_sum_bi_bref_x);

				linear_combination_three_vector(component_x_dot[i],
					minus_G_times_total_mass, component_x_dot_1st_term, 
					minus_G_times_total_mass, component_x_dot_2nd_term, 
					minus_G_times_total_mass, component_x_dot_3rd_term);
			}
		}
		else
		{
			for (int j = 0; j < number_of_bodies; j++)
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

					double x_relative_norm_seventh = pow(x_relative_norm, 7.0);
					double scaled_sum_bi_bj[9];
					linear_combination_square_matrix(scaled_sum_bi_bj,
						bodies[j].mass * bodies[i].I0, bodies[i].b,
						bodies[i].mass * bodies[j].I0, bodies[j].b);
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

					double j_component_x_dot[] = { 0.0, 0.0, 0.0 };
					linear_combination_three_vector(j_component_x_dot,
						G, component_x_dot_1st_term, 
						G, component_x_dot_2nd_term, 
						G, component_x_dot_3rd_term);

					linear_combination_vector(component_x_dot[i],
						1.0, component_x_dot[i],
						1.0, j_component_x_dot);

				} // end if (j != i)
			} // end summation on number of bodies
		} // end else for if (bodies[i].keplerian == false)

		// l component

		null_vector(component_l[i]);
		for (int j = 0; j < number_of_bodies; j++)
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
					(-3.0 * G * bodies[j].mass * bodies[i].I0) / x_relative_norm_fifth, x_cross_bx);

				linear_combination_vector(component_l[i],
					1.0, component_l[i],
					1.0, j_component_l);
			} 
		} 

		// b_eta and bk components

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

		// // q component
		// double half_omega[3];
		// scale_vector(half_omega, 0.5, bodies[i].omega);
		// double quaternion_half_omega[4];
		// quaternion_from_vector(quaternion_half_omega, half_omega);
		// quaternion_times_quaternion(component_q[i],
		// 	quaternion_half_omega, bodies[i].q);

		// Y component
		square_matrix_times_square_matrix(component_Y[i], 
			omega_hat, bodies[i].Y);

		// bs component
		double bs[9];
		construct_traceless_symmetric_matrix(bs, bodies[i].bs_me);
		double component_bs[] = { 0.0, 0.0, 0.0,
								  0.0, 0.0, 0.0,
								  0.0, 0.0, 0.0 };
		commutator(component_bs, omega_hat, bs);
		get_main_elements_traceless_symmetric_matrix(component_bs_me[i], 
			component_bs);

		// p component
		double p[9];
		construct_traceless_symmetric_matrix(p, bodies[i].p_me);
		double component_p[] = { 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0 };
		commutator(component_p, omega_hat, p);
		get_main_elements_traceless_symmetric_matrix(component_p_me[i], 
			component_p);

	} // end loop over bodies

	/* writing components */

	elements_counter = 0;
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_state_skip 
			= i * dim_state_per_body_without_elements 
			+ 5 * elements_counter;
		for (int j = 0; j < 3; j++)
		{
			f[0 + dim_state_skip + j] = component_x[i][j];
			f[3 + dim_state_skip + j] = component_x_dot[i][j];
			f[6 + dim_state_skip + j] = component_l[i][j];
		}
		for (int j = 0; j < 5; j++)
		{
			f[9 + dim_state_skip + j] = component_b_eta_me[i][j];
		}
		for (int k = 0; k < bodies[i].elements; k++)
		{
			for (int l = 0; l < 5; l++)
			{
				f[14 + dim_state_skip + 5*k + l]
					= component_bk_me[i][k][l];
			}
		}
		// for (int j = 0; j < 4; j++)
		// {
		// 	f[14 + 5 * bodies[i].elements + dim_state_skip + j]	
		// 		= component_q[i][j];
		// }
		for (int j = 0; j < 9; j++)
		{
			f[14 + 5 * bodies[i].elements + dim_state_skip + j]	
				= component_Y[i][j];
		}
		for (int j = 0; j < 5; j++)
		{
			f[23 + 5 * bodies[i].elements + dim_state_skip + j]	
				= component_bs_me[i][j];
		}
		for (int j = 0; j < 5; j++)
		{
			f[28 + 5 * bodies[i].elements + dim_state_skip + j]	
				= component_p_me[i][j];
		}
		elements_counter += bodies[i].elements;
	}

	for (int i = 0; i < number_of_bodies; i++)
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
	
	/* free Voigt elements */
	for (int i = 0; i < number_of_bodies; i++)
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

	return GSL_SUCCESS;
}
