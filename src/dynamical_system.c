#include "dynamical_system.h"

int
field_GV(double t, 
		 const double y[],
		 double f[],
       	 void *params)
{
	/* reinforce autonomous trait of the system */
	(void)(t);

	/* preparing variables and parameters */

	double 	*par = (double *)params;

	double 	G					= par[0];
	int		number_of_bodies	= (int) par[1];

	cltbdy 	*bodies;
	bodies = (cltbdy *) malloc (number_of_bodies * sizeof(cltbdy));

	int	dim_params_per_body_without_elements = 15;
	int	dim_state_per_body_without_elements = 23;
	int elements_total, elements_counter = 0; 
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		
		bodies[i].keplerian 			= (bool) par[2 + 0 + dim_params_skip];
		bodies[i].orbit_2body 			= (bool) par[2 + 1 + dim_params_skip];
		bodies[i].point_mass 			= (bool) par[2 + 2 + dim_params_skip];
		bodies[i].centrifugal 			= (bool) par[2 + 3 + dim_params_skip];
		bodies[i].tidal 				= (bool) par[2 + 4 + dim_params_skip];
		for (int j = 0; j < 3; j++)
		{
			bodies[i].omega[j] 			= par[2 + 5 + dim_params_skip + j];
		}
		bodies[i].mass 					= par[2 + 8 + dim_params_skip];
		bodies[i].I0 					= par[2 + 9 + dim_params_skip];
		bodies[i].gamma 				= par[2 + 10 + dim_params_skip];
		bodies[i].alpha 				= par[2 + 11 + dim_params_skip];
		bodies[i].eta					= par[2 + 12 + dim_params_skip];
		bodies[i].alpha_0 				= par[2 + 13 + dim_params_skip];
		bodies[i].elements				= (int) par[2 + 14 + dim_params_skip];
		if (bodies[i].elements > 0)
		{
			bodies[i].alpha_elements = (double *) malloc(bodies[i].elements * sizeof(double));
			for (int j = 0; j < bodies[i].elements; j++)
			{
				bodies[i].alpha_elements[j] 	= par[2 + 15 + dim_params_skip + j];
			}
			bodies[i].eta_elements = (double *) malloc(bodies[i].elements * sizeof(double));
			for (int j = 0; j < bodies[i].elements; j++)
			{
				bodies[i].eta_elements[j] 	= par[2 + 16 + dim_params_skip + j + bodies[i].elements - 1];
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
		// normalize_quaternion(bodies[i].q);

		elements_counter += bodies[i].elements;
	}
	elements_total = elements_counter;

	/* calculate omega and b for every body */
	for (int i = 0; i < number_of_bodies; i++)
	{
		// printf("omega of bodies %d before\n", i+1);
		// print_vector(bodies[i].omega);
		calculate_omega(i, bodies, number_of_bodies, G);
		// printf("omega of body %d after\n", i+1);
		// print_vector(bodies[i].omega);
		calculate_b(i, bodies, number_of_bodies, G);
	}

	/* for testing */
	// printf("omega inside = \n");
	// print_vector(omega);
	// printf("b = \n");
	// print_square_matrix(b);
	// exit(42);
	// null_matrix(b);
	// null_matrix(omega);
	// double b_me[5];
	// for (int i = 0; i < 5; i++)
	// {
	// 	b_me[i] = ((double) i) * 0.00000000001;
	// }
	// construct_traceless_symmetric_matrix(b, b_me);
	// printf("Y in field_GV = \n");
	// print_square_matrix(bodies[0].Y);

	// double **component_x;
	// double **component_x_dot;
	// double **component_l;
	// double **component_b0_me;
	// double **component_u_me;
	// double ***component_bk_me = *(&component_bk_me);

	// component_x	= (double **) malloc(number_of_bodies * sizeof(double *));
	// component_x_dot	= (double **) malloc(number_of_bodies * sizeof(double *));
	// component_l	= (double **) malloc(number_of_bodies * sizeof(double *));
	// component_b0_me	= (double **) malloc(number_of_bodies * sizeof(double *));
	// component_u_me	= (double **) malloc(number_of_bodies * sizeof(double *));
	// if (elements_total > 0)
	// {
	// 	component_bk_me	= (double ***) malloc(number_of_bodies * sizeof(double **));
	// }
	// for (int i = 0; i < number_of_bodies; i++)
	// {
	// 	component_x[i] = (double *) malloc(3 * sizeof(double));
	// 	component_x_dot[i] = (double *) malloc(3 * sizeof(double));
	// 	component_l[i] = (double *) malloc(3 * sizeof(double));
	// 	component_b0_me[i] = (double *) malloc(5 * sizeof(double));
	// 	component_u_me[i] = (double *) malloc(5 * sizeof(double));
	// 	if (bodies[i].elements > 0)
	// 	{
	// 		component_bk_me[i] = (double **) malloc(bodies[i].elements * sizeof(double *));

	// 		for (int j = 0; j < bodies[i].elements; j++)
	// 		{
	// 			component_bk_me[i][j] = (double *) malloc(5 * sizeof(double));
	// 		}
	// 	}
	// }

	double component_x[number_of_bodies][3];
	double component_x_dot[number_of_bodies][3];
	double component_l[number_of_bodies][3];
	double component_b0_me[number_of_bodies][5];
	double component_u_me[number_of_bodies][5];
	double ***component_bk_me = *(&component_bk_me);
	if (elements_total > 0)
	{
		component_bk_me	= (double ***) malloc(number_of_bodies * sizeof(double **));
	}
	for (int i = 0; i < number_of_bodies; i++)
	{
		if (bodies[i].elements > 0)
		{
			component_bk_me[i] = (double **) malloc(bodies[i].elements * sizeof(double *));

			for (int j = 0; j < bodies[i].elements; j++)
			{
				component_bk_me[i][j] = (double *) malloc(5 * sizeof(double));
			}
		}
	}
	double component_q[number_of_bodies][4];

	for (int i = 0; i < number_of_bodies; i++)
	{
		double omega_hat[9];
		hat_map(omega_hat, bodies[i].omega);

		/* calculating components */

		// x component

		copy_vector (component_x[i], bodies[i].x_dot);

		// x_dot component

		null_vector(component_x_dot[i]);
		if (bodies[i].keplerian == true)
		{
			if (i > 0)
			{
				double relative_to_ref_x[3];
				linear_combination_vector(relative_to_ref_x,
					1.0, bodies[i].x,
					-1.0, bodies[0].x);

				double x_relative_to_ref_norm		= norm_vector(relative_to_ref_x);
				double x_relative_to_ref_norm_cube	= pow(x_relative_to_ref_norm, 3.0);
				
				double minus_G_times_total_mass = -1.0 * G * (bodies[0].mass + bodies[i].mass);

				scale_vector (component_x_dot[i], 
					minus_G_times_total_mass / x_relative_to_ref_norm_cube, relative_to_ref_x);	
			}
		}
		else if (bodies[i].orbit_2body == true)
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
			} // end if (j != i)
		} // end summation on number of bodies

		// b0 component

		double b0[9];
		construct_traceless_symmetric_matrix(b0, bodies[i].b0_me);
		double component_b0[] = { 0.0, 0.0, 0.0,
								  0.0, 0.0, 0.0,
								  0.0, 0.0, 0.0 };
		commutator(component_b0, omega_hat, b0);
		get_main_elements_traceless_symmetric_matrix(component_b0_me[i], component_b0);

		// u and bk components

		double u[9];
		construct_traceless_symmetric_matrix(u, bodies[i].u_me);
		double lambda[9];
		linear_combination_square_matrix(lambda, 
			1.0, u, bodies[i].alpha, bodies[i].b);

		if (bodies[i].elements > 0)
		{
			double	bk_me_2d_array[bodies[i].elements][5];
			double	bk[bodies[i].elements][9];
			for (int k = 0; k < bodies[i].elements; k++)
			{
				for (int l = 0; l < 5; l++)
				{
					bk_me_2d_array[k][l] = bodies[i].bk_me[l + (k*5)];

					/* for testing */
					// printf("bk_me = %f\n",bk_me[i][j]);
				}
				construct_traceless_symmetric_matrix(bk[k], bk_me_2d_array[k]);

				linear_combination_square_matrix(lambda, 
					1.0, lambda, -1.0 * bodies[i].alpha, bk[k]);
			}

			for (int k = 0; k < bodies[i].elements; k++)
			{
				double omega_hat_comm_bk[9];
				commutator(omega_hat_comm_bk, omega_hat, bk[k]);
				double tau_elements = 
					bodies[i].eta_elements[k] / bodies[i].alpha_elements[k];
				double minus_bk_over_tau_elements[9];
				scale_square_matrix(minus_bk_over_tau_elements, 
					-1.0 / tau_elements, bk[k]);
				double lambda_over_eta_elements[9];
				scale_square_matrix(lambda_over_eta_elements,
					1.0 / bodies[i].eta_elements[k], lambda);

				/* for testing */
				// printf("\ntau_%d = %f\n", i, tau_elements[i]);
				// printf("\nminus_1_over_tau_%d = %f\n", i, -1.0 / tau_elements[i]);
				// printf("\nomega_hat_comm_bk = \n");
				// print_square_matrix(omega_hat_comm_bk);
				// printf("\nminus_bk_over_tau_elements = \n");
				// print_square_matrix(minus_bk_over_tau_elements);
				// printf("\nlambda_over_eta_elements = \n");
				// print_square_matrix(lambda_over_eta_elements);

				double component_bk[] = { 0.0, 0.0, 0.0,
										0.0, 0.0, 0.0,
										0.0, 0.0, 0.0 };
				linear_combination_three_square_matrix(component_bk,
					1.0, omega_hat_comm_bk,
					1.0, minus_bk_over_tau_elements,
					1.0, lambda_over_eta_elements);

				get_main_elements_traceless_symmetric_matrix(component_bk_me[i][k],
					component_bk);

				/* for testing */
				// printf("\nalpha_%d = %f\n", i, alpha_elements[i]);
				// printf("\ntau_%d = %f\n", i, tau_elements[i]);
				// printf("\ncomponent_bk = \n");
				// print_square_matrix(component_bk);
			}
		} // end bodies[i].elements > 0

		double omega_hat_comm_u[9];
		commutator(omega_hat_comm_u, omega_hat, u);
		double tau = bodies[i].eta / bodies[i].alpha;
		double component_u[] = { 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0 };
		linear_combination_square_matrix(component_u, 
			1.0, omega_hat_comm_u, -1.0 / tau, lambda);
		get_main_elements_traceless_symmetric_matrix(component_u_me[i], component_u);

		// q component
		double half_omega[3];
		scale_vector(half_omega, 0.5, bodies[i].omega);
		double quaternion_half_omega[4];
		quaternion_from_vector(quaternion_half_omega, half_omega);
		quaternion_times_quaternion(component_q[i],
			quaternion_half_omega, bodies[i].q);

	} // end loop over bodies

	/* for testing */
	// printf("\nomega_hat_comm_u = \n");
	// print_square_matrix(omega_hat_comm_u);
	// printf("\ntau = \n");
	// printf("%f\n", tau);
	// printf("\nlambda = \n");
	// print_square_matrix(lambda);
	// printf("\ncomponent_u = \n");
	// print_square_matrix(component_u);	
	// printf("\ncomponent_u_me = \n");
	// printf("%1.10e %1.10e %1.10e %1.10e %1.10e\n",
	// 	component_u_me[0], component_u_me[1], component_u_me[2],
	// 	component_u_me[3], component_u_me[4]);	
	// printf("\nu_me = \n");
	// printf("%1.10e %1.10e %1.10e %1.10e %1.10e\n",
	// 	u_me[0], u_me[1], u_me[2], u_me[3], u_me[4]);
	// printf("component_Y = \n");
	// print_square_matrix(component_Y[0]);
	// exit(42);

	/* writing components */

	elements_counter = 0;
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_state_skip = i * dim_state_per_body_without_elements + 5 * elements_counter;
		for (int j = 0; j < 3; j++)
		{
			f[0 + dim_state_skip + j] = component_x[i][j];
			f[3 + dim_state_skip + j] = component_x_dot[i][j];
			f[6 + dim_state_skip + j] = component_l[i][j];
		}
		for (int j = 0; j < 5; j++)
		{
			f[9 + dim_state_skip + j] = component_b0_me[i][j];
			f[14 + dim_state_skip + j] = component_u_me[i][j];
		}
		for (int k = 0; k < bodies[i].elements; k++)
		{
			for (int l = 0; l < 5; l++)
			{
				f[19 + dim_state_skip + 5*k + l] = component_bk_me[i][k][l];
			}
		}
		for (int j = 0; j < 4; j++)
		{
			f[19 + 5 * bodies[i].elements + dim_state_skip + j]	= component_q[i][j];
		}
		elements_counter += bodies[i].elements;

	} // end for (int i = 0; i < number_of_bodies; i++)

	/* for testing */
	// for (int i = 0; i < number_of_bodies * 19; i++)
	// {
	// 	printf("f[%d] = %1.5e\n", i, f[i]);
	// }

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

	/* for testing */
	// for (int i = 0; i < number_of_bodies; i++)
	// 		print_CelestialBody(bodies[i]);
	// exit(99);
	
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

	// printf("here again\n");

	return GSL_SUCCESS;
}
