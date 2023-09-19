#include "dynamical_system.h"

int
print_CelestialBody(cltbdy body)
{
	printf("name = ");
	printf("%s\n", body.name);

	printf("mass = ");
	printf("%1.5e\n", body.mass);
	printf("R = ");	
	printf("%1.5e\n", body.R);
	printf("lod = ");	
	printf("%1.5e\n", body.lod);

	printf("rg = ");	
	printf("%1.5e\n", body.rg);
	printf("J2 = ");	
	printf("%1.5e\n", body.J2);
	printf("C22 = ");	
	printf("%1.5e\n", body.C22);
	printf("I0 = ");	
	printf("%1.5e\n", body.I0);

	printf("obl = ");
	printf("%1.5e\n", body.obl);
	printf("psi = ");	
	printf("%1.5e\n", body.psi);
	printf("lib = ");	
	printf("%1.5e\n", body.lib);

	printf("a = ");
	printf("%1.5e\n", body.a);
	printf("e = ");	
	printf("%1.5e\n", body.e);
	printf("I = ");	
	printf("%1.5e\n", body.I);
	printf("M = ");	
	printf("%1.5e\n", body.M);
	printf("w = ");	
	printf("%1.5e\n", body.w);
	printf("Omega = ");	
	printf("%1.5e\n", body.Omega);

	printf("kf = ");	
	printf("%1.5e\n", body.kf);
	printf("Dt = ");	
	printf("%1.5e\n", body.Dt);
	printf("tau = ");	
	printf("%1.5e\n", body.tau);

	printf("gamma = ");	
	printf("%1.5e\n", body.gamma);
	printf("alpha = ");	
	printf("%1.5e\n", body.alpha);
	printf("eta = ");	
	printf("%1.5e\n", body.eta);
	printf("alpha_0 = ");	
	printf("%1.5e\n", body.alpha_0);
	printf("elements = ");	
	printf("%d\n", body.elements);
	for(int i = 0; i < body.elements; i++)
	{
		printf("alpha_%d = ", i+1);	
		printf("%1.5e\n", body.alpha_elements[i]);
		printf("eta_%d = ", i+1);	
		printf("%1.5e\n", body.eta_elements[i]);
	}

	printf("point mass = ");	
	printf("%d\n", body.point_mass);
	printf("centrifugal = ");	
	printf("%d\n", body.centrifugal);
	printf("tidal = ");	
	printf("%d\n", body.tidal);
	
	printf("x = ");
	print_vector(body.x);
	printf("x_dot = ");
	print_vector(body.x_dot);
	printf("l = ");
	print_vector(body.l);
	double b0[9];
	construct_traceless_symmetric_matrix(b0, body.b0_me);
	printf("b0 = \n");
	print_square_matrix(b0);
	double u[9];
	construct_traceless_symmetric_matrix(u, body.u_me);
	printf("u = \n");
	print_square_matrix(u);
	if (body.elements > 0)
	{
		double	bk_me_2d_array[body.elements][5];
		double	bk[9];
		for (int k = 0; k < body.elements; k++)
		{
			for (int l = 0; l < 5; l++)
			{
				bk_me_2d_array[k][l] = body.bk_me[l + (k*5)];
			}
			construct_traceless_symmetric_matrix(bk, bk_me_2d_array[k]);
			printf("b_%d = \n", k+1);	
			print_square_matrix(bk);
		}
	}

	printf("omega = ");
	print_vector(body.omega);
	printf("b = \n");
	print_square_matrix(body.b);

	return 0;
}

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

	int	dim_params_per_body_without_elements = 12;
	int	dim_state_per_body_without_elements = 19;
	int elements_total, elements_counter = 0; 
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		
		bodies[i].centrifugal 			= (bool) par[2 + 0 + dim_params_skip];
		bodies[i].tidal 					= (bool) par[2 + 1 + dim_params_skip];
		for (int j = 0; j < 3; j++)
		{
			bodies[i].omega[j] 			= par[2 + 2 + dim_params_skip + j];
		}
		bodies[i].mass 					= par[2 + 5 + dim_params_skip];
		bodies[i].I0 						= par[2 + 6 + dim_params_skip];
		bodies[i].gamma 					= par[2 + 7 + dim_params_skip];
		bodies[i].alpha 					= par[2 + 8 + dim_params_skip];
		bodies[i].eta						= par[2 + 9 + dim_params_skip];
		bodies[i].alpha_0 				= par[2 + 10 + dim_params_skip];
		bodies[i].elements				= (int) par[2 + 11 + dim_params_skip];
		if (bodies[i].elements > 0)
		{
			bodies[i].alpha_elements = (double *) malloc(bodies[i].elements * sizeof(double));
			for (int j = 0; j < bodies[i].elements; j++)
			{
				bodies[i].alpha_elements[j] 	= par[2 + 12 + dim_params_skip + j];
			}
			bodies[i].eta_elements = (double *) malloc(bodies[i].elements * sizeof(double));
			for (int j = 0; j < bodies[i].elements; j++)
			{
				bodies[i].eta_elements[j] 	= par[2 + 13 + dim_params_skip + j + bodies[i].elements - 1];
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

	// printf("here\n");

	/* for testing */	
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%s ", bodies[i].name);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].mass);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].lod);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].obl);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].psi);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].R);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].rg);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].J2);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].C22);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].lib);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].kf);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].Dt);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].tau);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].a);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].e);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].I);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].M);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].w);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].Omega);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", bodies[i].tidal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", bodies[i].centrifugal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// {
	// 	printf("Body %d\n", i+1);
	// 	printf("x = \n");
	// 	print_vector(bodies[i].x);
	// 	printf("x_dot = \n");
	// 	print_vector(bodies[i].x_dot);
	// 	printf("l = \n");
	// 	print_vector(bodies[i].l);
	// 	double b0_print[9], u_print[9];
	// 	construct_traceless_symmetric_matrix(b0_print, bodies[i].b0_me);
	// 	construct_traceless_symmetric_matrix(u_print, bodies[i].u_me);
	// 	printf("b0 = \n");
	// 	print_square_matrix(b0_print);
	// 	printf("u = \n");
	// 	print_square_matrix(u_print);
	// }
	// exit(99);

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

	for (int i = 0; i < number_of_bodies; i++)
	{
		double omega_hat[9];
		hat_map(omega_hat, bodies[i].omega);

		/* calculating components */

		// x component

		copy_vector (component_x[i], bodies[i].x_dot);

		// x_dot component and l component

		null_vector(component_x_dot[i]);
		null_vector(component_l[i]);

		for (int j = 0; j < number_of_bodies; j++)
		{
			if (j != i)
			{
				// x_dot component

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

				// l component

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

	}

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
	// 	printf("%s ", bodies[i].name);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].mass);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].lod);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].obl);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].psi);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].R);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].rg);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].J2);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].C22);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].lib);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].kf);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].Dt);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].tau);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].a);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].e);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].I);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].M);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].w);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%1.5e ", bodies[i].Omega);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", bodies[i].tidal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", bodies[i].centrifugal);
	// printf("\n");
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

int
hat_map(double x_hat[9], const double x[3])
{
	x_hat[0] = 0.0;
	x_hat[1] = -1.0 * x[2];
	x_hat[2] = x[1];
	x_hat[3] = x[2];
	x_hat[4] = 0.0;
	x_hat[5] = -1.0 * x[0];
	x_hat[6] = -1.0 * x[1];
	x_hat[7] = x[0];
	x_hat[8] = 0.0;

	return 0;
}

int
check_map(double x[3], const double x_hat[9])
{
	x[0] = -1.0 * x_hat[5];
	x[1] = x_hat[2];
	x[2] = -1.0 * x_hat[1];

	return 0;
}

int
construct_traceless_symmetric_matrix(double M[9], 
	const double M_main_elements[5])
{
	double M_11 = M_main_elements[0];
	double M_12 = M_main_elements[1];
	double M_13 = M_main_elements[2];
	double M_22 = M_main_elements[3];
	double M_23 = M_main_elements[4];

	M[0] = M_11;
	M[1] = M_12;
	M[2] = M_13;
	M[3] = M_12;
	M[4] = M_22;
	M[5] = M_23;
	M[6] = M_13;
	M[7] = M_23;
	M[8] = -1.0 * (M_11 + M_22);

	return 0;
}

int
get_main_elements_traceless_symmetric_matrix(double M_main_elements[5], 
	const double M[9])
{
	double M_11 = M[0];
	double M_12 = M[1];
	double M_13 = M[2];
	double M_22 = M[4];
	double M_23 = M[5];

	M_main_elements[0] = M_11;
	M_main_elements[1] = M_12;
	M_main_elements[2] = M_13;
	M_main_elements[3] = M_22;
	M_main_elements[4] = M_23;

	return 0;
}

double
parameter_gamma(const double G,	const double I0, 
	const double R, const double kf)
{
	return 3.0 * I0 * G / (pow(R, 5.0) * kf);
}

double
calculate_c(const cltbdy body)
{
	return body.gamma + body.alpha_0 + body.alpha;
}

int
calculate_f_tide(double f_tide[9],
				 const int id,
			 	 const cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G)
{
	/* initializing tidal force on body id as null */
	null_matrix(f_tide);

	/* loop over all bodies different from body id */
	for (int i = 0; i < number_of_bodies; i++)
	{
		if (i != id)
		{
			double relative_x[3];
			linear_combination_vector(relative_x,
				1.0, bodies[id].x,
				-1.0, bodies[i].x);
			double x_tensor_x[9];
			tensor_product(x_tensor_x, relative_x, relative_x);
			double Id[9];
			identity_matrix(Id);
			double scaled_id[9];
			scale_square_matrix(scaled_id, 
				norm_squared_vector(relative_x) / 3.0, Id);
			double x_norm_fifth = pow(norm_vector(relative_x), 5.0);
			double f_tide_component[9];
			linear_combination_square_matrix(f_tide_component, 
				 3.0 * G * bodies[i].mass / x_norm_fifth, x_tensor_x,
				-3.0 * G * bodies[i].mass / x_norm_fifth, scaled_id);
			linear_combination_square_matrix(f_tide, 
				1.0, f_tide,
				1.0, f_tide_component);
		}
	}

	return 0;
}

int
calculate_g	(double g[9], 
			 const int id,
			 const cltbdy *bodies,
			 const int number_of_bodies,
			 const double G)
{
	/* for testing */
	// null_matrix(g);
	// return 0;

	/* construct b0, and u matrices */
	double b0[9], u[9];
	construct_traceless_symmetric_matrix(b0, bodies[id].b0_me);
	construct_traceless_symmetric_matrix(u, bodies[id].u_me);

	/* prestress contribution */
	double alpha_0_b0[9];
	scale_square_matrix(alpha_0_b0, bodies[id].alpha_0, b0);

	/* calculate g without Voigt elements */
	linear_combination_square_matrix(g,
		 1.0, alpha_0_b0,
		-1.0, u);

	/* add tidal force if chosen */
	if (bodies[id].tidal == true)
	{
		/* calculate f_tide */
		double f_tide[9];
		calculate_f_tide(f_tide, id, bodies, number_of_bodies, G);
		linear_combination_square_matrix(g,
			1.0, g,
			1.0, f_tide);
	}

	/* add Voigt elements to g */
	double **bk_me_2d_array, **bk;
	if (bodies[id].elements > 0)
	{
		bk_me_2d_array 	= (double **) malloc(bodies[id].elements * sizeof(double *));
		bk 				= (double **) malloc(bodies[id].elements * sizeof(double *));
		for (int i = 0; i < bodies[id].elements; i++)
		{
			bk_me_2d_array[i] = (double *) malloc(5 * sizeof(double));
			for (int j = 0; j < 5; j++)
			{
				bk_me_2d_array[i][j] = bodies[id].bk_me[j + (i*5)];

				/* for testing */
				// printf("bk_me = %f\n",bk_me[i][j]);
			}
			bk[i] = (double *) malloc(9 * sizeof(double));
			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);

			/* for testing */
			// print_square_matrix(bk[i]);

			linear_combination_square_matrix(g,
							1.0, g,
				 bodies[id].alpha, bk[i]);
		}
	}

	/* for testing */
	// for (int i  = 0; i < elements; i++)
	// {
	// 	printf("\nbk_%d = \n", i+1);
	// 	print_square_matrix(bk[i]);
	// }

	/* freeing Voigt elements */
	if (bodies[id].elements > 0)
	{
		for (int i = 0; i < bodies[id].elements; i++) free(bk_me_2d_array[i]);
		free(bk_me_2d_array);
		for (int i = 0; i < bodies[id].elements; i++) free(bk[i]);
		free(bk);
	}

	return 0;
}

int
calculate_f_cent(double f_cent[9], const double omega[3])
{
	double omega_hat[9];
	hat_map(omega_hat, omega);
	double omega_hat_squared[9];
	square_matrix_times_square_matrix(omega_hat_squared,
		omega_hat, omega_hat);
	double trace_omega_hat_squared 
		= trace_square_matrix(omega_hat_squared);
	double Id[9];
	identity_matrix(Id);
	linear_combination_square_matrix(f_cent,
		-1.0, omega_hat_squared,
		trace_omega_hat_squared / 3.0, Id);

	return 0;
}

int
calculate_b	(const int id,
			 cltbdy *bodies,
			 const int number_of_bodies,
			 const double G)
{
	/* for testing */
	// null_matrix(b);
	// return 0;

	if (bodies[id].point_mass == true) /* if body is a point mass, b equals 0 */
	{
		null_matrix(bodies[id].b);
	}
	else if (bodies[id].centrifugal == false && bodies[id].tidal == false) /* if body is not deformable, b equals to b0 */
	{
		construct_traceless_symmetric_matrix(bodies[id].b, bodies[id].b0_me);
	}
	else
	{
		/* calculate g */
		double g[9];
		calculate_g(g, id, bodies, number_of_bodies, G);

		/* calculate c */
		double c = calculate_c(bodies[id]);

		/* calculate b */
		scale_square_matrix(bodies[id].b, 1.0 / c, g);

		/* add centrifugal force if chosen */
		if (bodies[id].centrifugal == true)
		{
			/* calculate centrifugal force */
			double f_cent[9];
			calculate_f_cent(f_cent, bodies[id].omega);
			linear_combination_square_matrix(bodies[id].b,
				1.0, bodies[id].b,
				1.0 / c, f_cent);
		}
	}

	/* for testing */
	// printf("\nb = \n");
	// print_square_matrix(b);
	// printf("\ntilde_x = \n");
	// print_vector(tilde_x);
	// printf("\nb0 = \n");
	// print_square_matrix(b0);
	// printf("\nu = \n");
	// print_square_matrix(u);
	// printf("\nG = %f\n", G);
	// printf("\nm2 = %f\n", m2);
	// printf("\ngamma = %f\n", gamma);
	// printf("\nalpha = %f\n", alpha);
	// printf("\nalpha_0 = %f\n", alpha_0);
	// exit(42);

	// printf("%d\n", id);
	// print_CelestialBody(bodies[id]);

	return 0;
}

int
calculate_inertia_tensor(double I[9], const double I0, const double b[9])
{
	double Id[9];
	identity_matrix(Id);
	linear_combination_square_matrix(I, I0, Id, -1.0 * I0, b);
	return 0;
}

int
calculate_l	(const int id,
			 cltbdy *bodies,
			 const int number_of_bodies,
			 const double G)
{
	double I[9];
	calculate_inertia_tensor(I, bodies[id].I0, bodies[id].b);
	square_matrix_times_vector(bodies[id].l, I, bodies[id].omega);
	
	return 0;
}

int
calculate_omega	(const int id,
				 cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G)
{
	/* for testing */
	// scale_vector(omega, 1.0 / I0, l);
	// return 0;

	if (bodies[id].point_mass == true)
	{
		return 0;
	}

	/* method parameters */
	int 	number_of_iterates = 5;
	double 	max_error = 1e-8;
	double	error = 1.0;
	double 	omega[3], previous_omega[3];
	copy_vector(omega, bodies[id].omega);

	for (int i = 0; i < number_of_iterates; i++)
	// while (error > max_error)
	{
		/* for testing */
		// printf("%d\n", id);
		// print_vector(omega);

		/* store omega previous value */
		copy_vector(previous_omega, omega);

		/* calculate g */
		double g[9];
		calculate_g(g, id, bodies, number_of_bodies, G);

		/* calculate c */
		double c = calculate_c(bodies[id]);

		/* calculate H = 0 */
		double Id[9];
		identity_matrix(Id);
		double aux_H_first_term[9];
		linear_combination_square_matrix(aux_H_first_term,
			 1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (bodies[id].centrifugal == true)
		{
			linear_combination_square_matrix(aux_H_first_term,
				1.0, aux_H_first_term,
				2.0 * norm_squared_vector(omega) / (3.0 * c), Id);
		}
		double H_first_term[3];
		square_matrix_times_vector(H_first_term, aux_H_first_term, omega);
		double H[3];
		linear_combination_vector(H, 
							1.0, H_first_term,
			 -1.0 / bodies[id].I0, bodies[id].l);
		double minus_H[3];
		scale_vector(minus_H, -1.0, H);

		/* calculate DH */
		double DH[9];
		linear_combination_square_matrix(DH,
			1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (bodies[id].centrifugal == true)
		{
			double omega_tensor_omega[9];
			tensor_product(omega_tensor_omega, omega, omega);
			double DH_third_term[9];
			linear_combination_square_matrix(DH_third_term,
				2.0 * norm_squared_vector(omega) / (3.0 * c), Id,
				4.0 / (3.0 * c), omega_tensor_omega);
			linear_combination_square_matrix(DH,
				1.0, DH,
				1.0, DH_third_term);
		}

		/* solving linear equation m*x=b using LU decomposition */
		gsl_matrix_view m
			= gsl_matrix_view_array (DH, 3, 3);
		gsl_vector_view b
			= gsl_vector_view_array (minus_H, 3);
		double omega_minus_previous_omega[] = { 0.0, 0.0, 0.0 };
		gsl_vector_view x
			= gsl_vector_view_array (omega_minus_previous_omega, 3);

		int s;
		gsl_permutation * p = gsl_permutation_alloc (3);
		gsl_linalg_LU_decomp (&m.matrix, p, &s);
		gsl_linalg_LU_solve (&m.matrix, p, &b.vector, &x.vector);

		linear_combination_vector(omega, 
			1.0, omega_minus_previous_omega,
			1.0, previous_omega);

		// error 
		// 	= norm_vector(omega_minus_previous_omega) / norm_vector(previous_omega);
		
		error = norm_vector(omega_minus_previous_omega);

		/* for testing */
		// printf("iter = %d error = %1.5e\n", i+1, error);

		gsl_permutation_free (p);

	}

	if (error > max_error)
	{
		fprintf(stderr, "Error: error higher than the max allowed\n");
		fprintf(stderr, "for omega calculation.\n");
		exit(99);
	}

	copy_vector(bodies[id].omega, omega);

	/* for testing */
	// printf("\nomega = \n");
	// print_vector(omega);
	// exit(42);

	return 0;
}

int
calculate_total_angular_momentum(double l_total[3],
				 		 		 const cltbdy *bodies,
			 	 		 		 const int number_of_bodies,
			 	 		 		 const double G)
{
	null_vector(l_total);

	for (int i = 0; i < number_of_bodies; i++)
	{
		double x_cross_x_dot[3];
		cross_product(x_cross_x_dot, bodies[i].x, bodies[i].x_dot);
		double l_total_component[3];
		linear_combination_vector(l_total_component,
			1.0, bodies[i].l,
			bodies[i].mass, x_cross_x_dot);
		linear_combination_vector(l_total,
			1.0, l_total,
			1.0, l_total_component);
	}

	return 0;
}

int
calculate_center_of_mass(double center_of_mass[3],
				 		 const cltbdy *bodies,
			 	 		 const int number_of_bodies,
			 	 		 const double G)
{
	double sum_of_masses = 0.0;
	double sum_of_masses_times_positions[] = {0.0,0.0,0.0};

	for (int i = 0; i < number_of_bodies; i++)
	{
		sum_of_masses += bodies[i].mass;
		double mass_times_position[3];
		scale_vector(mass_times_position, bodies[i].mass, bodies[i].x);
		linear_combination_vector(sum_of_masses_times_positions,
			1.0, sum_of_masses_times_positions,
			1.0, mass_times_position);
	}

	scale_vector(center_of_mass, 
		1.0 / sum_of_masses, sum_of_masses_times_positions);

	return 0;
}

double
calculate_J2(const double m, const double R, const double I[9])
{
	double J2;

	double I_11 = I[0];
	double I_22 = I[4];
	double I_33 = I[8];

	J2 = (2.0 * I_33 - I_11 - I_22) / (2.0 * m * R * R);

	return J2;
}

double
calculate_C22(const double m, const double R, const double I[9])
{
	double C22;

	double I_11 = I[0];
	double I_22 = I[4];

	C22 = (I_22 - I_11) / (4.0 * m * R * R);

	return C22;
}

double
calibrate_Imk2(const double rate, const double dist, 
	const double m1, const double m2, const double I0, 
	const double R,	const double omega_z, const double G)
{
	double Imk2;

	double alpha = rate;
	double r = dist;
	double a = R;
	double Omega = omega_z;

	double M = m1 + m2;
	double m = (m1 * m2) / M;

	double r2 = pow(r, 2.0);
	double r3 = pow(r, 3.0);
	double r5 = pow(r, 5.0);

	double n = sqrt((G * M) / r3);

	double n2 = pow(n, 2.0);

	double n_p = -(3.0 / 2.0) * sqrt((G * M) / r5);
	double Omega_p = -(m / (2.0 * I0)) * sqrt((G * M) / r);

	double h_1 = m * n * n_p * r2;
	double h_2 = m * n2 * r;
	double h_3 = I0 * Omega * Omega_p;
	double h_4 = (G * m1 * m2) / r2;

	double h = h_1 + h_2 + h_3 + h_4;

	double N_2 = sqrt(5.0 / (4.0 * M_PI * 24.0));

	double a_hat = (G * m2) / (2.0 * N_2 * r3);

	double a5 = pow(a, 5.0);
	double a_hat2 = pow(a_hat, 2.0);

	double omega_SD = 2.0 * (Omega - n); // Semi-diurnal

	double beta = -(5.0 * a5 * a_hat2 * omega_SD) / (32.0 * M_PI * G);

	Imk2 = -(alpha * h) / beta;

	return Imk2;
}

int
calculate_tau_v_and_tau(double tau_v_pair[2], double tau_pair[2],
	const double nu, const double Imk2, const double dist,
	const double m1, const double m2, const double kf, 
	const double omega_z, const double G)
{
	double tau_local_minus;
	double tau_local_plus;
	double tau_v_local_minus;
	double tau_v_local_plus;

	double r = dist;
	double b = Imk2;
	double Omega = omega_z;

	double M = m1 + m2;

	double r3 = pow(r, 3.0);

	double n = sqrt((G * M) / r3);

	double omega_tilde = 2.0 * (Omega - n); // Semi-diurnal

	double kf2 = pow(kf, 2.0);
	double nu2 = pow(nu, 2.0);
	double b2 = pow(b, 2.0);

	tau_local_minus = (-kf * nu) / (2.0 * b * omega_tilde) 
		- (1.0 / (2.0 * omega_tilde)) * sqrt(((kf2 * nu2) / b2) - 4.0);
	tau_local_plus = (-kf * nu) / (2.0 * b * omega_tilde) 
		+ (1.0 / (2.0 * omega_tilde)) * sqrt(((kf2 * nu2) / b2) - 4.0);

	tau_v_local_minus = nu * tau_local_minus; 
	tau_v_local_plus = nu * tau_local_plus; 

	tau_v_pair[0] = tau_v_local_minus;
	tau_v_pair[1] = tau_v_local_plus;
	tau_pair[0] = tau_local_minus;
	tau_pair[1] = tau_local_plus;

	return 0;
}

int
calculate_k2(double *re, double *im, const double sigma, 
	const double kf, const double tau_v, const double tau)
{
	double tau_e = tau - tau_v;
	double sigma2 = sigma * sigma;
	double tau2 = tau * tau;

	*re =  kf * ((1.0 + sigma2 * tau_e * tau) / (1.0 + sigma2 * tau2));
	*im = -kf * ((sigma * tau_v) / (1.0 + sigma2 * tau2));

	return 0;
}
