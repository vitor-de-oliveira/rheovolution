#include "dynamical_system.h"

// int
// field_1EB1PM(double t, const double y[], double f[],
//        		 void *params)
// {
// 	(void)(t);

// 	/* preparing parameters */

// 	double 	*par = (double *)params;

// 	double	omega_seed[] 	= { par[0], par[1], par[2] };
// 	double 	G				= par[3];
// 	double 	m1		 		= par[4];
// 	double 	m2		 		= par[5];
// 	double 	I0 		 		= par[6];
// 	double	gamma	 		= par[7];
// 	double 	alpha 	 		= par[8];
// 	double 	eta 	 		= par[9];
// 	double	alpha_0	 		= par[10];
// 	int 	elements 		= (int) par[11];
// 	double	*alpha_elements, *eta_elements;
// 	if (elements > 0)
// 	{
// 		alpha_elements 	= (double *) malloc(elements * sizeof(double));
// 		eta_elements 	= (double *) malloc(elements * sizeof(double));
// 		for (int i = 0; i < elements; i++)
// 		{
// 			alpha_elements[i] 	= par[12 + (2*i)];
// 			eta_elements[i]		= par[13 + (2*i)];

// 			/* for testing */
// 			// printf("alpha_%d = %f\n", i, alpha_elements[i]);
// 			// printf("eta_%d = %f\n", i, eta_elements[i]);
// 		}
// 	}
// 	bool	centrifugal		= (bool) par[12 + (elements * 2)];
// 	bool	tidal			= (bool) par[13 + (elements * 2)];

// 	/* preparing variables */

// 	double tilde_x[3], tilde_x_dot[3], l[3];
// 	double b0_me[9], u_me[9];
// 	double *bk_me = *(&bk_me), **bk_me_2d_array;

// 	for (int i = 0; i < 3; i++)	tilde_x[i] 		= y[0 + i];
// 	for (int i = 0; i < 3; i++) tilde_x_dot[i] 	= y[3 + i];
// 	for (int i = 0; i < 3; i++) l[i] 			= y[6 + i];
// 	for (int i = 0; i < 5; i++) b0_me[i] 		= y[9 + i];
// 	for (int i = 0; i < 5; i++) u_me[i] 		= y[14 + i];
// 	if (elements > 0)
// 	{
// 		bk_me = (double *) malloc(elements * 5 * sizeof(double));
// 		for (int i  = 0; i < elements * 5; i++)
// 		{
// 			bk_me[i] = y[19 + i];
// 		}
// 		bk_me_2d_array = (double **) malloc(elements * sizeof(double));
// 		for (int i = 0; i < elements; i++)
// 		{
// 			bk_me_2d_array[i] = (double *) malloc(5 * sizeof(double));
// 			for (int j = 0; j < 5; j++)
// 			{
// 				bk_me_2d_array[i][j] = y[19 + j + (i*5)];

// 				/* for testing */
// 				// printf("bk_me_2d_array = %f\n",bk_me_2d_array[i][j]);
// 			}
// 		}
// 	}

// 	double b0[9], u[9];
// 	construct_traceless_symmetric_matrix(b0, b0_me);
// 	construct_traceless_symmetric_matrix(u, u_me);

// 	double **bk = *(&(*(&bk)));
// 	if (elements > 0)
// 	{
// 		bk = (double **) malloc(elements * sizeof(double));
// 		for (int i = 0; i < elements; i++)
// 		{
// 			bk[i] = (double *) malloc(9 * sizeof(double));
// 			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);

// 			/* for testing */
// 			// print_square_matrix(bk[i]);
// 		}
// 	}

// 	/* calculate omega and b */
// 	double omega[3], b[9];
// 	// copy_vector(omega, omega_seed); // for testing
// 	calculate_omega(omega, omega_seed, G, m2, I0, gamma, alpha_0, 
// 		alpha, tilde_x, l, b0_me, u_me, elements, bk_me, centrifugal, tidal);
// 	calculate_b(b, G, m2, gamma, alpha_0, alpha,
// 		tilde_x, omega, b0_me, u_me, elements, bk_me,
// 		centrifugal, tidal);
	
// 	/* for testing */
// 	// printf("omega inside = \n");
// 	// print_vector(omega);
// 	// printf("b = \n");
// 	// print_square_matrix(b);
// 	// exit(42);
// 	// null_matrix(b);
// 	// null_matrix(omega);
// 	// double b_me[5];
// 	// for (int i = 0; i < 5; i++)
// 	// {
// 	// 	b_me[i] = ((double) i) * 0.00000000001;
// 	// }
// 	// construct_traceless_symmetric_matrix(b, b_me);

// 	double omega_hat[9];
// 	hat_map(omega_hat, omega);

// 	/* useful definitions */

// 	double minus_G_times_total_mass = -1.0 * G * (m1 + m2);

// 	double tilde_x_norm 		= norm_vector(tilde_x);
// 	double tilde_x_norm_cube 	= pow(tilde_x_norm, 3.0);
// 	double tilde_x_norm_fifth 	= pow(tilde_x_norm, 5.0);

// 	double tau = eta / alpha;
// 	double *tau_elements;
// 	if (elements > 0)
// 	{
// 		tau_elements = (double *) malloc(elements * sizeof(double));
// 		for (int i = 0; i < elements; i++)
// 		{
// 			tau_elements[i] = eta_elements[i] / alpha_elements[i];
// 		}			
// 	}

// 	double bx[3];
// 	square_matrix_times_vector(bx, b, tilde_x);
// 	double x_cross_bx[3];
// 	cross_product(x_cross_bx, tilde_x, bx);

// 	double lambda[9];
// 	linear_combination_square_matrix(lambda, 1.0, u, alpha, b);
// 	for (int i = 0; i < elements; i++)
// 	{
// 		linear_combination_square_matrix(lambda, 1.0, lambda, -1.0*alpha, bk[i]);
// 	}

// 	/* calculating components */

// 	// tilde_x component

// 	double component_tilde_x[] = { 0.0, 0.0, 0.0 };
// 	copy_vector (component_tilde_x, tilde_x_dot);

// 	// tilde_x_dot component

// 	double component_tilde_x_dot_1st_term[] = { 0.0, 0.0, 0.0 };
// 	scale_vector (component_tilde_x_dot_1st_term, 
// 		1.0 / tilde_x_norm_cube, tilde_x);

// 	double component_tilde_x_dot_2nd_term[] = { 0.0, 0.0, 0.0 };
// 	double tilde_x_norm_seventh = pow(tilde_x_norm, 7.0);
// 	double bx_dot_x = dot_product(bx, tilde_x);
// 	scale_vector (component_tilde_x_dot_2nd_term, 
// 		(15. * I0 * bx_dot_x) / (2. * m1 * tilde_x_norm_seventh), tilde_x);

// 	double component_tilde_x_dot_3rd_term[] = { 0.0, 0.0, 0.0 };
// 	scale_vector (component_tilde_x_dot_3rd_term, 
// 		(-3.0 * I0) / (m1 * tilde_x_norm_fifth), bx);

// 	double component_tilde_x_dot[] = { 0.0, 0.0, 0.0 };
// 	linear_combination_three_vector(component_tilde_x_dot,
// 		minus_G_times_total_mass, component_tilde_x_dot_1st_term, 
// 		minus_G_times_total_mass, component_tilde_x_dot_2nd_term, 
// 		minus_G_times_total_mass, component_tilde_x_dot_3rd_term);

// 	// l component

// 	double component_l[] = { 0.0, 0.0, 0.0 };
// 	scale_vector (component_l, 
// 		-3.0 * G * m2 * I0 / tilde_x_norm_fifth, x_cross_bx);

// 	// b0 component

// 	double component_b0[] = { 0.0, 0.0, 0.0,
// 							  0.0, 0.0, 0.0,
// 							  0.0, 0.0, 0.0 };
// 	commutator(component_b0, omega_hat, b0);
// 	double component_b0_me[] = { 0.0, 0.0, 0.0,
// 							          0.0, 0.0 };
// 	get_main_elements_traceless_symmetric_matrix(component_b0_me, component_b0);

// 	// u component

// 	double omega_hat_comm_u[9];
// 	commutator(omega_hat_comm_u, omega_hat, u);
// 	double component_u[] = { 0.0, 0.0, 0.0,
// 							 0.0, 0.0, 0.0,
// 							 0.0, 0.0, 0.0 };
// 	linear_combination_square_matrix(component_u, 
// 		1.0, omega_hat_comm_u, -1.0 / tau, lambda);
// 	double component_u_me[] = { 0.0, 0.0, 0.0,
// 							         0.0, 0.0 };
// 	get_main_elements_traceless_symmetric_matrix(component_u_me, component_u);

// 	/* for testing */
// 	// printf("\nomega_hat_comm_u = \n");
// 	// print_square_matrix(omega_hat_comm_u);
// 	// printf("\ntau = \n");
// 	// printf("%f\n", tau);
// 	// printf("\nlambda = \n");
// 	// print_square_matrix(lambda);
// 	// printf("\ncomponent_u = \n");
// 	// print_square_matrix(component_u);	
// 	// printf("\ncomponent_u_me = \n");
// 	// printf("%1.10e %1.10e %1.10e %1.10e %1.10e\n",
// 	// 	component_u_me[0], component_u_me[1], component_u_me[2],
// 	// 	component_u_me[3], component_u_me[4]);	
// 	// printf("\nu_me = \n");
// 	// printf("%1.10e %1.10e %1.10e %1.10e %1.10e\n",
// 	// 	u_me[0], u_me[1], u_me[2], u_me[3], u_me[4]);	
// 	// exit(42);

// 	// bk components

// 	double **component_bk_me = *(&(*(&component_bk_me)));
// 	if (elements > 0)
// 	{
// 		component_bk_me = (double **) malloc(elements * sizeof(double));
// 		for (int i = 0; i < elements; i++)
// 		{
// 			component_bk_me[i] = (double *) malloc(5 * sizeof(double));

// 			double omega_hat_comm_bk[9];
// 			commutator(omega_hat_comm_bk, omega_hat, bk[i]);
// 			double minus_bk_over_tau_elements[9];
// 			scale_square_matrix(minus_bk_over_tau_elements, 
// 				-1.0 / tau_elements[i], bk[i]);
// 			double lambda_over_eta_elements[9];
// 			scale_square_matrix(lambda_over_eta_elements,
// 				1.0 / eta_elements[i], lambda);

// 			/* for testing */
// 			// printf("\ntau_%d = %f\n", i, tau_elements[i]);
// 			// printf("\nminus_1_over_tau_%d = %f\n", i, -1.0 / tau_elements[i]);
// 			// printf("\nomega_hat_comm_bk = \n");
// 			// print_square_matrix(omega_hat_comm_bk);
// 			// printf("\nminus_bk_over_tau_elements = \n");
// 			// print_square_matrix(minus_bk_over_tau_elements);
// 			// printf("\nlambda_over_eta_elements = \n");
// 			// print_square_matrix(lambda_over_eta_elements);

// 			double component_bk[] = { 0.0, 0.0, 0.0,
// 							 		  0.0, 0.0, 0.0,
// 							          0.0, 0.0, 0.0 };
// 			linear_combination_three_square_matrix(component_bk,
// 				1.0, omega_hat_comm_bk,
// 				1.0, minus_bk_over_tau_elements,
// 				1.0, lambda_over_eta_elements);

// 			get_main_elements_traceless_symmetric_matrix(component_bk_me[i],
// 				component_bk);

// 			/* for testing */
// 			// printf("\nalpha_%d = %f\n", i, alpha_elements[i]);
// 			// printf("\ntau_%d = %f\n", i, tau_elements[i]);
// 			// printf("\ncomponent_bk = \n");
// 			// print_square_matrix(component_bk);
// 		}
// 	}

// 	/* writing components */	

// 	for (int i = 0; i < 3; i++) f[i] 		= component_tilde_x[i];
// 	for (int i = 0; i < 3; i++) f[3 + i] 	= component_tilde_x_dot[i];
// 	for (int i = 0; i < 3; i++) f[6 + i] 	= component_l[i];
// 	for (int i = 0; i < 5; i++) f[9 + i]	= component_b0_me[i];
// 	for (int i = 0; i < 5; i++) f[14 + i]	= component_u_me[i];
// 	for (int i = 0; i < elements; i++)
// 	{
// 		for (int j = 0; j < 5; j++)
// 		{
// 			f[19 + j + (i*5)] = component_bk_me[i][j];
// 		}
// 	}

// 	/* for testing */
// 	// printf("tilde_x = \n");
// 	// print_vector(tilde_x);
// 	// printf("\ntilde_x_dot = \n");
// 	// print_vector(tilde_x_dot);
// 	// printf("\nl = \n");
// 	// print_vector(l);
// 	// printf("\nb0 = \n");
// 	// print_square_matrix(b0);
// 	// printf("\nb = \n");
// 	// print_square_matrix(b);
// 	// printf("\nu = \n");
// 	// print_square_matrix(u);
// 	// for (int i  = 0; i < elements; i++)
// 	// {
// 	// 	printf("\nbk_%d = \n", i+1);
// 	// 	print_square_matrix(bk[i]);
// 	// }
// 	// printf("\nomega = \n");
// 	// print_vector(omega);
// 	// printf("\nG = %e\n", G);
// 	// printf("\nm1 = %e\n", m1);
// 	// printf("\nm2 = %e\n", m2);
// 	// printf("\nI0 = %e\n", I0);
// 	// printf("\ngamma = %e\n", gamma);
// 	// printf("\nalpha = %e\n", alpha);
// 	// printf("\neta = %e\n", eta);
// 	// printf("\nalpha_0 = %e\n", alpha_0);
// 	// for (int i  = 0; i < elements; i++)
// 	// {
// 	// 	printf("\nalpha_%d = %e\n", i+1, alpha_elements[i]);
// 	// 	printf("\neta%d = %e\n", i+1, eta_elements[i]);
// 	// }
// 	// printf("\ncomponent_tilde_x = \n");
// 	// print_vector(component_tilde_x);
// 	// printf("\ncomponent_tilde_x_dot = \n");
// 	// print_vector(component_tilde_x_dot);
// 	// printf("\ncomponent_l = \n");
// 	// print_vector(component_l);
// 	// printf("\ncomponent_u = \n");
// 	// print_square_matrix(component_u);	
// 	// printf("\ncomponent_b0 = \n");
// 	// print_square_matrix(component_b0);
// 	// printf("\nlambda = \n");
// 	// print_square_matrix(lambda);
// 	// exit(42);
	
// 	/* freeing Voigt elements */
// 	if (elements > 0)
// 	{
// 		free(alpha_elements);
// 		free(eta_elements);
// 		free(tau_elements);
// 		free(bk_me);
// 		for (int i = 0; i < elements; i++) free(bk_me_2d_array[i]);
// 		free(bk_me_2d_array);
// 		for (int i = 0; i < elements; i++) free(bk[i]);
// 		free(bk);
// 		for (int i = 0; i < elements; i++) free(component_bk_me[i]);
// 		free(component_bk_me);
// 	}

// 	return GSL_SUCCESS;
// }

int
field_GV(double t, 
		 const double y[],
		 double f[],
       	 void *params)
{
	/* reinforce autonomous trait of system */
	(void)(t);

	/* preparing variables and parameters */

	double 	*par = (double *)params;

	// for (int i = 0; i < 19 * 2; i++)
	// {
	// 	printf("y[%d] = %1.10e\n", i, y[i]);
	// }
	// for (int i = 0; i < 12 * 2 + 2; i++)
	// {
	// 	printf("%f\n", par[i]);
	// }
	// exit(98);

	double 	G					= par[0];
	int		number_of_bodies	= (int) par[1];

	cltbdy 	*body;
	body = malloc (number_of_bodies * sizeof(cltbdy));

	int	dim_params_per_body_without_elements = 12;
	int	elements_counter = 0; 
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_params_skip = i * dim_params_per_body_without_elements + 2 * elements_counter;
		
		body[i].centrifugal 			= (bool) par[2 + 0 + dim_params_skip];
		body[i].tidal 					= (bool) par[2 + 1 + dim_params_skip];
		for (int j = 0; j < 3; j++)
		{
			body[i].omega[j] 			= par[2 + 2 + dim_params_skip + j];
		}
		body[i].mass 					= par[2 + 5 + dim_params_skip];
		body[i].I0 						= par[2 + 6 + dim_params_skip];
		body[i].gamma 					= par[2 + 7 + dim_params_skip];
		body[i].alpha 					= par[2 + 8 + dim_params_skip];
		body[i].eta						= par[2 + 9 + dim_params_skip];
		body[i].alpha_0 				= par[2 + 10 + dim_params_skip];
		body[i].elements				= (int) par[2 + 11 + dim_params_skip];
		if (body[i].elements > 0)
		{
			body[i].alpha_elements = malloc(body[i].elements * sizeof(double));
			for (int j = 0; j < body[i].elements; j++)
			{
				body[i].alpha_elements[j] 	= par[2 + 12 + dim_params_skip + j];
			}
			body[i].eta_elements = malloc(body[i].elements * sizeof(double));
			for (int j = 0; j < body[i].elements; j++)
			{
				body[i].eta_elements[j] 	= par[2 + 13 + dim_params_skip + j + body[i].elements - 1];
			}
		}
		elements_counter += body[i].elements;
	}

	int	dim_state_per_body_without_elements = 19;
	elements_counter = 0; 
	for (int i = 0; i < number_of_bodies; i++)
	{
		int	dim_state_skip = i * dim_state_per_body_without_elements + 5 * elements_counter;
		
		for (int j = 0; j < 3; j++)
		{
			body[i].x[j] 		= y[0 + dim_state_skip + j];
			body[i].x_dot[j] 	= y[3 + dim_state_skip + j];
			body[i].l[j] 		= y[6 + dim_state_skip + j];
		}
		for (int j = 0; j < 5; j++)
		{
			body[i].b0_me[j] 	= y[9 + dim_state_skip + j];
			body[i].u_me[j] 	= y[14 + dim_state_skip + j];
		}
		if (body[i].elements > 0)
		{
			body[i].bk_me = malloc(5 * body[i].elements * sizeof(double));
			for (int j = 0; j < 5 * body[i].elements; j++)
			{
				body[i].bk_me[j] = y[19 + dim_state_skip + j];
			}
		}

		elements_counter += body[i].elements;
	}

	int elements_total = elements_counter;

	/* calculate omega and b for every body */
	for (int i = 0; i < number_of_bodies; i++)
	{
		// printf("omega of body %d before\n", i+1);
		// print_vector(body[i].omega);
		calculate_omega(i, body, number_of_bodies, G);
		// printf("omega of body %d after\n", i+1);
		// print_vector(body[i].omega);
		calculate_b(i, body, number_of_bodies, G);
	}

	// printf("here\n");

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
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", body[i].tidal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", body[i].centrifugal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// {
	// 	printf("Body %d\n", i+1);
	// 	printf("x = \n");
	// 	print_vector(body[i].x);
	// 	printf("x_dot = \n");
	// 	print_vector(body[i].x_dot);
	// 	printf("l = \n");
	// 	print_vector(body[i].l);
	// 	double b0_print[9], u_print[9];
	// 	construct_traceless_symmetric_matrix(b0_print, body[i].b0_me);
	// 	construct_traceless_symmetric_matrix(u_print, body[i].u_me);
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

	double **component_x;
	double **component_x_dot;
	double **component_l;
	double **component_b0_me;
	double **component_u_me;
	double ***component_bk_me = *(&component_bk_me);

	component_x	= (double **) malloc(number_of_bodies * sizeof(double));
	component_x_dot	= (double **) malloc(number_of_bodies * sizeof(double));
	component_l	= (double **) malloc(number_of_bodies * sizeof(double));
	component_b0_me	= (double **) malloc(number_of_bodies * sizeof(double));
	component_u_me	= (double **) malloc(number_of_bodies * sizeof(double));
	if (elements_total > 0)
	{
		component_bk_me	= (double ***) malloc(number_of_bodies * sizeof(double));
	}
	for (int i = 0; i < number_of_bodies; i++)
	{
		component_x[i] = (double *) malloc(3 * sizeof(double));
		component_x_dot[i] = (double *) malloc(3 * sizeof(double));
		component_l[i] = (double *) malloc(3 * sizeof(double));
		component_b0_me[i] = (double *) malloc(5 * sizeof(double));
		component_u_me[i] = (double *) malloc(5 * sizeof(double));
		if (body[i].elements > 0)
		{
			component_bk_me[i] = (double **) malloc(body[i].elements * sizeof(double));

			for (int j = 0; j < body[i].elements; j++)
			{
				component_bk_me[i][j] = (double *) malloc(5 * sizeof(double));
			}
		}
	}

	for (int i = 0; i < number_of_bodies; i++)
	{
		double omega_hat[9];
		hat_map(omega_hat, body[i].omega);

		/* useful definitions */

		// double minus_G_times_total_mass = -1.0 * G * (m1 + m2);

		// double tilde_x_norm 		= norm_vector(tilde_x);
		// double tilde_x_norm_cube 	= pow(tilde_x_norm, 3.0);
		// double tilde_x_norm_fifth 	= pow(tilde_x_norm, 5.0);

		// double tau = eta / alpha;
		// double *tau_elements;
		// if (elements > 0)
		// {
		// 	tau_elements = (double *) malloc(elements * sizeof(double));
		// 	for (int i = 0; i < elements; i++)
		// 	{
		// 		tau_elements[i] = eta_elements[i] / alpha_elements[i];
		// 	}			
		// }

		// double bx[3];
		// square_matrix_times_vector(bx, b, tilde_x);
		// double x_cross_bx[3];
		// cross_product(x_cross_bx, tilde_x, bx);

		/* calculating components */

		// x component

		copy_vector (component_x[i], body[i].x_dot);

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
					1.0, body[i].x,
					-1.0, body[j].x);

				double x_relative_norm 		= norm_vector(relative_x);
				double x_relative_norm_cube = pow(x_relative_norm, 3.0);
				double component_x_dot_1st_term[] = { 0.0, 0.0, 0.0 };
				scale_vector (component_x_dot_1st_term, 
					-1.0 * body[j].mass / x_relative_norm_cube, relative_x);

				double x_relative_norm_seventh = pow(x_relative_norm, 7.0);
				double scaled_sum_bi_bj[9];
				linear_combination_square_matrix(scaled_sum_bi_bj,
					body[j].mass * body[i].I0, body[i].b,
					body[i].mass * body[j].I0, body[j].b);
				double scaled_sum_bi_bj_x[3];
				square_matrix_times_vector(scaled_sum_bi_bj_x, 
					scaled_sum_bi_bj, relative_x);
				double scaled_sum_bi_bj_x_dot_x = 
					dot_product(scaled_sum_bi_bj_x, relative_x);
				double component_x_dot_2nd_term[] = { 0.0, 0.0, 0.0 };
				scale_vector (component_x_dot_2nd_term, 
					(-15.0 * scaled_sum_bi_bj_x_dot_x) / (2. * body[i].mass * x_relative_norm_seventh), 
					relative_x);

				double x_relative_norm_fifth = pow(x_relative_norm, 5.0);
				double component_x_dot_3rd_term[] = { 0.0, 0.0, 0.0 };
				scale_vector (component_x_dot_3rd_term, 
					3.0 / (body[i].mass * x_relative_norm_fifth), scaled_sum_bi_bj_x);

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
				square_matrix_times_vector(bx, body[i].b, relative_x);
				double x_cross_bx[3];
				cross_product(x_cross_bx, relative_x, bx);
				double j_component_l[] = { 0.0, 0.0, 0.0 };
				scale_vector (j_component_l, 
					(-3.0 * G * body[j].mass * body[i].I0) / x_relative_norm_fifth, x_cross_bx);

				linear_combination_vector(component_l[i],
					1.0, component_l[i],
					1.0, j_component_l);
			}
		} // end summation on number of bodies

		// b0 component

		double b0[9];
		construct_traceless_symmetric_matrix(b0, body[i].b0_me);
		double component_b0[] = { 0.0, 0.0, 0.0,
								  0.0, 0.0, 0.0,
								  0.0, 0.0, 0.0 };
		commutator(component_b0, omega_hat, b0);
		get_main_elements_traceless_symmetric_matrix(component_b0_me[i], component_b0);

		// u component

		double u[9];
		construct_traceless_symmetric_matrix(u, body[i].u_me);
		double omega_hat_comm_u[9];
		commutator(omega_hat_comm_u, omega_hat, u);
		double tau = body[i].eta / body[i].alpha;
		double lambda[9];
		linear_combination_square_matrix(lambda, 
			1.0, u, body[i].alpha, body[i].b);

		double **bk_me_2d_array, **bk = *(&bk);
		if (body[i].elements > 0)
		{
			bk_me_2d_array 	= (double **) malloc(body[i].elements * sizeof(double));
			bk 				= (double **) malloc(body[i].elements * sizeof(double));
			for (int k = 0; k < body[i].elements; k++)
			{
				bk_me_2d_array[k] = (double *) malloc(5 * sizeof(double));
				for (int l = 0; l < 5; l++)
				{
					bk_me_2d_array[k][l] = body[i].bk_me[l + (k*5)];

					/* for testing */
					// printf("bk_me = %f\n",bk_me[i][j]);
				}
				bk[k] = (double *) malloc(9 * sizeof(double));
				construct_traceless_symmetric_matrix(bk[k], bk_me_2d_array[k]);
				free(bk_me_2d_array[k]);
			}
			free(bk_me_2d_array);
		}

		for (int k = 0; k < body[i].elements; k++)
		{
			linear_combination_square_matrix(lambda, 
				1.0, lambda, -1.0 * body[i].alpha, bk[k]);
		}
		double component_u[] = { 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0,
								 0.0, 0.0, 0.0 };
		linear_combination_square_matrix(component_u, 
			1.0, omega_hat_comm_u, -1.0 / tau, lambda);
		get_main_elements_traceless_symmetric_matrix(component_u_me[i], component_u);

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

		// bk components

		for (int k = 0; k < body[i].elements; k++)
		{
			double omega_hat_comm_bk[9];
			commutator(omega_hat_comm_bk, omega_hat, bk[k]);
			double tau_elements = 
				body[i].eta_elements[k] / body[i].alpha_elements[k];
			double minus_bk_over_tau_elements[9];
			scale_square_matrix(minus_bk_over_tau_elements, 
				-1.0 / tau_elements, bk[k]);
			double lambda_over_eta_elements[9];
			scale_square_matrix(lambda_over_eta_elements,
				1.0 / body[i].eta_elements[k], lambda);

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

		if (body[i].elements > 0)
		{
			for (int k = 0; k < body[i].elements; k++)
			{
				free(bk[k]);
			}
			free(bk);
		}
	}

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
		for (int k = 0; k < body[i].elements; k++)
		{
			for (int l = 0; l < 5; l++)
			{
				f[19 + dim_state_skip + 5*k + l] = component_bk_me[i][k][l];
			}
		}
		elements_counter += body[i].elements;
	}

	/* for testing */
	// for (int i = 0; i < number_of_bodies * 19; i++)
	// {
	// 	printf("f[%d] = %1.5e\n", i, f[i]);
	// }

	for (int i = 0; i < number_of_bodies; i++)
	{
		free(component_x[i]);
		free(component_x_dot[i]);
		free(component_l[i]);
		free(component_b0_me[i]);
		free(component_u_me[i]);
		if (body[i].elements > 0)
		{
			for (int j = 0; j < body[i].elements; j++)
			{
				free(component_bk_me[i][j]);
			}
			free(component_bk_me[i]);
		}
	}
	free(component_x);
	free(component_x_dot);
	free(component_l);
	free(component_b0_me);
	free(component_u_me);
	if (elements_total > 0)
	{
		free(component_bk_me);
	}

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
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", body[i].tidal);
	// printf("\n");
	// for (int i = 0; i < number_of_bodies; i++)
	// 	printf("%d ", body[i].centrifugal);
	// printf("\n");
	// exit(99);
	
	/* free Voigt elements */
	for (int i = 0; i < number_of_bodies; i++)
	{
		if (body[i].elements > 0)
		{
			free(body[i].alpha_elements);
			free(body[i].eta_elements);	
			free(body[i].bk_me);
		}
	}

	/* free array of celestial bodies */
	free(body);

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
			 	 const cltbdy *body,
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
				1.0, body[id].x,
				-1.0, body[i].x);
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
				 3.0 * G * body[i].mass / x_norm_fifth, x_tensor_x,
				-3.0 * G * body[i].mass / x_norm_fifth, scaled_id);
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
			 const cltbdy *body,
			 const int number_of_bodies,
			 const double G)
{
	/* for testing */
	// null_matrix(g);
	// return 0;

	/* construct b0, and u matrices */
	double b0[9], u[9];
	construct_traceless_symmetric_matrix(b0, body[id].b0_me);
	construct_traceless_symmetric_matrix(u, body[id].u_me);

	/* prestress contribution */
	double alpha_0_b0[9];
	scale_square_matrix(alpha_0_b0, body[id].alpha_0, b0);

	/* calculate g without Voigt elements */
	linear_combination_square_matrix(g,
		 1.0, alpha_0_b0,
		-1.0, u);

	/* add tidal force if chosen */
	if (body[id].tidal == true)
	{
		/* calculate f_tide */
		double f_tide[9];
		calculate_f_tide(f_tide, id, body, number_of_bodies, G);
		linear_combination_square_matrix(g,
			1.0, g,
			1.0, f_tide);
	}

	/* add Voigt elements to g */
	double **bk_me_2d_array, **bk;
	if (body[id].elements > 0)
	{
		bk_me_2d_array 	= (double **) malloc(body[id].elements * sizeof(double));
		bk 				= (double **) malloc(body[id].elements * sizeof(double));
		for (int i = 0; i < body[id].elements; i++)
		{
			bk_me_2d_array[i] = (double *) malloc(5 * sizeof(double));
			for (int j = 0; j < 5; j++)
			{
				bk_me_2d_array[i][j] = body[id].bk_me[j + (i*5)];

				/* for testing */
				// printf("bk_me = %f\n",bk_me[i][j]);
			}
			bk[i] = (double *) malloc(9 * sizeof(double));
			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);

			/* for testing */
			// print_square_matrix(bk[i]);

			linear_combination_square_matrix(g,
							1.0, g,
				 body[id].alpha, bk[i]);
		}
	}

	/* for testing */
	// for (int i  = 0; i < elements; i++)
	// {
	// 	printf("\nbk_%d = \n", i+1);
	// 	print_square_matrix(bk[i]);
	// }

	/* freeing Voigt elements */
	if (body[id].elements > 0)
	{
		for (int i = 0; i < body[id].elements; i++) free(bk_me_2d_array[i]);
		free(bk_me_2d_array);
		for (int i = 0; i < body[id].elements; i++) free(bk[i]);
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
			 cltbdy *body,
			 const int number_of_bodies,
			 const double G)
{
	/* for testing */
	// null_matrix(b);
	// return 0;

	if (body[id].point_mass == true) /* if body is a point mass, b equals 0 */
	{
		null_matrix(body[id].b);
	}
	else if (body[id].centrifugal == false && body[id].tidal == false) /* if body is not deformable, b equals to b0 */
	{
		construct_traceless_symmetric_matrix(body[id].b, body[id].b0_me);
	}
	else
	{
		/* calculate g */
		double g[9];
		calculate_g(g, id, body, number_of_bodies, G);

		/* calculate c */
		double c = calculate_c(body[id]);

		/* calculate b */
		scale_square_matrix(body[id].b, 1.0 / c, g);

		/* add centrifugal force if chosen */
		if (body[id].centrifugal == true)
		{
			/* calculate centrifugal force */
			double f_cent[9];
			calculate_f_cent(f_cent, body[id].omega);
			linear_combination_square_matrix(body[id].b,
				1.0, body[id].b,
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
			 cltbdy *body,
			 const int number_of_bodies,
			 const double G)
{
	double I[9];
	calculate_inertia_tensor(I, body[id].I0, body[id].b);
	square_matrix_times_vector(body[id].l, I, body[id].omega);
	
	return 0;
}

int
calculate_omega	(const int id,
				 cltbdy *body,
			 	 const int number_of_bodies,
			 	 const double G)
{
	/* for testing */
	// scale_vector(omega, 1.0 / I0, l);
	// return 0;

	if (body[id].point_mass == true)
	{
		return 0;
	}

	/* method parameters */
	int 	number_of_iterates = 5;
	double 	max_error = 1e-8;
	double	error = 1.0;
	double 	omega[3], previous_omega[3];
	copy_vector(omega, body[id].omega);

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
		calculate_g(g, id, body, number_of_bodies, G);

		/* calculate c */
		double c = calculate_c(body[id]);

		/* calculate H = 0 */
		double Id[9];
		identity_matrix(Id);
		double aux_H_first_term[9];
		linear_combination_square_matrix(aux_H_first_term,
			 1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (body[id].centrifugal == true)
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
			 -1.0 / body[id].I0, body[id].l);
		double minus_H[3];
		scale_vector(minus_H, -1.0, H);

		/* calculate DH */
		double DH[9];
		linear_combination_square_matrix(DH,
			1.0, Id,
			-1.0 / c, g);
		/* add centrifugal term if chosen */
		if (body[id].centrifugal == true)
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

	copy_vector(body[id].omega, omega);

	/* for testing */
	// printf("\nomega = \n");
	// print_vector(omega);
	// exit(42);

	return 0;
}

int
total_angular_momentum	(double l_total[3],
				 		 const cltbdy *body,
			 	 		 const int number_of_bodies,
			 	 		 const double G)
{
	null_vector(l_total);

	for (int i = 0; i < number_of_bodies; i++)
	{
		double x_cross_x_dot[3];
		cross_product(x_cross_x_dot, body[i].x, body[i].x_dot);
		double l_total_component[3];
		linear_combination_vector(l_total_component,
			1.0, body[i].l,
			body[i].mass, x_cross_x_dot);
		linear_combination_vector(l_total,
			1.0, l_total,
			1.0, l_total_component);
	}

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