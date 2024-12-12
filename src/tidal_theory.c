#include "tidal_theory.h"

int
calculate_inertia_tensor(double I[9], const double I0, const double b[9])
{
	double Id[9];
	identity_matrix(Id);
	linear_combination_square_matrix(I, I0, Id, -1.0 * I0, b);
	return 0;
}

int
body_frame_deformation_from_stokes_coefficients	(double B[9],
												 const cltbdy body)
{
	double div = 3.0 * body.rg - 2.0 * body.J2;

	B[0] = (body.J2 + 6.0 * body.C22) / div;
	B[1] = (6.0 * body.S22) / div;
	B[2] = (3.0 * body.C21) / div;
	B[3] = B[1];
	B[4] = (body.J2 - 6.0 * body.C22) / div;
	B[5] = (3.0 * body.S21) / div;
	B[6] = B[2];
	B[7] = B[5];
	B[8] = -1.0 * (B[0] + B[4]); // B[8] = (-2.0 * body.J2) / div

	return 0;
}

double
parameter_gamma_0	(const double G,
					 const double I0, 
					 const double R,
					 const double k0)
{
	return 3.0 * I0 * G / (pow(R, 5.0) * k0);
}

int
calculate_Y_and_Y_transpose(cltbdy *body)
{
	rotation_matrix_from_quaternion((*body).Y, (*body).q);
	transpose_square_matrix((*body).Y_trans, (*body).Y);
	return 0;
}

int
calculate_bs_me(cltbdy *body)
{
	double B_stokes[9];
	construct_traceless_symmetric_matrix(B_stokes, (*body).Bs_me);
	double b_stokes[9];
	square_matrix_times_square_matrix(b_stokes, (*body).Y, B_stokes);
	square_matrix_times_square_matrix(b_stokes, b_stokes, (*body).Y_trans);
	get_main_elements_traceless_symmetric_matrix((*body).bs_me, b_stokes);

	return 0;
}

int
calculate_p_me(cltbdy *body)
{
	double P[9];
	construct_traceless_symmetric_matrix(P, (*body).P_me);
	double p[9];
	square_matrix_times_square_matrix(p, (*body).Y, P);
	square_matrix_times_square_matrix(p, p, (*body).Y_trans);
	get_main_elements_traceless_symmetric_matrix((*body).p_me, p);

	return 0;
}

double
calculate_c(const cltbdy body)
{
	return body.gamma_0 + body.alpha;
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
			double scaled_Id[9];
			scale_square_matrix(scaled_Id, 
				norm_squared_vector(relative_x) / 3.0, Id);
			double x_norm_fifth = pow(norm_vector(relative_x), 5.0);
			double f_tide_component[9];
			linear_combination_square_matrix(f_tide_component, 
				 3.0 * G * bodies[i].mass / x_norm_fifth, x_tensor_x,
				-3.0 * G * bodies[i].mass / x_norm_fifth, scaled_Id);
			linear_combination_square_matrix(f_tide, 
				1.0, f_tide,
				1.0, f_tide_component);
		}
	}

	return 0;
}

int
calculate_f_ps	(double f_ps[9],
			 	 const cltbdy body)
{
	construct_traceless_symmetric_matrix(f_ps, body.p_me);

	return 0;
}

int
calculate_f_rheo(double f_rheo[9],
			 	 const cltbdy body)
{
	double b_eta[9];
	construct_traceless_symmetric_matrix(b_eta, body.b_eta_me);
	scale_square_matrix(f_rheo, body.alpha, b_eta);
	if (body.elements > 0)
	{
		double bk_me_2d_array[body.elements][5];
		double bk[body.elements][9];
		for (int i = 0; i < body.elements; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				bk_me_2d_array[i][j] = body.bk_me[j + (i*5)];
			}
			construct_traceless_symmetric_matrix(bk[i], bk_me_2d_array[i]);
			linear_combination_square_matrix(f_rheo,
				1.0, f_rheo,
				body.alpha, bk[i]);
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
	double f_rheo[9], f_ps[9], f_tide[9];
	null_matrix(f_rheo);
	null_matrix(f_ps);
	null_matrix(f_tide);
	
	/* f_rheo (rheology) */
	if (bodies[id].deformable == true)
	{
		calculate_f_rheo(f_rheo, bodies[id]);
	}

	/* f_ps (prestress) */
	if (bodies[id].prestress == true)
	{
		calculate_f_ps(f_ps, bodies[id]);
	}

	/* f_tide (tides) */
	if (bodies[id].tidal == true)
	{
		calculate_f_tide(f_tide, id, bodies, number_of_bodies, G);
	}

	/* calculate g */
	linear_combination_three_square_matrix(g,
		1.0, f_rheo,
		1.0, f_ps,
		1.0, f_tide);

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
calculate_F_cent_mean	(double F_cent_mean[9],
				 	  	 const double mean_omega)
{
	double a = (mean_omega * mean_omega) / 3.0;
	double M[] = {1.0, 0.0, 0.0,
				  0.0, 1.0, 0.0,
				  0.0, 0.0, -2.0};
	scale_square_matrix(F_cent_mean, a, M);

	return 0;
}

int
calculate_b	(const int id,
			 cltbdy *bodies,
			 const int number_of_bodies,
			 const double G)
{
	/**
	 * if body is a point mass, b is NAN 
	 * if body is not deformable, b equals to bs
	**/
	if (bodies[id].point_mass == true)
	{
		nan_matrix(bodies[id].b);
	}
	else if (bodies[id].deformable == false)
	{
		construct_traceless_symmetric_matrix(bodies[id].b, 
			bodies[id].bs_me);
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

	return 0;
}

int
calculate_l	(cltbdy *body)
{
	double I[9];
	calculate_inertia_tensor(I, (*body).I0, (*body).b);
	square_matrix_times_vector((*body).l, I, (*body).omega);
	
	return 0;
}

int
calculate_omega	(const int id,
				 cltbdy *bodies,
			 	 const int number_of_bodies,
			 	 const double G)
{
	/**
	 * if body is a point mass, omega is constant
	**/
	if (bodies[id].point_mass == true)
	{
		return 0;
	}

	/* method parameters */
	int		counter_iterates = 0;
	int		max_number_of_iterates = 100;
	double 	max_error = 1.0e-8;
	double	error = 1.0;
	double 	omega[3], previous_omega[3];
	copy_vector(omega, bodies[id].omega);

	/* auxiliary variable */
	double Id[9];
	identity_matrix(Id);

	// for (int i = 0; i < max_number_of_iterates; i++) // option 1
	while (error > max_error) // option 2
	{
		counter_iterates++;
		if (counter_iterates == max_number_of_iterates) // option 2
		{
			fprintf(stderr, "Error: maximum number of iterates reached\n");
			fprintf(stderr, "for omega calculation in body %d.\n", id + 1);
			exit(99);			
		}

		/* store previous value of omega */
		copy_vector(previous_omega, omega);

		double	H[3], DH[9];

		if (bodies[id].deformable == false)
		{
			double bs_i[9];
			construct_traceless_symmetric_matrix(bs_i, 
				bodies[id].bs_me);

			/* calculate H = 0 */
			double aux_H_first_term[9];
			linear_combination_square_matrix(aux_H_first_term,
				 1.0, Id,
				-1.0, bs_i);
			double H_first_term[3];
			square_matrix_times_vector(H_first_term, 
				aux_H_first_term, omega);
			linear_combination_vector(H, 
								 1.0, H_first_term,
				-1.0 / bodies[id].I0, bodies[id].l);

			/* calculate DH */
			linear_combination_square_matrix(DH,
				 1.0, Id,
				-1.0, bs_i);
		}
		else
		{
			/* calculate g */
			double g[9];
			calculate_g(g, id, bodies, number_of_bodies, G);

			/* calculate c */
			double c = calculate_c(bodies[id]);

			/* calculate H = 0 */
			double aux_H_first_term[9];
			linear_combination_square_matrix(aux_H_first_term,
				1.0, Id,
				-1.0 / c, g);

			/* add centrifugal term in H if chosen */
			if (bodies[id].centrifugal == true)
			{
				linear_combination_square_matrix(aux_H_first_term,
					1.0, aux_H_first_term,
					2.0 * norm_squared_vector(omega) / (3.0 * c), Id);
			}
			double H_first_term[3];
			square_matrix_times_vector(H_first_term, 
				aux_H_first_term, omega);
			linear_combination_vector(H, 
								1.0, H_first_term,
				-1.0 / bodies[id].I0, bodies[id].l);

			/* calculate DH */
			linear_combination_square_matrix(DH,
				1.0, Id,
				-1.0 / c, g);

			/* add centrifugal term to DH if chosen */
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
		}

		double minus_H[3];
		scale_vector(minus_H, -1.0, H);

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

		error = norm_vector(omega_minus_previous_omega) / 
					norm_vector(previous_omega);
		
		// error = norm_vector(omega_minus_previous_omega);

		/* for testing */
		// printf("iter = %d error = %1.5e\n", counter_iterates, error);

		gsl_permutation_free (p);

	}

	// if (error > max_error) // option 1
	// {
	// 	fprintf(stderr, "Error: error value higher than the max allowed\n");
	// 	fprintf(stderr, "for omega calculation in body %d.\n", id + 1);
	// 	fprintf(stderr, "Value = %1.5e.\n", error);
	// 	fprintf(stderr, "Max value = %1.5e.\n", max_error);
	// 	exit(99);
	// }

	copy_vector(bodies[id].omega, omega);

	return 0;
}

int
calculate_total_angular_momentum(double l_total[3],
				 		 		 const cltbdy *bodies,
			 	 		 		 const int number_of_bodies)
{
	null_vector(l_total);

	for (int i = 0; i < number_of_bodies; i++)
	{
		double x_cross_x_dot[3];
		cross_product(x_cross_x_dot, bodies[i].x, bodies[i].x_dot);
		double l_total_component[3];
		scale_vector(l_total_component, bodies[i].mass, x_cross_x_dot);
		if (bodies[i].point_mass == false)
		{
			linear_combination_vector(l_total_component,
				1.0, l_total_component,
				1.0, bodies[i].l);
		}
		linear_combination_vector(l_total,
			1.0, l_total,
			1.0, l_total_component);
	}

	return 0;
}

double
total_energy_without_deformation(const cltbdy *bodies,
			 	 			 	 const int number_of_bodies,
							 	 const double G)
{
	double E_translational = 0.0;
	if (number_of_bodies != 1)
	{
		for (int i = 0; i < number_of_bodies; i++)
		{
			E_translational += bodies[i].mass 
				* norm_squared_vector(bodies[i].x_dot);
		}
		E_translational *= 0.5;
	}

	double E_gravitational = 0.0;
	if (number_of_bodies > 1)
	{
		for (int i = 0; i < number_of_bodies - 1; i++)
		{
			for (int j = i + 1; j < number_of_bodies; j++)
			{
				double x_ij[3];
				linear_combination_vector(x_ij,
					1.0, bodies[i].x,
					-1.0, bodies[j].x);
				E_gravitational += bodies[i].mass * bodies[j].mass
					/ norm_vector(x_ij);
			}
		}
		E_gravitational *= -1.0 * G;
	}

	double E_rotational = 0.0;
	for (int i = 0; i < number_of_bodies; i++)
	{
		if (bodies[i].point_mass == false)
		{
			double omega_hat[9];
			hat_map(omega_hat, bodies[i].omega);
			E_rotational += 0.5 * bodies[i].I0 
				* norm_squared_square_matrix(omega_hat);
				// I HAVE TO ADD DEFORMATION HERE AS WELL
		}
	}

	return E_translational + E_gravitational + E_rotational;
}

double
dissipation_function(const cltbdy body)
{
	double D = 0.0;

	double b_eta[9];
	construct_traceless_symmetric_matrix(b_eta, body.b_eta_me);
	double lambda[9];
	linear_combination_square_matrix(lambda, 
		body.alpha, body.b, -1.0 * body.alpha, b_eta);

	if (body.elements > 0)
	{
		double	bk_me_2d_array[body.elements][5];
		double	bk[body.elements][9];
		for (int k = 0; k < body.elements; k++)
		{
			for (int l = 0; l < 5; l++)
			{
				bk_me_2d_array[k][l] = body.bk_me[l + (k*5)];

			}
			construct_traceless_symmetric_matrix(bk[k], bk_me_2d_array[k]);

			// update lambda
			linear_combination_square_matrix(lambda, 
				1.0, lambda, -1.0 * body.alpha, bk[k]);
		}

		for (int k = 0; k < body.elements; k++)
		{
			double minus_bk_over_tau_elements[9];
			scale_square_matrix(minus_bk_over_tau_elements, 
				-1.0 * body.alpha_elements[k] / body.eta_elements[k],
				bk[k]);
			double lambda_over_eta_elements[9];
			scale_square_matrix(lambda_over_eta_elements,
				1.0 / body.eta_elements[k], lambda);

			/* Y * \dot{B}_{e_k} * Y^T */
			double Y_B_k_dot_T_trans[9];
			linear_combination_square_matrix(Y_B_k_dot_T_trans,
				1.0, minus_bk_over_tau_elements,
				1.0, lambda_over_eta_elements);

			/* ||\dot{B}_{e_k}||^2 = ||Y \dot{B}_{e_k} Y^T||^2 */
			double norm_squared_B_k_dot 
				= norm_squared_square_matrix(Y_B_k_dot_T_trans);

			D += body.I0 * (0.5 * body.eta_elements[k] * norm_squared_B_k_dot);
		}

	} // end if body.elements > 0

	/* Y * \dot{B}_eta * Y^T */
	double Y_B_eta_dot_T_trans[9];
	scale_square_matrix(Y_B_eta_dot_T_trans,
		1.0 / body.eta, lambda);

	/* ||\dot{B}_eta||^2 = ||Y \dot{B}_eta Y^T||^2 */
	double norm_squared_B_eta_dot 
		= norm_squared_square_matrix(Y_B_eta_dot_T_trans);

	D += body.I0 * (0.5 * body.eta * norm_squared_B_eta_dot);

	return D;
}
