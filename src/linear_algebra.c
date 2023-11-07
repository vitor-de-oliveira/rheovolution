#include "linear_algebra.h"

int
null_vector(double x[])
{
	for (int i = 0; i < DIM_3; i++)
	{
		x[i] = 0.0;
	}
	return 0;
}

int
print_vector(const double x[])
{
	for (int i = 0; i < DIM_3; i++)
	{
		printf("%1.10e ", x[i]);
	}
	printf("\n");
	return 0;
}

int
copy_vector(double y[], const double x[])
{
	for (int i = 0; i < DIM_3; i++)
	{
		y[i] = x[i];
	}
	return 0;
}

double
norm_squared_vector(const double x[])
{
	double result;
	double sum = 0.0;
	for (int i = 0; i < DIM_3; i++)
	{
		sum += x[i] * x[i];
	}
	result = sum;
	return result;
}

double
norm_vector(const double x[])
{
	return sqrt(norm_squared_vector(x));
}

int
scale_vector(double ax[], const double a, const double x[])
{
	double ax_local[DIM_3];

	for (int i = 0; i < DIM_3; i++)
	{
		ax_local[i] = a * x[i];
	}

	for (int i = 0; i < DIM_3; i++)
	{
		ax[i] = ax_local[i];
	}

	return 0;
}

int
vector_from_spherical_coordinates   (double x[],
                                     const double r,
                                     const double theta,
                                     const double phi)
{
	x[0] = r * cos(theta) * sin(phi);
	x[1] = r * sin(theta) * sin(phi);	
	x[2] = r * cos(phi);

	return 0;
}

double
dot_product(const double x[], const double y[])
{
	double result = 0.0;
	for (int i = 0; i < DIM_3; i++)
	{
		result += x[i] * y[i];
	}
	return result;
}

int
cross_product(double z[], const double x[], const double y[])
{
	double z_local[DIM_3];
	
	if (DIM_3 == 3)
	{
		z_local[0] = x[1] * y[2] - x[2] * y[1];
		z_local[1] = x[2] * y[0] - x[0] * y[2];
		z_local[2] = x[0] * y[1] - x[1] * y[0];
	}
	else
	{
		printf("Error: cross product for dimensions lower than 3");
		printf(" not implemented yet.\n");
		exit(4);
	}

	for (int i = 0; i < DIM_3; i++)
	{
		z[i] = z_local[i];
	}

	return 0;
}

double
angle_between_two_vectors(const double x[], const double y[])
{
	double omega;

	double frac = dot_product(x,y)/(norm_vector(x)*norm_vector(y));

	// tries to correct numerical errors
	if (frac > 1.0)
	{
		omega = 0.0;
	}
	else if (frac < -1.0)
	{
		omega = M_PI;
	}
	else
	{
		omega = acos(frac);
	}

	return omega;
}

int
linear_combination_vector(double z[], const double a, const double x[], 
    const double b, const double y[])
{
	double ax[DIM_3];
	scale_vector(ax, a, x);
	double by[DIM_3];
	scale_vector(by, b, y);

	for (int i = 0; i < DIM_3; i++)
	{
		z[i] = ax[i] + by[i];
	}

	return 0;
}

int
linear_combination_three_vector(double v[], const double a, const double x[], 
    const double b, const double y[], const double c, const double z[])
{
	double ax[DIM_3];
	scale_vector(ax, a, x);
	double by[DIM_3];
	scale_vector(by, b, y);
	double cz[DIM_3];
	scale_vector(cz, c, z);

	for (int i = 0; i < DIM_3; i++)
	{
		v[i] = ax[i] + by[i] + cz[i];
	}

	return 0;
}

int
null_matrix(double N[])
{
	for (int i = 0; i < DIM_3*DIM_3; i++)
	{
		N[i] = 0.0;
	}
	return 0;
}

int
identity_matrix(double I[])
{
	for (int i = 0; i < DIM_3*DIM_3; i++)
	{
		if (i % (DIM_3 + 1) == 0)
		{
			I[i] = 1.0;
		}
		else 
		{
			I[i] = 0.0;
		}
	}
	return 0;
}

int
rotation_matrix_3d_x(double R[], const double phi)
{
	R[0] = 1.0; R[1] = 0.0; 	 R[2] = 0.0;
	R[3] = 0.0; R[4] = cos(phi); R[5] = -1.0 * sin(phi);
	R[6] = 0.0; R[7] = sin(phi); R[8] = cos(phi);
	return 0;
}

int
rotation_matrix_3d_y(double R[], const double phi)
{
	R[0] = cos(phi); 		R[1] = 0.0; R[2] = sin(phi);
	R[3] = 0.0; 			R[4] = 1.0; R[5] = 0.0;
	R[6] = -1.0 * sin(phi); R[7] = 0.0; R[8] = cos(phi);
	return 0;
}

int
rotation_matrix_3d_z(double R[], const double phi)
{
	R[0] = cos(phi); R[1] = -1.0 * sin(phi); R[2] = 0.0;
	R[3] = sin(phi); R[4] = cos(phi); 		 R[5] = 0.0;
	R[6] = 0.0; 	 R[7] = 0.0; 			 R[8] = 1.0;
	return 0;
}

int
print_square_matrix(const double M[])
{
	for (int i = 0; i < DIM_3 * DIM_3; i++)
	{
		printf("%1.10e ", M[i]);
		if ((i+1) % DIM_3 == 0) printf("\n");
	}
	return 0;
}

int
copy_square_matrix(double N[], const double M[])
{
	for (int i = 0; i < DIM_3*DIM_3; i++)
	{
		N[i] = M[i];
	}
	return 0;
}

int
transpose_square_matrix(double Mt[], const double M[])
{
	for (int i = 0; i < DIM_3; i++)
	{
		for (int j = 0; j < DIM_3; j++)
		{
			Mt[i*DIM_3 + j] = M[i + j*DIM_3];
		}
	}
	return 0;
}

double
trace_square_matrix(const double M[])
{
	double trace = 0.0;
	for (int i = 0; i < DIM_3*DIM_3; i = i + (DIM_3 + 1))
	{
		trace += M[i];
	}
	return trace;
}

double
norm_squared_square_matrix(const double M[])
{
	return inner_product_square_matrix(M, M);
}

double
norm_square_matrix(const double M[])
{
	return sqrt(norm_squared_square_matrix(M));
}

int
scale_square_matrix(double aM[], const double a, const double M[])
{
	double aM_local[DIM_3*DIM_3];

	for (int i = 0; i < DIM_3; i++)
	{
		for (int j = 0; j < DIM_3; j++)
		{
			aM_local[(DIM_3)*i + j] = a * M[(DIM_3)*i + j];
		}
	}

	for (int i = 0; i < DIM_3 * DIM_3; i++)
	{
		aM[i] = aM_local[i];
	}

	return 0;
}

int
calculate_eigenvectors_matrix(double M_eig[], const double M[])
{

	double M_copy[9];
	copy_square_matrix(M_copy, M);

	gsl_matrix_view M_gsl
		= gsl_matrix_view_array (M_copy, 3, 3);

	gsl_vector *eval = gsl_vector_alloc (3);
	gsl_matrix *evec = gsl_matrix_alloc (3, 3);

	gsl_eigen_symmv_workspace * w =	gsl_eigen_symmv_alloc (3);

	gsl_eigen_symmv (&M_gsl.matrix, eval, evec, w);

	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

	// choose the direction of the last eigenvector
	// as having a positive last entry
	for (int i = 0; i < 9; i++)
	{
		M_eig[i] = evec->data[i];
	}
	if (M_eig[8] < 0.0)
	{
		M_eig[2] *= -1.0;
		M_eig[5] *= -1.0;
		M_eig[8] *= -1.0;
	}

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_eigen_symmv_free(w);

	return 0;
}

int
calculate_square_matrix_inverse(double M_inv[], const double M[])
{
	double M_copy[9];
	copy_square_matrix(M_copy, M);

	gsl_matrix_view M_gsl
		= gsl_matrix_view_array (M_copy, 3, 3);

	int s;
	gsl_permutation * p = gsl_permutation_alloc (3);
	gsl_linalg_LU_decomp (&M_gsl.matrix, p, &s);

	gsl_matrix *M_inv_gsl = gsl_matrix_alloc (3, 3);

	gsl_linalg_LU_invert(&M_gsl.matrix, p, M_inv_gsl);

	for (int i = 0; i < 9; i++)
	{
		M_inv[i] = M_inv_gsl->data[i];
	}

	gsl_permutation_free (p);
	gsl_matrix_free (M_inv_gsl);

	return 0;
}

int
calculate_diagonalized_square_matrix(double M_diag[], const double M[])
{
	/* M_diag = P_inv * M * P */
	
	double P[9];
	calculate_eigenvectors_matrix(P, M);

	double P_inv[9];
	calculate_square_matrix_inverse(P_inv, P);

	square_matrix_times_square_matrix(M_diag, P_inv, M);
	square_matrix_times_square_matrix(M_diag, M_diag, P);

	return 0;
}

int
square_matrix_times_vector(double y[], const double M[], const double x[])
{
	double y_local[DIM_3];
	double result;

	for (int i = 0; i < DIM_3; i++)
	{
		result = 0.0;
		for (int j = 0; j < DIM_3; j++)
		{
			result += M[(DIM_3)*i + j] * x[j];
		}
		y_local[i] = result;
	}

	for (int i = 0; i < DIM_3; i++)
	{
		y[i] = y_local[i];
	}

	return 0;
}

int
linear_combination_vector_square_matrix(double z[], const double a, 
    const double M[], const double x[], const double b, const double y[])
{
	double Mx[DIM_3];
	square_matrix_times_vector(Mx, M, x);
	double aMx[DIM_3];
	scale_vector(aMx, a, Mx);
	double by[DIM_3];
	scale_vector(by, b, y);

	for (int i = 0; i < DIM_3; i++)
	{
		z[i] = aMx[i] + by[i];
	}

	return 0;
}

int
square_matrix_times_square_matrix(double MN[], 
    const double M[], const double N[])
{
	double MN_local[DIM_3*DIM_3];
	double result;

	for (int i = 0; i < DIM_3; i++)
	{
		for (int j = 0; j < DIM_3; j++)
		{
			result = 0.0;
			for (int k = 0; k < DIM_3; k++)
			{
				result += M[DIM_3*i + k] * N[DIM_3*k + j];
			}
			MN_local[i*DIM_3 + j] = result;
		}
	}

	for (int i = 0; i < DIM_3 * DIM_3; i++)
	{
		MN[i] = MN_local[i];
	}

	return 0;
}

int
linear_combination_square_matrix(double O[], const double a, const double M[], 
    const double b, const double N[])
{
	double aM[DIM_3*DIM_3];
	scale_square_matrix(aM, a, M);
	double bN[DIM_3*DIM_3];
	scale_square_matrix(bN, b, N);

	for (int i = 0; i < DIM_3 * DIM_3; i++)
	{
		O[i] = aM[i] + bN[i];
	}
	
	return 0;
}

int
linear_combination_three_square_matrix(double P[], const double a, const double M[], 
    const double b, const double N[], const double c, const double O[])
{
	double aM[DIM_3*DIM_3];
	scale_square_matrix(aM, a, M);
	double bN[DIM_3*DIM_3];
	scale_square_matrix(bN, b, N);
	double cO[DIM_3*DIM_3];
	scale_square_matrix(cO, c, O);

	for (int i = 0; i < DIM_3 * DIM_3; i++)
	{
		P[i] = aM[i] + bN[i] + cO[i];
	}
	
	return 0;
}

int
commutator(double O[], const double M[], const double N[])
{
	double MN[DIM_3*DIM_3];
	square_matrix_times_square_matrix(MN, M, N);
	double NM[DIM_3*DIM_3];
	square_matrix_times_square_matrix(NM, N, M);

	for (int i = 0; i < DIM_3 * DIM_3; i++)
	{
		O[i] = MN[i] - NM[i];
	}

	return 0;
}

double
inner_product_square_matrix(const double M[], const double N[])
{
	double N_transpose[9];
	transpose_square_matrix(N_transpose, N);
	double M_times_N_transpose[9];
	square_matrix_times_square_matrix(M_times_N_transpose,
		M, N_transpose);

	return 0.5 * trace_square_matrix(M_times_N_transpose);
}

int
tensor_product(double M[], const double x[], const double y[])
{
	for (int i = 0; i < DIM_3; i++)
	{
		for (int j = 0; j < DIM_3; j++)
		{
			M[(DIM_3)*i + j] = x[i] * y[j];
		}
	}

	return 0;
}

int
print_quaternion(const double q[])
{
	for (int i = 0; i < 4; i++)
	{
		printf("%1.10e ", q[i]);
	}
	printf("\n");
	return 0;
}

int
copy_quaternion(double qc[], const double q[])
{
	for (int i = 0; i < 4; i++)
	{
		qc[i] = q[i];
	}
	return 0;
}

double
norm_quaternion(const double q[])
{
	return sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
}

int
normalize_quaternion(double q_to_normalize[])
{
	double q_norm = norm_quaternion(q_to_normalize);

	if (q_norm > 1.0e-14)	// q not null
	{
		for (int i = 0; i < 4; i++)
		{
			q_to_normalize[i] /= q_norm;
		}
	}

	return 0;
}

int
identity_quaternion(double qI[])
{
	qI[0] = 1.0;
	qI[1] = 0.0;
	qI[2] = 0.0;
	qI[3] = 0.0;
	return 0;
}

int
quaternion_from_vector(double qv[4], const double v[3])
{
	qv[0] = 0.0;
	qv[1] = v[0];
	qv[2] = v[1];
	qv[3] = v[2];
	return 0;
}

int
conjugate_quaternion(double qc[], double q[])
{
	qc[0] = q[0];
	qc[1] = -1.0 * q[1];
	qc[2] = -1.0 * q[2];
	qc[3] = -1.0 * q[3];
	return 0;
}

int
quaternion_times_quaternion(double t[], const double r[], const double s[])
{
	// I am defining everything by hand for now
	// would have to extend this lib to 4d otherwise
	double t_local[4];

	t_local[0] = r[0]*s[0] - r[1]*s[1] - r[2]*s[2] - r[3]*s[3];
	t_local[1] = r[1]*s[0] + r[0]*s[1] - r[3]*s[2] + r[2]*s[3];
	t_local[2] = r[2]*s[0] + r[3]*s[1] + r[0]*s[2] - r[1]*s[3];
	t_local[3] = r[3]*s[0] - r[2]*s[1] + r[1]*s[2] + r[0]*s[3];  

	for (int i = 0; i < 4; i++)
	{
		t[i] = t_local[i];
	}

	return 0;
}

int
rotation_quaternion(double qr[4], const double alpha, const double u[3])
{
	double scaled_u[3];
	scale_vector(scaled_u, sin(0.5 * alpha), u);

	qr[0] = cos(0.5 * alpha);
	qr[1] = scaled_u[0];
	qr[2] = scaled_u[1];
	qr[3] = scaled_u[2];

	return 0;
}

int
rotation_quaternion_x(double qr[4], const double alpha)
{
	double u[] = {1.0, 0.0, 0.0};

	rotation_quaternion(qr, alpha, u);

	return 0;
}

int
rotation_quaternion_y(double qr[4], const double alpha)
{
	double u[] = {0.0, 1.0, 0.0};

	rotation_quaternion(qr, alpha, u);

	return 0;
}

int
rotation_quaternion_z(double qr[4], const double alpha)
{
	double u[] = {0.0, 0.0, 1.0};

	rotation_quaternion(qr, alpha, u);

	return 0;
}

int
rotate_vector_with_quaternion(double v_rot[3], const double q[4], const double v[3])
{
	// I am defining everything by hand for now
	// this should be the same as [0,v_rot] = q x [0,v] x q_conjugated
	double M[9];

	M[0] = q[0]*q[0] + q[1]*q[1] - 0.5;
	M[1] = q[1]*q[2] - q[0]*q[3];	
	M[2] = q[1]*q[3] + q[0]*q[2];
	M[3] = q[1]*q[2] + q[0]*q[3];
	M[4] = q[0]*q[0] + q[2]*q[2] - 0.5;
	M[5] = q[2]*q[3] - q[0]*q[1];
	M[6] = q[1]*q[3] - q[0]*q[2];
	M[7] = q[2]*q[3] + q[0]*q[1];
	M[8] = q[0]*q[0] + q[3]*q[3] - 0.5;

	double M_times_v[3];
	square_matrix_times_vector(M_times_v, M, v);

	scale_vector(v_rot, 2.0, M_times_v);

	return 0;
}

int
rotation_matrix_from_quaternion(double R[9], const double q[4])
{
	// I am defining everything by hand for now
	// should be a way to get here from standard operations

	// R[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
	// R[1] = 2.0 * (q[1]*q[2] - q[0]*q[3]);
	// R[2] = 2.0 * (q[1]*q[3] + q[0]*q[2]);
	// R[3] = 2.0 * (q[1]*q[2] + q[0]*q[3]);
	// R[4] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
	// R[5] = 2.0 * (q[2]*q[3] - q[0]*q[1]);
	// R[6] = 2.0 * (q[1]*q[3] - q[0]*q[2]);
	// R[7] = 2.0 * (q[2]*q[3] + q[0]*q[1]);
	// R[8] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];

	// from wikipedia
	double a = q[0];
	double b = q[1];
	double c = q[2];
	double d = q[3];

	double s = 2.0 / (a*a + b*b + c*c + d*d);
	double bs = b * s;
	double cs = c * s;
	double ds = d * s;
	double ab = a * bs;
	double ac = a * cs;
	double ad = a * ds;
	double bb = b * bs;
	double bc = b * cs;
	double bd = b * ds;
	double cc = c * cs;
	double cd = c * ds;
	double dd = d * ds;

	R[0] = 1.0 - cc - dd;
	R[1] = bc - ad;
	R[2] = bd + ac;
	R[3] = bc + ad;
	R[4] = 1.0 - bb - dd;
	R[5] = cd - ab;
	R[6] = bd - ac;
	R[7] = cd + ab;
	R[8] = 1.0 - bb - cc;

	return 0;
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
