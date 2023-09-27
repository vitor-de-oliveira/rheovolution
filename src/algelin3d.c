#include "algelin3d.h"

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
