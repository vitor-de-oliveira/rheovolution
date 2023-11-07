/**
 * Library for simple linear algebra operations
 * It is currently implemented for 3d but can be easily extended
 * It is currently only implemented for double precision
 * 
 * Author: Vitor M. de Oliveira
 * Date: 12 june 2023
**/

#ifndef LIN_ALG_H
#define LIN_ALG_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>

#define DIM_3 3 // makes it easier to update this lib later
                // to work on other dimensions

/* vector */

int
null_vector(double x[]);

int
print_vector(const double x[]);

int
copy_vector(double y[], const double x[]);

double
norm_squared_vector(const double x[]);

double
norm_vector(const double x[]);

int
scale_vector(double ax[], const double a, const double x[]);

int
vector_from_spherical_coordinates   (double x[],
                                     const double r,
                                     const double theta,
                                     const double phi);

/* vector vector */

double
dot_product(const double x[], const double y[]);

int
cross_product(double z[], const double x[], const double y[]);

double
angle_between_two_vectors(const double x[], const double y[]);

int
linear_combination_vector(double z[], const double a, const double x[], 
    const double b, const double y[]);

int
linear_combination_three_vector(double v[], const double a, const double x[], 
    const double b, const double y[], const double c, const double z[]);

/* square matrix */

int
null_matrix(double N[]);

int
identity_matrix(double I[]);

int
rotation_matrix_3d_x(double R[], const double phi);

int
rotation_matrix_3d_y(double R[], const double phi);

int
rotation_matrix_3d_z(double R[], const double phi);

int
print_square_matrix(const double M[]);

int
copy_square_matrix(double N[], const double M[]);

int
transpose_square_matrix(double Mt[], const double M[]);

double
trace_square_matrix(const double M[]);

double
norm_squared_square_matrix(const double M[]);

double
norm_square_matrix(const double M[]);

int
scale_square_matrix(double aM[], const double a, const double M[]);

// only defined for real symmetric matrices
int
calculate_eigenvectors_matrix(double M_eig[], const double M[]);

// only defined for real symmetric matrices
int
calculate_square_matrix_inverse(double M_inv[], const double M[]);

// only defined for real symmetric matrices
int
calculate_diagonalized_square_matrix(double M_diag[], const double M[]);

/* square matrix vector */

int
square_matrix_times_vector(double y[], const double M[], const double x[]);

int
linear_combination_vector_square_matrix(double z[], const double a, 
    const double M[], const double x[], const double b, const double y[]);

/* square matrix square matrix */

int
square_matrix_times_square_matrix(double MN[], 
    const double M[], const double N[]);

int
linear_combination_square_matrix(double O[], const double a, const double M[], 
    const double b, const double N[]);

int
linear_combination_three_square_matrix(double P[], const double a, const double M[], 
    const double b, const double N[], const double c, const double O[]);
    
int
commutator(double O[], const double M[], const double N[]);

double
inner_product_square_matrix(const double M[], const double N[]);

/* tensor */

int
tensor_product(double M[], const double x[], const double y[]);

/* quaternion */

int
print_quaternion(const double q[]);

int
copy_quaternion(double qc[], const double q[]);

double
norm_quaternion(const double q[]);

int
normalize_quaternion(double q_to_normalize[]);

int
identity_quaternion(double qI[]);

int
quaternion_from_vector(double qv[4], const double v[3]);

int
conjugate_quaternion(double qc[], double q[]);

int
quaternion_times_quaternion(double t[], const double r[], const double s[]);

int
rotation_quaternion(double qr[4], const double alpha, const double u[3]);

int
rotation_quaternion_x(double qr[4], const double alpha);

int
rotation_quaternion_y(double qr[4], const double alpha);

int
rotation_quaternion_z(double qr[4], const double alpha);

int
rotate_vector_with_quaternion(double v_rot[3], const double q[4], const double v[3]);

int
rotation_matrix_from_quaternion(double R[9], const double q[4]);

/* additional functions */

int
hat_map(double x_hat[9], const double x[3]);

int
construct_traceless_symmetric_matrix(double M[9], 
	const double M_main_elements[5]);

int
get_main_elements_traceless_symmetric_matrix(double M_main_elements[5], 
	const double M[9]);

#endif
