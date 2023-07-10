/**
 * Library for simple linear algebra operations
 * It is currently implemented for 3d but can be easily extended
 * It is currently only implemented for double precision
 * 
 * Author: Vitor M. de Oliveira
 * Date: 12 june 2023
**/

#ifndef ALGELIN3D_H
#define ALGELIN3D_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

/* vector vector */

double
dot_product(const double x[], const double y[]);

int
cross_product(double z[], const double x[], const double y[]);

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

double
trace_square_matrix(const double M[]);

int
scale_square_matrix(double aM[], const double a, const double M[]);

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

/* tensor */

int
tensor_product(double M[], const double x[], const double y[]);

#endif