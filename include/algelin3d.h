#ifndef ALGELIN3D_H
#define ALGELIN3D_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIM_3 3 // makes it easier to update this lib later
                // to work on other dimensions

/* vector */

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
copy_square_matrix(double N[], const double M[]);

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
commutator(double O[], const double M[], const double N[]);

/* tensor */

int
tensor_product(double M[], const double x[], const double y[]);

#endif