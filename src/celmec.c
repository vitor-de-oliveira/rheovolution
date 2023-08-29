#include "celmec.h"

double
kepler_period(double m1, double m2, double G, double a)
{
	return 2.0 * M_PI * sqrt((a * a * a) / (G * (m1 + m2)));
}

double
kepler_period_only_m1(double m1, double G, double a)
{
	return 2.0 * M_PI * sqrt((a * a * a) / (G * m1));
}

double
root_function_kepler(double E, void *params)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	double M = p->M;
	return E - e * sin(E) - M;
}

double
root_derivative_kepler(double E, void *params)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	return 1.0 - e * cos(E);
}

void
root_fdf_kepler(double E, void *params, double *y, 
                    double *dy)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	double M = p->M;
	*y = E - e * sin(E) - M;
	*dy = 1.0 - e * cos(E);
}

double
kepler_equation(const double e, const double M)
{
	int status, iter = 0, max_iter = 100;
	double E0, E;
	struct root_params_kepler params_root = {e, M};

    /* initial guess (works up to e = 0.99) */
    E = M + e * sin(M);

    /* stonger initial guess */
    // E = M + e * sin(M) + 0.5 * e * e * sin(2.0 * M)
        // + (e * e * e / 8.0) * (3.0 * sin (3.0 * M) 
        // - sin(M));

    const gsl_root_fdfsolver_type *T 
        = gsl_root_fdfsolver_steffenson;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
    gsl_function_fdf FDF;
    FDF.f = &root_function_kepler;
    FDF.df = &root_derivative_kepler;
    FDF.fdf = &root_fdf_kepler;
    FDF.params = &params_root;
    gsl_root_fdfsolver_set(s, &FDF, E);

    do
    {
        iter++;
        gsl_root_fdfsolver_iterate(s);
        E0 = E;
        E = gsl_root_fdfsolver_root(s);
        status = gsl_root_test_delta(E, E0, 1e-15, 0);
        if (iter == max_iter)
        {
            printf("Warning: reached maximum iterate\n");
            return 1;
        }
    } while (status == GSL_CONTINUE);

    gsl_root_fdfsolver_free(s);

    return E;
}

double
calculate_eccentricity  (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[])
{
    // standard gravitational parameter
    double mu = G * (m1 + m2);

    // double x_norm = norm_vector(x);
    // double v_norm_squared = norm_squared_vector(v);
    // double x_dot_v = dot_product(x, v);

    // double term_1[] = {0.0, 0.0, 0.0};
    // scale_vector(term_1, v_norm_squared - (mu / x_norm), x);
    // double term_2[] = {0.0, 0.0, 0.0};
    // scale_vector(term_2, x_dot_v, v);

    // double e_vec[] = {0.0, 0.0, 0.0};
    // linear_combination_vector(e_vec, 
    //     1.0 / mu, term_1, -1.0 / mu, term_2);

    // orbital momentum vector
    double h[3];
    cross_product(h, x, v);

    double v_cross_h[3];
    cross_product(v_cross_h, v, h);
    double term_1[] = {0.0, 0.0, 0.0};
    scale_vector(term_1, 1.0 / mu, v_cross_h);
    double term_2[] = {0.0, 0.0, 0.0};
    scale_vector(term_2, -1.0 / norm_vector(x), x);

    double e_vec[] = {0.0, 0.0, 0.0};
    linear_combination_vector(e_vec, 
        1.0, term_1, 1.0, term_2);

    return norm_vector(e_vec);
}

double
calculate_semi_major_axis   (const double G,
                             const double m1,
                             const double m2,
                             const double x[],
                             const double v[])
{
    // standard gravitational parameter
    double mu = G * (m1 + m2);

    // double x_norm = norm_vector(x);
    // double v_norm_squared = norm_squared_vector(v);

    // // specific mechanical energy
    // double E = v_norm_squared / 2.0 - mu / x_norm;

    // return -1.0 * mu / (2.0 * E);

    return 1.0 / ((2.0 / norm_vector(x)) - (norm_squared_vector(v) / mu));
}
