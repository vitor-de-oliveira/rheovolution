#include "celmec.h"

double
kepler_period(double m1, double m2, double G, double a)
{
	return sqrt(4.0 * M_PI * M_PI * a * a * a / (G * (m1 + m2)));
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

