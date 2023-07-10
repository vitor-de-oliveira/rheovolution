#include "celmec.h"

double
kepler_period(double m1, double m2, double G, double a)
{
	return sqrt(4.0 * M_PI * M_PI * a * a * a / (G * (m1 + m2)));
}

double
root_function_kepler(double u, void *params)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	double t = p->t;
	return u - e * sin(u) - t;
}

double
root_derivative_kepler(double u, void *params)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	return 1.0 - e * cos(u);
}

void
root_fdf_kepler(double u, void *params, double *y, 
                    double *dy)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	double t = p->t;
	*y = u - e * sin(u) - t;
	*dy = 1.0 - e * cos(u);
}

double
kepler_equation(const double e, const double t)
{
	int status, iter = 0, max_iter = 100;
	double u0, u;
	struct root_params_kepler params_root = {e, t};

    /* initial guess (works up to e = 0.99) */
    u = t + e * sin(t);

    /* stonger initial guess */
    // u = t + e * sin(t) + 0.5 * e * e * sin(2.0 * t)
        // + (e * e * e / 8.0) * (3.0 * sin (3.0 * t) 
        // - sin(t));

    const gsl_root_fdfsolver_type *T 
        = gsl_root_fdfsolver_steffenson;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
    gsl_function_fdf FDF;
    FDF.f = &root_function_kepler;
    FDF.df = &root_derivative_kepler;
    FDF.fdf = &root_fdf_kepler;
    FDF.params = &params_root;
    gsl_root_fdfsolver_set(s, &FDF, u);

    do
    {
        iter++;
        gsl_root_fdfsolver_iterate(s);
        u0 = u;
        u = gsl_root_fdfsolver_root(s);
        status = gsl_root_test_delta(u, u0, 1e-15, 0);
        if (iter == max_iter)
        {
            printf("Warning: reached maximum iterate\n");
            return 1;
        }
    } while (status == GSL_CONTINUE);

    gsl_root_fdfsolver_free(s);

    return u;
}

