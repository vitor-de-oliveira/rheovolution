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
calculate_semi_major_axis   (const double G,
                             const double m1,
                             const double m2,
                             const double x[],
                             const double v[])
{
    // standard gravitational parameter
    double mu = G * (m1 + m2);

    // semi-major axis
    double a = 1.0 / ((2.0 / norm_vector(x)) - (norm_squared_vector(v) / mu));

    return a;
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

    // orbital momentum vector
    double h[3];
    cross_product(h, x, v);

    // eccentricity vector
    double v_cross_h[3];
    cross_product(v_cross_h, v, h);
    double term_1[] = {0.0, 0.0, 0.0};
    scale_vector(term_1, 1.0 / mu, v_cross_h);
    double term_2[] = {0.0, 0.0, 0.0};
    scale_vector(term_2, -1.0 / norm_vector(x), x);
    double e_vec[] = {0.0, 0.0, 0.0};
    linear_combination_vector(e_vec, 
        1.0, term_1, 1.0, term_2);

    // eccentricity
    double e = norm_vector(e_vec);

    return e;
}

double
calculate_inclination   (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[])
{
    // orbital momentum vector
    double h[3];
    cross_product(h, x, v);

    // inclination
    double I = acos(h[2] / norm_vector(h));

    return I;
}

double
calculate_true_anomaly  (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[])
{
    // standard gravitational parameter
    double mu = G * (m1 + m2);

    // orbital momentum vector
    double h[3];
    cross_product(h, x, v);

    // eccentricity vector
    double v_cross_h[3];
    cross_product(v_cross_h, v, h);
    double term_1[] = {0.0, 0.0, 0.0};
    scale_vector(term_1, 1.0 / mu, v_cross_h);
    double term_2[] = {0.0, 0.0, 0.0};
    scale_vector(term_2, -1.0 / norm_vector(x), x);
    double e_vec[] = {0.0, 0.0, 0.0};
    linear_combination_vector(e_vec, 
        1.0, term_1, 1.0, term_2);

    // inclination
    double I = calculate_inclination(G, m1, m2, x, v);

    // eccentricity
    double e = calculate_eccentricity(G, m1, m2, x, v);

    // true anomaly
    double nu;

    if (e < 1e-15)
    {
        if (I < 1e-15)
        {
            nu = acos(x[0] / norm_vector(x));
            if (v[0] > 0.0)
            {
                nu = 2.0 * M_PI - nu;
            }
        }
        else
        {
            // vector pointing towards the ascending node
            double z_vec[] = {0.0, 0.0, 1.0};
            double n[3];
            cross_product(n, z_vec, h);

            nu = acos(dot_product(n, x) / (norm_vector(n) * norm_vector(x)));
            if (dot_product(n, v) > 0.0)
            {
                nu = 2.0 * M_PI - nu;
            }
        }
    }
    else
    {
        nu = acos(dot_product(e_vec, x) / (norm_vector(e_vec) *  norm_vector(x)));
        if (dot_product(x, v) < 0.0)
        {
            nu = 2.0 * M_PI - nu;
        } 
    }

    return nu;
}

double
calculate_eccentric_anomaly (const double G,
                             const double m1,
                             const double m2,
                             const double x[],
                             const double v[])
{
    // true anomaly
    double nu = calculate_true_anomaly(G, m1, m2, x, v);

    // eccentricity
    double e = calculate_eccentricity(G, m1, m2, x, v);

    // eccentric anomaly
    double E = 2.0 * atan2( tan(0.5*nu), sqrt((1.0 + e)/(1.0 - e)) );

    return E;
}

double
calculate_mean_anomaly  (const double G,
                         const double m1,
                         const double m2,
                         const double x[],
                         const double v[])
{
    // true anomaly
    double E = calculate_eccentric_anomaly(G, m1, m2, x, v);

    // eccentricity
    double e = calculate_eccentricity(G, m1, m2, x, v);

    // mean anomaly
    double M = E - e * sin(E);

    return M;
}


double
calculate_argument_of_periapsis (const double G,
                                 const double m1,
                                 const double m2,
                                 const double x[],
                                 const double v[])
{
    // standard gravitational parameter
    double mu = G * (m1 + m2);
    
    // orbital momentum vector
    double h[3];
    cross_product(h, x, v);

    // eccentricity vector
    double v_cross_h[3];
    cross_product(v_cross_h, v, h);
    double term_1[] = {0.0, 0.0, 0.0};
    scale_vector(term_1, 1.0 / mu, v_cross_h);
    double term_2[] = {0.0, 0.0, 0.0};
    scale_vector(term_2, -1.0 / norm_vector(x), x);
    double e_vec[] = {0.0, 0.0, 0.0};
    linear_combination_vector(e_vec, 
        1.0, term_1, 1.0, term_2);

    // inclination
    double I = calculate_inclination(G, m1, m2, x, v);

    // eccentricity
    double e = calculate_eccentricity(G, m1, m2, x, v);

    // argument of periapsis
    double w;

    if (e < 1e-15)
    {
        w = 0.0;
    }
    else
    {
        if (I < 1e-15)
        {
            w = acos(e_vec[0] / e);
        }
        else
        {
            // vector pointing towards the ascending node
            double z_vec[] = {0.0, 0.0, 1.0};
            double n[3];
            cross_product(n, z_vec, h);

            w = acos(dot_product(n, e_vec) / (norm_vector(n) * norm_vector(e_vec)));
        }
        if (e_vec[2] < 0.0)
        {
            w = 2.0 * M_PI - w;
        }
    }

    return w;
}

double
calculate_longitude_of_the_ascending_node   (const double G,
                                             const double m1,
                                             const double m2,
                                             const double x[],
                                             const double v[])
{
    // inclination
    double I = calculate_inclination(G, m1, m2, x, v);

    // longitude of the ascending node
    double Omega;

    if (I < 1e-15)
    {
        Omega = 0.0;
    }
    else
    {
        // orbital momentum vector
        double h[3];
        cross_product(h, x, v);

        // vector pointing towards the ascending node
        double z_vec[] = {0.0, 0.0, 1.0};
        double n[3];
        cross_product(n, z_vec, h);

        Omega = acos(n[0] / norm_vector(n));
        if (n[1] < 0.0)
        {
            Omega = 2.0 * M_PI - Omega;
        }
    }

    return Omega;
}

int
print_CelestialBody(cltbdy body)
{
	printf("name = ");
	printf("%s\n", body.name);

	printf("mass = ");
	printf("%1.5e\n", body.mass);
	printf("R = ");	
	printf("%1.5e\n", body.R);
	printf("lod = ");	
	printf("%1.5e\n", body.lod);

	printf("rg = ");	
	printf("%1.5e\n", body.rg);
	printf("J2 = ");	
	printf("%1.5e\n", body.J2);
	printf("C22 = ");	
	printf("%1.5e\n", body.C22);
	printf("I0 = ");	
	printf("%1.5e\n", body.I0);

	printf("obl = ");
	printf("%1.5e\n", body.obl);
	printf("psi = ");	
	printf("%1.5e\n", body.psi);
	printf("lib = ");	
	printf("%1.5e\n", body.lib);

	printf("a = ");
	printf("%1.5e\n", body.a);
	printf("e = ");	
	printf("%1.5e\n", body.e);
	printf("I = ");	
	printf("%1.5e\n", body.I);
	printf("M = ");	
	printf("%1.5e\n", body.M);
	printf("w = ");	
	printf("%1.5e\n", body.w);
	printf("Omega = ");	
	printf("%1.5e\n", body.Omega);

	printf("kf = ");	
	printf("%1.5e\n", body.kf);
	printf("Dt = ");	
	printf("%1.5e\n", body.Dt);
	printf("tau = ");	
	printf("%1.5e\n", body.tau);

	printf("gamma = ");	
	printf("%1.5e\n", body.gamma);
	printf("alpha = ");	
	printf("%1.5e\n", body.alpha);
	printf("eta = ");	
	printf("%1.5e\n", body.eta);
	printf("alpha_0 = ");	
	printf("%1.5e\n", body.alpha_0);
	printf("elements = ");	
	printf("%d\n", body.elements);
	for(int i = 0; i < body.elements; i++)
	{
		printf("alpha_%d = ", i+1);	
		printf("%1.5e\n", body.alpha_elements[i]);
		printf("eta_%d = ", i+1);	
		printf("%1.5e\n", body.eta_elements[i]);
	}

	printf("point mass = ");	
	printf("%d\n", body.point_mass);
	printf("centrifugal = ");	
	printf("%d\n", body.centrifugal);
	printf("tidal = ");	
	printf("%d\n", body.tidal);
	
	printf("x = ");
	print_vector(body.x);
	printf("x_dot = ");
	print_vector(body.x_dot);
	printf("l = ");
	print_vector(body.l);
	double b0[9];
	construct_traceless_symmetric_matrix(b0, body.b0_me);
	printf("b0 = \n");
	print_square_matrix(b0);
	double u[9];
	construct_traceless_symmetric_matrix(u, body.u_me);
	printf("u = \n");
	print_square_matrix(u);
	if (body.elements > 0)
	{
		double	bk_me_2d_array[body.elements][5];
		double	bk[9];
		for (int k = 0; k < body.elements; k++)
		{
			for (int l = 0; l < 5; l++)
			{
				bk_me_2d_array[k][l] = body.bk_me[l + (k*5)];
			}
			construct_traceless_symmetric_matrix(bk, bk_me_2d_array[k]);
			printf("b_%d = \n", k+1);	
			print_square_matrix(bk);
		}
	}

	printf("omega = ");
	print_vector(body.omega);
	printf("b = \n");
	print_square_matrix(body.b);

	return 0;
}

int
calculate_orbital_elements  (cltbdy *body, 
                             const cltbdy body_ref,
                             const double G)
{
	double relative_x[3];
	linear_combination_vector(relative_x,
		1.0, (*body).x,
		-1.0, body_ref.x);
	double relative_x_dot[3];
	linear_combination_vector(relative_x_dot,
		1.0, (*body).x_dot,
		-1.0, body_ref.x_dot);

    double m1, m2, x[3], v[3];

    m1 = body_ref.mass;
    m2 = (*body).mass;
    copy_vector(x, relative_x);
    copy_vector(v, relative_x_dot);

    // standard gravitational parameter
    double mu = G * (m1 + m2);

    // semi-major axis
    double a = 1.0 / ((2.0 / norm_vector(x)) - (norm_squared_vector(v) / mu));

    // orbital momentum vector
    double h[3];
    cross_product(h, x, v);

    // eccentricity vector
    double v_cross_h[3];
    cross_product(v_cross_h, v, h);
    double term_1[] = {0.0, 0.0, 0.0};
    scale_vector(term_1, 1.0 / mu, v_cross_h);
    double term_2[] = {0.0, 0.0, 0.0};
    scale_vector(term_2, -1.0 / norm_vector(x), x);
    double e_vec[] = {0.0, 0.0, 0.0};
    linear_combination_vector(e_vec, 
        1.0, term_1, 1.0, term_2);

    // node vector (vector pointing towards the ascending node)
    double z_vec[] = {0.0, 0.0, 1.0};
    double n[3];
    cross_product(n, z_vec, h);

    // eccentricity
    double e = norm_vector(e_vec);

    // inclination
    double I = acos(h[2] / norm_vector(h));

    // true anomaly
    double nu;

    if (e < 1e-15)
    {
        if (I < 1e-15)
        {
            nu = acos(x[0] / norm_vector(x));
            if (v[0] > 0.0)
            {
                nu = 2.0 * M_PI - nu;
            }
        }
        else
        {
            nu = acos(dot_product(n, x) / (norm_vector(n) * norm_vector(x)));
            if (dot_product(n, v) > 0.0)
            {
                nu = 2.0 * M_PI - nu;
            }
        }
    }
    else
    {
        nu = acos(dot_product(e_vec, x) / (norm_vector(e_vec) *  norm_vector(x)));
        if (dot_product(x, v) < 0.0)
        {
            nu = 2.0 * M_PI - nu;
        } 
    }

    // eccentric anomaly
    double E = 2.0 * atan2( tan(0.5*nu), sqrt((1.0 + e)/(1.0 - e)) );

    // mean anomaly
    double M = E - e * sin(E);

     // argument of periapsis
    double w;

    if (e < 1e-15)
    {
        w = 0.0;
    }
    else
    {
        if (I < 1e-15)
        {
            w = acos(e_vec[0] / e);
        }
        else
        {
            w = acos(dot_product(n, e_vec) / (norm_vector(n) * norm_vector(e_vec)));
        }
        if (e_vec[2] < 0.0)
        {
            w = 2.0 * M_PI - w;
        }
    }

    // longitude of the ascending node
    double Omega;

    if (I < 1e-15)
    {
        Omega = 0.0;
    }
    else
    {
        Omega = acos(n[0] / norm_vector(n));
        if (n[1] < 0.0)
        {
            Omega = 2.0 * M_PI - Omega;
        }
    }

    (*body).a = a;
    (*body).e = e;
    (*body).I = I;
    (*body).M = M;
    (*body).w = w;
    (*body).Omega = Omega;

    return 0;
}
