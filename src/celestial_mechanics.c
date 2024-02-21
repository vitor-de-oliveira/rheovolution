#include "celestial_mechanics.h"

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
calculate_rg(const double m, const double R, const double I[9])
{
	double rg;

	double I_33 = I[8];

	rg = I_33 / (m * R * R);

	return rg;
}


double
calculate_J2(const double m, const double R, const double I[9])
{
	double J2;

	double I_11 = I[0];
	double I_22 = I[4];
	double I_33 = I[8];

	J2 = (2.0 * I_33 - I_11 - I_22) / (2.0 * m * R * R);

	return J2;
}

double
calculate_C22(const double m, const double R, const double I[9])
{
	double C22;

	double I_11 = I[0];
	double I_22 = I[4];

	C22 = (I_22 - I_11) / (4.0 * m * R * R);

	return C22;
}

double
calculate_S22(const double m, const double R, const double I[9])
{
	double S22;

	double I_12 = I[1];

	S22 = -1.0 * I_12 / (2.0 * m * R * R);

	return S22;
}

double
calculate_C21(const double m, const double R, const double I[9])
{
	double C21;

	double I_13 = I[2];

	C21 = -1.0 * I_13 / (m * R * R);

	return C21;
}

double
calculate_S21(const double m, const double R, const double I[9])
{
	double S21;

	double I_23 = I[5];

	S21 = -1.0 * I_23 / (m * R * R);

	return S21;
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
copy_CelestialBody	(cltbdy *body_dest,
					 const cltbdy body_src)
{
    if (body_dest->elements != body_src.elements)
    {
        fprintf(stderr, "Error: cltbdy structs with ");
        fprintf(stderr, "different sizes being copied.\n");
        return -1;
    }

    strcpy(body_dest->name, body_src.name);
    
    body_dest->mass = body_src.mass;
    body_dest->R = body_src.R;
    body_dest->lod = body_src.lod;
    body_dest->azi = body_src.azi;
    body_dest->pol = body_src.pol;

    body_dest->I0 = body_src.I0;
    body_dest->rg = body_src.rg;
    body_dest->J2 = body_src.J2;
    body_dest->C22 = body_src.C22;
    body_dest->S22 = body_src.S22;
    body_dest->C21 = body_src.C21;
    body_dest->S21 = body_src.S21;

    body_dest->obl = body_src.obl;
    body_dest->psi = body_src.psi;
    body_dest->lib = body_src.lib;

    body_dest->a = body_src.a;
    body_dest->e = body_src.e;
    body_dest->I = body_src.I;
    body_dest->M = body_src.M;
    body_dest->w = body_src.w;
    body_dest->Omega = body_src.Omega;

    body_dest->k0 = body_src.k0;
    body_dest->Dt = body_src.Dt;
    body_dest->tau = body_src.tau;

    body_dest->gamma_0 = body_src.gamma_0;
    body_dest->alpha = body_src.alpha;
    body_dest->eta = body_src.eta;
    for (int i = 0; i < body_dest->elements; i++)
    {
        body_dest->alpha_elements[i] = body_src.alpha_elements[i];
        body_dest->eta_elements[i] = body_src.eta_elements[i];
    }

    body_dest->point_mass = body_src.point_mass;
    body_dest->deformable = body_src.deformable;
    body_dest->prestress = body_src.prestress;
    body_dest->centrifugal = body_src.centrifugal;
    body_dest->tidal = body_src.tidal;

    for (int i = 0; i < 3; i++)
    {
        body_dest->x[i] = body_src.x[i];
        body_dest->x_dot[i] = body_src.x_dot[i];
        body_dest->l[i] = body_src.l[i];
    }
    for (int i = 0; i < 5; i++)
    {
        body_dest->p_me[i] = body_src.p_me[i];
        body_dest->b_eta_me[i] = body_src.b_eta_me[i];
    }
    for (int i = 0; i < body_dest->elements * 5; i++)
    {
        body_dest->bk_me[i] = body_src.bk_me[i];
    }

    for (int i = 0; i < 5; i++)
    {
        body_dest->bs_me[i] = body_src.bs_me[i];
    }

    for (int i = 0; i < 4; i++)
    {
        body_dest->q[i] = body_src.q[i];
    }

    for (int i = 0; i < 3; i++)
    {
        body_dest->omega[i] = body_src.omega[i];
    }
    for (int i = 0; i < 9; i++)
    {
        body_dest->b[i] = body_src.b[i];
    }

    for (int i = 0; i < 3; i++)
    {
        body_dest->relative_x[i] = body_src.relative_x[i];
        body_dest->relative_x_dot[i] = body_src.relative_x_dot[i];
    }

    return 0;
}

cltbdy
create_and_copy_CelestialBody(const cltbdy body_src)
{
    cltbdy body_copy;

    body_copy.elements = body_src.elements;

    if (body_copy.elements > 0)
	{
		body_copy.alpha_elements 
			= (double *) malloc(body_copy.elements * sizeof(double));
		body_copy.eta_elements 	
			= (double *) malloc(body_copy.elements * sizeof(double));
        body_copy.bk_me 	
			= (double *) malloc(body_copy.elements * 5 * sizeof(double));
    }

    copy_CelestialBody(&body_copy, body_src);

    return body_copy;
}
                    

int
print_CelestialBody(cltbdy body)
{
	printf("name = ");
	printf("%s\n", body.name);

	printf("mass = ");
	printf("%1.10e\n", body.mass);
	printf("R = ");	
	printf("%1.10e\n", body.R);
	printf("lod = ");	
	printf("%1.10e\n", body.lod);

	printf("I0 = ");	
	printf("%1.10e\n", body.I0);
	printf("rg = ");	
	printf("%1.10e\n", body.rg);
	printf("J2 = ");	
	printf("%1.10e\n", body.J2);
	printf("C22 = ");	
	printf("%1.10e\n", body.C22);
	printf("S22 = ");	
	printf("%1.10e\n", body.S22);
	printf("C21 = ");	
	printf("%1.10e\n", body.C21);
	printf("S21 = ");	
	printf("%1.10e\n", body.S21);

	printf("obl = ");
	printf("%1.10e\n", body.obl);
	printf("psi = ");	
	printf("%1.10e\n", body.psi);
	printf("lib = ");	
	printf("%1.10e\n", body.lib);

	printf("a = ");
	printf("%1.10e\n", body.a);
	printf("e = ");	
	printf("%1.10e\n", body.e);
	printf("I = ");	
	printf("%1.10e\n", body.I);
	printf("M = ");	
	printf("%1.10e\n", body.M);
	printf("w = ");	
	printf("%1.10e\n", body.w);
	printf("Omega = ");	
	printf("%1.10e\n", body.Omega);

	printf("k0 = ");	
	printf("%1.10e\n", body.k0);
	printf("Dt = ");	
	printf("%1.10e\n", body.Dt);
	printf("tau = ");	
	printf("%1.10e\n", body.tau);

	printf("gamma_0 = ");	
	printf("%1.10e\n", body.gamma_0);
	printf("alpha = ");	
	printf("%1.10e\n", body.alpha);
	printf("eta = ");	
	printf("%1.10e\n", body.eta);
	printf("elements = ");	
	printf("%d\n", body.elements);
	for(int i = 0; i < body.elements; i++)
	{
		printf("alpha_%d = ", i+1);	
		printf("%1.10e\n", body.alpha_elements[i]);
		printf("eta_%d = ", i+1);	
		printf("%1.10e\n", body.eta_elements[i]);
	}

	printf("point mass = ");	
	printf("%d\n", body.point_mass);
	printf("deformable = ");	
	printf("%d\n", body.deformable);
	printf("prestress = ");	
	printf("%d\n", body.prestress);
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
	double p[9];
	construct_traceless_symmetric_matrix(p, body.p_me);
	printf("p = \n");
	print_square_matrix(p);
	double b_eta[9];
	construct_traceless_symmetric_matrix(b_eta, body.b_eta_me);
	printf("b_eta = \n");
	print_square_matrix(b_eta);
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
	
    double bs[9];
	construct_traceless_symmetric_matrix(bs, body.bs_me);
	printf("bs = \n");
	print_square_matrix(bs);

	printf("q = \n");
	print_quaternion(body.q);

	printf("omega = ");
	print_vector(body.omega);
	printf("b = \n");
	print_square_matrix(body.b);

    printf("relative_x = ");
	print_vector(body.relative_x);
	printf("relative_x_dot = ");
	print_vector(body.relative_x_dot);

	return 0;
}

int
print_state_variables(cltbdy body)
{
	printf("x = ");
	print_vector(body.x);
	printf("x_dot = ");
	print_vector(body.x_dot);
	printf("l = ");
	print_vector(body.l);
	double p[9];
	construct_traceless_symmetric_matrix(p, body.p_me);
	printf("p = \n");
	print_square_matrix(p);
	double b_eta[9];
	construct_traceless_symmetric_matrix(b_eta, body.b_eta_me);
	printf("b_eta = \n");
	print_square_matrix(b_eta);
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
	double bs[9];
	construct_traceless_symmetric_matrix(bs, body.bs_me);
	printf("bs = \n");
	print_square_matrix(bs);

	return 0;
}

int
calculate_orbital_elements  (cltbdy *body, 
                             const cltbdy body_ref,
                             const double G)
{
    double a, e, I, M, w, Omega; // orbital elements

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

    /* separate functions */
    // // semi-major axis
    // a = calculate_semi_major_axis(G, m1, m2, x, v);

    // // eccentricity
    // e = calculate_eccentricity(G, m1, m2, x, v);

    // // inclination
    // I = calculate_inclination(G, m1, m2, x, v);

    // // mean anomaly
    // M = calculate_mean_anomaly(G, m1, m2, x, v);

    // // argument of periapsis
    // w = calculate_argument_of_periapsis(G, m1, m2, x, v);

    // // longitude of the ascending node
    // Omega = calculate_longitude_of_the_ascending_node(G, m1, m2, x, v);

    /* inline */
    // // standard gravitational parameter
    // double mu = G * (m1 + m2);

    // // orbital momentum vector
    // double h[3];
    // cross_product(h, x, v);

    // // eccentricity vector
    // double v_cross_h[3];
    // cross_product(v_cross_h, v, h);
    // double term_1[] = {0.0, 0.0, 0.0};
    // scale_vector(term_1, 1.0 / mu, v_cross_h);
    // double term_2[] = {0.0, 0.0, 0.0};
    // scale_vector(term_2, -1.0 / norm_vector(x), x);
    // double e_vec[] = {0.0, 0.0, 0.0};
    // linear_combination_vector(e_vec, 
    //     1.0, term_1, 1.0, term_2);

    // // node vector (vector pointing towards the ascending node)
    // double z_vec[] = {0.0, 0.0, 1.0};
    // double n[3];
    // cross_product(n, z_vec, h);

    // // semi-major axis
    // a = 1.0 / ((2.0 / norm_vector(x)) - (norm_squared_vector(v) / mu));

    // // eccentricity
    // e = norm_vector(e_vec);

    // // inclination
    // I = acos(h[2] / norm_vector(h));

    // // true anomaly
    // double nu;

    // if (e < 1e-15)
    // {
    //     if (I < 1e-15)
    //     {
    //         nu = acos(x[0] / norm_vector(x));
    //         if (v[0] > 0.0)
    //         {
    //             nu = 2.0 * M_PI - nu;
    //         }
    //     }
    //     else
    //     {
    //         nu = acos(dot_product(n, x) / (norm_vector(n) * norm_vector(x)));
    //         if (dot_product(n, v) > 0.0)
    //         {
    //             nu = 2.0 * M_PI - nu;
    //         }
    //     }
    // }
    // else
    // {
    //     nu = acos(dot_product(e_vec, x) / (e *  norm_vector(x)));
    //     if (dot_product(x, v) < 0.0)
    //     {
    //         nu = 2.0 * M_PI - nu;
    //     } 
    // }

    // // eccentric anomaly
    // double E = 2.0 * atan2( tan(0.5*nu), sqrt((1.0 + e)/(1.0 - e)) );

    // // mean anomaly
    // M = E - e * sin(E);

    //  // argument of periapsis
    // if (e < 1e-15)
    // {
    //     w = 0.0;
    // }
    // else
    // {
    //     if (I < 1e-15)
    //     {
    //         w = acos(e_vec[0] / e);
    //     }
    //     else
    //     {
    //         w = acos(dot_product(n, e_vec) / (norm_vector(n) * e));
    //     }
    //     if (e_vec[2] < 0.0)
    //     {
    //         w = 2.0 * M_PI - w;
    //     }
    // }

    // // longitude of the ascending node
    // if (I < 1e-15)
    // {
    //     Omega = 0.0;
    // }
    // else
    // {
    //     Omega = acos(n[0] / norm_vector(n));
    //     if (n[1] < 0.0)
    //     {
    //         Omega = 2.0 * M_PI - Omega;
    //     }
    // }

    /* Murray */

    // // semi-major axis
    // a = 1.0 / ((2.0 / norm_vector(x)) - (norm_squared_vector(v) / mu));

    // // eccentricity
    // e = sqrt(1.0 - norm_squared_vector(h) / (mu *a));
    
    // // inclination
    // I = acos(h[2]/norm_vector(h));

    // // longitude of the ascending node
    // double ac_Omega;
    // if (h[2] > 0.0)
    // {
    //     Omega = asin(h[0]/(norm_vector(h)*sin(I)));
    //     ac_Omega = acos(-h[1]/(norm_vector(h)*sin(I)));    
    // }
    // else
    // {
    //     Omega = asin(-h[0]/(norm_vector(h)*sin(I)));
    //     ac_Omega = acos(h[1]/(norm_vector(h)*sin(I)));
    // }
    // if (ac_Omega > M_PI / 2.0)
    // {
    //     if (Omega > 0.0)
    //     {
    //         Omega += 2.0 * (M_PI/2.0 - Omega);
    //     }
    //     else
    //     {
    //         Omega -= 2.0 * (-M_PI/2.0 - Omega);
    //     }
    // }

    // // true anomaly
    // double R_dot = sqrt(norm_squared_vector(v)-norm_squared_vector(h)/norm_squared_vector(x));
    // double x_dot_v = dot_product(x, v);
    // if (x_dot_v) R_dot *= -1.0;
    // double f, ac_f;
    // f = asin((a*(1.0-e*e)*R_dot)/(norm_vector(h)*e));
    // ac_f = acos((1.0/e)*((a*(1.0-e*e))/norm_vector(x)-1.0));
    // if (ac_f > M_PI / 2.0)
    // {
    //     if (f > 0.0)
    //     {
    //         f += 2.0 * (M_PI/2.0 - f);
    //     }
    //     else
    //     {
    //         f -= 2.0 * (-M_PI/2.0 - f);
    //     }
    // }

    // // argument of periapsis
    // double w_plus_f, ac_w_plus_f;
    // w_plus_f = asin(x[2]/(norm_vector(x)*sin(I)));
    // ac_w_plus_f = acos((1.0/cos(Omega))*((x[0]/norm_vector(x))+sin(Omega)*(x[2]/(norm_vector(x)*sin(I)))*cos(I)));
    // if (ac_w_plus_f > M_PI / 2.0)
    // {
    //     if (w_plus_f > 0.0)
    //     {
    //         w_plus_f += 2.0 * (M_PI/2.0 - w_plus_f);
    //     }
    //     else
    //     {
    //         w_plus_f -= 2.0 * (-M_PI/2.0 - w_plus_f);
    //     }
    // }
    // w = w_plus_f - f;

    // // eccentric anomaly
    // double E = 2.0 * atan2( tan(0.5*f), sqrt((1.0 + e)/(1.0 - e)) );

    // // mean anomaly
    // M = E - e * sin(E);

    /* Mercury */
    // standard gravitational parameter
    double mu = G * (m1 + m2);

    // orbital momentum vector
    double h_vec[3];
    cross_product(h_vec, x, v);

    // auxiliary variables
    double r = norm_vector(x);
    double h = norm_vector(h_vec);
    double v2 = norm_squared_vector(v);
    double h2 = norm_squared_vector(h_vec);
    double rv = dot_product(x, v);
    double s = h2 / mu;

    // semi-major axis
    a = (mu *r) / (2.0 * mu - r * v2);

    // inclination and node
    double ci = h_vec[2] / h;
    if (fabs(ci) < 1.0)
    {
        I = acos(ci);
        Omega = atan2(h_vec[0],-h_vec[1]);
        if (Omega < 0.0) Omega += 2.0 * M_PI;
    }
    else
    {
        if (ci < 0.0)
        {
            I = M_PI;
        }
        else
        {
            I = 0.0;
        }
        Omega = 0.0;
    }

    // eccentricity
    double temp = 1.0 + s * (v2 / mu  -  2.0 / r);
    if (temp > 0.0)
    {
        e = sqrt(temp);
    }
    else
    {
        e = 0.0;
    }

    // true longitude
    double to;
    double temp2;
    double true_long;
    if (fabs(h_vec[1]) > 1e-15)
    {
        to = -h_vec[0]/h_vec[1];
        temp = (1.0 - ci) * to;
        temp2 = temp * temp;
        true_long = atan2((x[1]*(1.0+temp2*ci)-x[0]*temp),(x[0]*(temp2+ci)-x[1]*temp));
    }
    else
    {
        true_long = atan2(x[1] * ci, x[0]);
    }
    if (ci < 0.0) true_long += M_PI;

    double p; // longitude of perihelion
    double ce;
    double bige;
    double cf;
    double f;
    if (e < 3e-8)
    {
        p = 0.0;
        M = true_long;
    }
    else
    {
        ce = (v2 * r - mu) / (e * mu);
        if (e < 1.0) // mean anomaly for ellipse
        {
            if (fabs(ce) > 1.0)
            {
                if (ce < 0.0)
                {
                    ce = -1.0;
                }
                else
                {
                    ce = 1.0;
                }
            }
            bige = acos(ce);
            if (rv < 0.0) bige = 2.0 * M_PI - bige;
            M = bige - e * sin(bige);
        }
        else // mean anomaly for hyperbola
        {
            if (ce < 1.0) ce = 1.0;
            bige = log ( ce + sqrt(ce * ce - 1.0) );
            if (rv < 0.0) bige *= -1.0;
            M = e * sinh(bige) - bige;
        }
        // longitude of perihelion
        cf = (s - r) / (e * r);
        if (fabs(cf) > 1.0)
        {
            if (cf < 0.0)
            {
                cf = -1.0;
            }
            else
            {
                cf = 1.0;
            }
        }
        f = acos(cf);
        if (rv < 0.0) f = 2.0 * M_PI - f;
        p = true_long - f;
        p = fmod(p + 4.0 * M_PI, 2.0 * M_PI);
    }

    w = p - Omega;

    if (M < 0.0) M += 2.0 * M_PI;
    if (M > 2.0 * M_PI) M = fmod(M, 2.0 * M_PI);

    // writting orbital elements to body (in radians)
    (*body).a = a;
    (*body).e = e;
    (*body).I = I;
    (*body).M = M;
    (*body).w = w;
    (*body).Omega = Omega;

    return 0;
}

int
initialize_angular_velocity_on_z_axis(cltbdy *body)
{
    double omega_on_body[] = {0.0, 0.0, 2.0 * M_PI / (*body).lod};
    rotate_vector_with_quaternion((*body).omega,
        (*body).q, omega_on_body);

    return 0;
}

int
initialize_angular_velocity(cltbdy *body)
{
    double omega_on_body_spherical[] = 
        {2.0 * M_PI / (*body).lod, 
         (*body).azi,
         (*body).pol};
    double omega_on_body[] = {0.0, 0.0, 0.0};
    spherical_to_cartesian_coordinates(omega_on_body, 
        omega_on_body_spherical);
    rotate_vector_with_quaternion((*body).omega,
        (*body).q, omega_on_body);

    return 0;
}

int
calculate_obliquity_free_body_from_angular_velocity(cltbdy *body)
{
    double e_z[] = {0.0, 0.0, 1.0};

    body->obl = angle_between_two_vectors(body->omega, e_z);

    return 0;
}

int
calculate_obliquity_on_orbit_from_angular_velocity(cltbdy *body)
{
    double h[3]; // orbital momentum vector
    cross_product(h, body->relative_x, body->relative_x_dot);

    body->obl = angle_between_two_vectors(body->omega, h);

    return 0;
}

int
calculate_obliquity_free_body_from_figure_axis_of_solid_frame(cltbdy *body)
{
    double e_z[] = {0.0, 0.0, 1.0};

    double figure_axis_on_body[] = {0.0, 0.0, 1.0};
    double figure_axis[3];
    rotate_vector_with_quaternion(figure_axis,
        body->q, figure_axis_on_body);

    body->obl = angle_between_two_vectors(figure_axis, e_z);

    return 0;
}

int
calculate_obliquity_on_orbit_from_figure_axis_of_solid_frame(cltbdy *body)
{
    double h[3]; // orbital momentum vector
    cross_product(h, body->relative_x, body->relative_x_dot);

    double figure_axis_on_body[] = {0.0, 0.0, 1.0};
    double figure_axis[3];
    rotate_vector_with_quaternion(figure_axis,
        body->q, figure_axis_on_body);

    body->obl = angle_between_two_vectors(figure_axis, h);

    return 0;
}

double
angle_between_spin_axis_and_figure_axis_of_solid_frame(const cltbdy body)
{
    double figure_axis_on_body[] = {0.0, 0.0, 1.0};
    double figure_axis[3];
    rotate_vector_with_quaternion(figure_axis,
        body.q, figure_axis_on_body);

	return angle_between_two_vectors(body.omega, figure_axis);
}

double
angle_between_spin_axis_and_figure_axis(const cltbdy body)
{
	double P[9];
	calculate_eigenvectors_matrix(P, body.b);
    double figure_axis_on_body[] = {0.0, 0.0, 1.0};
    double figure_axis[3];
    square_matrix_times_vector(figure_axis,
        P, figure_axis_on_body);

	return angle_between_two_vectors(body.omega, figure_axis);
}

double
angle_between_spin_axis_and_angular_momentum(const cltbdy body)
{
	return angle_between_two_vectors(body.omega, body.l);
}

int
calculate_center_of_mass(double center_of_mass[3],
				 		 const cltbdy *bodies,
			 	 		 const int number_of_bodies,
			 	 		 const double G)
{
	double sum_of_masses = 0.0;
	double sum_of_masses_times_positions[] = {0.0,0.0,0.0};

	for (int i = 0; i < number_of_bodies; i++)
	{
		sum_of_masses += bodies[i].mass;
		double mass_times_position[3];
		scale_vector(mass_times_position, bodies[i].mass, bodies[i].x);
		linear_combination_vector(sum_of_masses_times_positions,
			1.0, sum_of_masses_times_positions,
			1.0, mass_times_position);
	}

	scale_vector(center_of_mass, 
		1.0 / sum_of_masses, sum_of_masses_times_positions);

	return 0;
}

double
largest_time_scale	(const cltbdy *bodies,
			 	 	 const int number_of_bodies,
			 	 	 const double G)
{
    double largest = 0.0;
    double largest_orbital_period = 0.0;
    double largest_rotational_period = 0.0;
    double largest_rheology_time_scale = 0.0;

    for (int i = 0; i < number_of_bodies; i++)
    {
        if (i > 0)
        {
            double body_orbital_period = 
                kepler_period(bodies[0].mass, bodies[i].mass, 
                    G, bodies[i].a);
            if (body_orbital_period > largest_orbital_period)
            {
                largest_orbital_period = body_orbital_period;
            }
        }
        double body_rotational_period = bodies[i].lod;
        if (body_rotational_period > largest_rotational_period)
        {
            largest_rotational_period = body_rotational_period;
        }
        double body_rheology_max_time = bodies[i].tau;
        if (body_rheology_max_time > largest_rheology_time_scale)
        {
            largest_rheology_time_scale = body_rheology_max_time;
        }
    }

    if (largest_orbital_period > largest)
    {
        largest = largest_orbital_period;
    }
    if (largest_rotational_period > largest)
    {
        largest = largest_rotational_period;
    }
    if (largest_rheology_time_scale > largest)
    {
        largest = largest_rheology_time_scale;
    }

    return largest;
}
