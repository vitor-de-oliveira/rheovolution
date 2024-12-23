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
    body_dest->orb = body_src.orb;
    body_dest->rot = body_src.rot;
    body_dest->rot_ini = body_src.rot_ini;
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
	printf("orb = ");	
	printf("%1.10e\n", body.orb);
	printf("rot = ");	
	printf("%1.10e\n", body.rot);
	printf("rot_ini = ");	
	printf("%1.10e\n", body.rot_ini);
	printf("azi = ");	
	printf("%1.10e\n", body.azi);
	printf("pol = ");	
	printf("%1.10e\n", body.pol);

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
	
    double Bs[9];
	construct_traceless_symmetric_matrix(Bs, body.Bs_me);
	printf("Bs = \n");
	print_square_matrix(Bs);
    double P[9];
	construct_traceless_symmetric_matrix(P, body.P_me);
	printf("P = \n");
	print_square_matrix(P);
    double bs[9];
	construct_traceless_symmetric_matrix(bs, body.bs_me);
	printf("bs = \n");
	print_square_matrix(bs);
	double p[9];
	construct_traceless_symmetric_matrix(p, body.p_me);
	printf("p = \n");
	print_square_matrix(p);

	printf("q = \n");
	print_quaternion(body.q);
	printf("Y = \n");
	print_square_matrix(body.Y);
	printf("Y_trans = \n");
	print_square_matrix(body.Y_trans);

	printf("omega = ");
	print_vector(body.omega);
	printf("b = \n");
	print_square_matrix(body.b);

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

    // semi-major axis
    a = calculate_semi_major_axis(G, m1, m2, x, v);

    // eccentricity
    e = calculate_eccentricity(G, m1, m2, x, v);

    // inclination
    I = calculate_inclination(G, m1, m2, x, v);

    // mean anomaly
    M = calculate_mean_anomaly(G, m1, m2, x, v);

    // argument of periapsis
    w = calculate_argument_of_periapsis(G, m1, m2, x, v);

    // longitude of the ascending node
    Omega = calculate_longitude_of_the_ascending_node(G, m1, m2, x, v);
    
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
    double omega_on_body[] = {0.0, 0.0, 2.0 * M_PI / (*body).rot_ini};
    rotate_vector_with_quaternion((*body).omega,
        (*body).q, omega_on_body);

    return 0;
}

int
initialize_angular_velocity(cltbdy *body)
{
    double omega_on_body_spherical[] = 
        {2.0 * M_PI / (*body).rot_ini, 
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
calculate_obliquity_on_orbit_from_angular_velocity  (cltbdy *body,
                                                     const cltbdy body_ref)
{
    double relative_x[3];
	linear_combination_vector(relative_x,
		1.0, (*body).x,
		-1.0, body_ref.x);
	double relative_x_dot[3];
	linear_combination_vector(relative_x_dot,
		1.0, (*body).x_dot,
		-1.0, body_ref.x_dot);

    double h[3]; // orbital momentum vector
    cross_product(h, relative_x, relative_x_dot);

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
calculate_obliquity_on_orbit_from_figure_axis_of_solid_frame(cltbdy *body,
                                                             const cltbdy body_ref)
{
    double relative_x[3];
	linear_combination_vector(relative_x,
		1.0, (*body).x,
		-1.0, body_ref.x);
	double relative_x_dot[3];
	linear_combination_vector(relative_x_dot,
		1.0, (*body).x_dot,
		-1.0, body_ref.x_dot);

    double h[3]; // orbital momentum vector
    cross_product(h, relative_x, relative_x_dot);

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

double
angle_between_relative_x_and_I1(const cltbdy body,
                                const cltbdy body_ref)
{
    double relative_x[3];
	linear_combination_vector(relative_x,
		1.0, body.x,
		-1.0, body_ref.x);

	double P[9];
	calculate_eigenvectors_matrix(P, body.b);
    double I1_on_body[] = {1.0, 0.0, 0.0};
    double I1[3];
    square_matrix_times_vector(I1, P, I1_on_body);

	return angle_between_two_vectors(relative_x, I1);
}

int
calculate_center_of_mass(double location_of_center_of_mass[3],
                         double velocity_of_center_of_mass[3],
				 		 const cltbdy *bodies,
			 	 		 const int number_of_bodies)
{
	double sum_of_masses_times_positions[] = {0.0,0.0,0.0};
    double sum_of_masses_times_velocities[] = {0.0,0.0,0.0};
	double total_mass = 0.0;

	for (int i = 0; i < number_of_bodies; i++)
	{
		linear_combination_vector(sum_of_masses_times_positions,
			1.0, sum_of_masses_times_positions,
			bodies[i].mass, bodies[i].x);
		linear_combination_vector(sum_of_masses_times_velocities,
			1.0, sum_of_masses_times_velocities,
			bodies[i].mass, bodies[i].x_dot);
        total_mass += bodies[i].mass;
    }

	scale_vector(location_of_center_of_mass, 
		1.0 / total_mass, sum_of_masses_times_positions);
	scale_vector(velocity_of_center_of_mass, 
		1.0 / total_mass, sum_of_masses_times_velocities);

	return 0;
}
