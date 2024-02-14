/**
 * Code designed to calibrate the Love number k2 in various situations
 * and plot this paramter as a function of perturbing frequency
*/

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_rng.h>

#include "../include/tidal_theory.h"

/* semi-major axis growth in the Earth-Moon system */

double
calibrate_Imk2(const double rate, const double dist, 
	const double m1, const double m2, const double I0, 
	const double R,	const double omega_z, const double G);

int
calculate_tau_v_and_tau(double tau_v_pair[2], double tau_pair[2],
	const double nu, const double Imk2, const double dist,
	const double m1, const double m2, const double kf, 
	const double omega_z, const double G);

int
calculate_tau_v_and_tau_from_Rek2(double *tau_v, 
	double *tau, const double Rek2, const double nu,
	const double kf, const double sigma);

int
calculate_k2(double *re, double *im, const double sigma, 
	const double kf, const double tau_v, const double tau);

/***** Implementing gV ******/

typedef struct gV_rheology
{
	int		m;
	double	gamma;
	double 	alpha;
	double  eta;
	double *alpha_k;
	double *eta_k;
	double *sigma;
	double *tau_a;
	double *tau_b;
} gvrheo;

double*
convert_parameters_gV_summations_C_and_D(gvrheo gV,
										 double freq);

int
convert_parameters_gV_f	(const gsl_vector * x,
						 void *params,
            			 gsl_vector * f);

int
convert_parameters_gV_df(const gsl_vector * x, 
						 void *params,
               			 gsl_matrix * J);

int
convert_parameters_gV_fdf 	(const gsl_vector * x, 
							 void *params,
                			 gsl_vector * f,
							 gsl_matrix * J);

int
print_state_f (size_t iter,
			   gsl_multiroot_fsolver * s);

int
print_state_fdf (size_t iter,
			 	 gsl_multiroot_fdfsolver * s);

int
convert_parameters_gV	(gvrheo *gV);

int
calculate_tau_a_tau_b	(double *tau_a,
						 double *tau_b,
						 double k0,
						 double sigma,
						 double Re_k2,
						 double Im_k2);

int
main(int argc, char *argv[])
{
	// double G = 4.0 * M_PI * M_PI;
	// double m1 = 3.00317e-6;
	// double m2 = 3.69464e-8;
	// double a = 2.56955e-3;
	// double I0 = 0.3308;
	// double R = 6.3710e+3 / 1.495978707e8;
	// double omega_z = (2.0 * M_PI) / (0.9972708 / 365.25);
	// double kf = 0.933;

	/* calibration */
	// double rate = 2.54014310646e-13; // 3.8 cm/yr in AU/yr
	// double Imk2;
	// Imk2 = calibrate_Imk2(rate, a, m1, m2, I0, R, omega_z, G);
	// printf("b = %1.15e\n", Imk2);
	// double nu = 0.680; // tau_v = nu * tau
	// double tau_v_pair[2];
	// double tau_pair[2];
	// calculate_tau_v_and_tau(tau_v_pair, tau_pair, 
	// 	nu, Imk2, a, m1, m2, kf, omega_z, G);
	// printf("tau_plus = %1.15e tau_v_plus = %1.15e = %1.15e s\n", tau_pair[1],
	// 	tau_v_pair[1], tau_v_pair[1] * (365.25 * 24.0 * 60.0 * 60.0));
	// printf("tau_minus = %1.15e tau_v_minus = %1.15e = %1.15e s\n", tau_pair[0],
	// 	tau_v_pair[0], tau_v_pair[0] * (365.25 * 24.0 * 60.0 * 60.0));


	// double M = m1 + m2;
	// double n = sqrt((G * M) / pow(a, 3.0));
	// double omega_SD = 2.0 * (omega_z - n); // Semi-diurnal
	// double real, imag;
	// calculate_k2(&real, &imag, omega_SD, kf, tau_v_pair[1], tau_pair[1]);
	// printf("real = %1.5e imag = %1.15e\n", real, imag);

	/* Love number as a function of frequency */
	// for (double sigma_loop = 1e-10; sigma_loop < 1e20; sigma_loop *= 1.01)
	// {
	// 	double real_loop, imag_loop;
	// 	/* measure k2 for lower tau */
	// 	calculate_k2(&real_loop, &imag_loop, sigma_loop, kf, 
	// 		tau_v_pair[0], tau_pair[0]);
	// 	printf("%1.15e %1.15e %1.15e %1.15e ", 
	// 		sigma_loop, real_loop, -1.0 * imag_loop,
	// 		sqrt(real_loop * real_loop + imag_loop * imag_loop));
	// 	/* measure k2 for higher tau */
	// 	calculate_k2(&real_loop, &imag_loop, sigma_loop, kf, 
	// 		tau_v_pair[1], tau_pair[1]);
	// 	printf("%1.15e %1.15e %1.15e\n", 
	// 		real_loop, -1.0 * imag_loop,
	// 		sqrt(real_loop * real_loop + imag_loop * imag_loop));
	// }

	/* Chandler wobble */
	// double sigma_chandler = (2.0 * M_PI) / (433.0 / 365.25); // rad/yr
	// printf("Chandler wobble:\n");
	// printf("Chandler wobble frequency Cwf = %1.5e\n", sigma_chandler);
	// double real, imag;
	// calculate_k2(&real, &imag, sigma_chandler, kf, 
	// 	tau_v_pair[0], tau_pair[0]);
	// printf("tau = %1.5e\n", tau_pair[0]);
	// printf("Re(k2) = %1.5e Im(k2) = %1.5e |k_2| = %1.5e\n", 
	// 	real, imag, sqrt(real*real+imag*imag));
	// calculate_k2(&real, &imag, sigma_chandler, kf, 
	// 	tau_v_pair[1], tau_pair[1]);
	// printf("tau = %1.5e\n", tau_pair[1]);
	// printf("Re(k2) = %1.5e Im(k2) = %1.5e |k_2| = %1.5e\n", 
	// 	real, imag, sqrt(real*real+imag*imag));

	// double Rek2 = 0.358;
	// double kf = 0.933;
	// double nu = 0.680;
	// double sigma = (2.0 * M_PI) / (433.0 / 365.25);

	// double tau, tau_v;
	// calculate_tau_v_and_tau_from_Rek2(&tau_v, &tau,
	// 	Rek2, nu, kf, sigma);
	// printf("tau = %1.5e tau_v = %1.5e = %1.5e s\n", tau,
	// 	tau_v , tau_v * (365.25 * 24.0 * 60.0 * 60.0));

	// double real, imag;
	// calculate_k2(&real, &imag, sigma, kf, tau_v, tau);
	// printf("real = %1.5e imag = %1.15e\n", real, imag);

	// double real_k2, imag_k2;
	// double sigma_CWob = (2.0 * M_PI) / (433.0 / 365.25);
	// double k_0 = 0.35761431522;
	// double tau_v = 3.953e-4;
	// double tau = 1.707e-3;

	// calculate_k2(&real_k2, &imag_k2, sigma_CWob, k_0, 
	// 		tau_v, tau);
	// printf("real = %1.5e imag = %1.15e\n", real_k2, imag_k2);
	// exit(45);

	/* generalised Voigt */

	// double tau_a, tau_b;

	// double Omega 	= (2.0 * M_PI) / (0.9972708 / 365.25);
	// double m1 		= 3.00317e-6;
	// double m2 		= 3.69464e-8;
	// double r 		= 2.56955e-3;
	// double M 		= m1 + m2;
	// double G 		= 4.0 * M_PI * M_PI;
	// double r3		= pow(r, 3.0);
	// double n 		= sqrt((G * M) / r3);

	// double sigma 	= 2.0 * (Omega - n); // Semi-diurnal
	// double k0 		= 0.933;
	// double Im_k2 	= -2.546376867228876e-02;
	// double Re_k2 	= k0 * (1.0 + 
	// 	sigma*sigma*(5.609836775924446e-03-3.814689007628624e-03)*5.609836775924446e-03)
	// 	/ (1.0 + sigma*sigma*5.609836775924446e-03*5.609836775924446e-03);

	// calculate_tau_a_tau_b(&tau_a, &tau_b, k0, sigma, Re_k2, Im_k2);

	double k0 = 0.933;
	// k0 = 0.35761431522;

	// my calibration for semi-diurnal freq
	double sigma_SD = 4434.21;
	double Im_k2_SD = -2.54638e-2;
	double Re_k2_SD = 0.299584;

	double tau_a_SD, tau_b_SD;
	
	// tau_a_SD = 5.60984e-3 - 3.81469e-3;
	// tau_b_SD = 5.60984e-3;

	calculate_tau_a_tau_b(&tau_a_SD, &tau_b_SD, k0, 
		sigma_SD, Re_k2_SD, Im_k2_SD);

	// M2(L) from Table 1 on Ragazzo 2017
	double sigma_M2L 	= (2.0 * M_PI) / (12.421 / (365.25 * 24.0));
	double Im_k2_M2L 	= -0.02496;
	double Re_k2_M2L 	= 0.2811;

	double tau_a_M2L, tau_b_M2L;
	
	calculate_tau_a_tau_b(&tau_a_M2L, &tau_b_M2L, k0, 
		sigma_M2L, Re_k2_M2L, Im_k2_M2L);

	// Chandler Wobble approx from Fig. 10 Ragazzo 2022
	double sigma_CW 	= (2.0 * M_PI) / (433.0 / 365.25);
	double Im_k2_CW 	= -0.002; // taken by eye from plot
	double Re_k2_CW 	= 0.358;
	
	double tau_a_CW, tau_b_CW;

	calculate_tau_a_tau_b(&tau_a_CW, &tau_b_CW, k0, 
		sigma_CW, Re_k2_CW, Im_k2_CW);

	// Chandler Wobble from Chen 2023
	// double sigma_CW 	= (2.0 * M_PI) / (430.4 / 365.25);
	// double Im_k2_CW 	= -0.00226238;
	// double Re_k2_CW 	= 0.35010616;
	
	// double tau_a_CW, tau_b_CW;

	// calculate_tau_a_tau_b(&tau_a_CW, &tau_b_CW, k0, 
	// 	sigma_CW, Re_k2_CW, Im_k2_CW);
	
	// printf ("tau_a_CW = %1.5e\n", tau_a_CW);
	// printf ("tau_b_CW = %1.5e\n", tau_b_CW);
	// printf ("tau_e_CW / tau_CW  = %1.5e\n", tau_a_CW / tau_b_CW);
	// double year_to_s = 365.25 * 24.0 * 60.0 * 60.0;
	// printf ("tau_v_CW = %1.5e s\n", (tau_b_CW - tau_a_CW) * year_to_s);
	// printf ("tau_CW = %1.5e yr\n", tau_b_CW);
	// exit(12);

	//  N2(L) from Table 1 on Ragazzo 2017
	double sigma_N2L = (2.0 * M_PI) / (12.658 / (365.25 * 24.0));
	double Im_k2_N2L = -0.03214;
	double Re_k2_N2L = 0.2825;

	double tau_a_N2L, tau_b_N2L;

	calculate_tau_a_tau_b(&tau_a_N2L, &tau_b_N2L, k0, 
		sigma_N2L, Re_k2_N2L, Im_k2_N2L);

	/*** converting gV ***/

	gvrheo gV;

	double gamma = 1.629106922962863e+09; // result for parameter_gamma() in data_processing.c for Earth
	// gamma = (1.629106922962863e+09+2.6211612974e9); // accounting for mu0 as in Ragazzo 2022

	/* maxwell model */
	// int	   m = 0;
	// double sigma[] = {sigma_M2L};
	// double tau_a[] = {tau_a_M2L};
	// double tau_b[] = {tau_b_M2L};
	// double Re_k2[] = {Re_k2_M2L};
	// double Im_k2[] = {Im_k2_M2L};

	/* burgers model */
	int	   m = 1;
	double sigma[] = {sigma_M2L, sigma_CW};
	double tau_a[] = {tau_a_M2L, tau_a_CW};
	double tau_b[] = {tau_b_M2L, tau_b_CW};
	double Re_k2[] = {Re_k2_M2L, Re_k2_CW};
	double Im_k2[] = {Im_k2_M2L, Im_k2_CW};
	// double sigma[] = {sigma_SD, sigma_CW};
	// double tau_a[] = {tau_a_SD, tau_a_CW};
	// double tau_b[] = {tau_b_SD, tau_b_CW};
	// double Re_k2[] = {Re_k2_SD, Re_k2_CW};
	// double Im_k2[] = {Im_k2_SD, Im_k2_CW};

	/* two elements model */
	// int	  m = 2;
	// double sigma[] = {sigma_M2L, sigma_CW, sigma_N2L};
	// double tau_a[] = {tau_a_M2L, tau_a_CW, tau_a_N2L};
	// double tau_b[] = {tau_b_M2L, tau_b_CW, tau_b_N2L};
	// double Re_k2[] = {Re_k2_M2L, Re_k2_CW, Re_k2_N2L};
	// double Im_k2[] = {Im_k2_M2L, Im_k2_CW, Im_k2_N2L};

	gV.m 	 = m;
	gV.gamma = gamma;
	gV.sigma = sigma;
	gV.tau_a = tau_a;
	gV.tau_b = tau_b;

	if (gV.m > 0)
	{
		gV.alpha_k 	= (double *) malloc(gV.m * sizeof(double));
		gV.eta_k	= (double *) malloc(gV.m * sizeof(double));
	}

	convert_parameters_gV(&gV);

	printf ("gamma = %1.10e\n", gV.gamma);
	printf ("alpha = %1.10e\n", gV.alpha);
	printf ("eta = %1.10e\n", gV.eta);
	for (int i = 0; i < gV.m; i++)
	{
		printf ("alpha_%d = %1.10e\n", i+1, gV.alpha_k[i]);
		printf ("eta_%d = %1.10e\n", i+1, gV.eta_k[i]);
	}

	double sigma_reference = 1.0;
	sigma_reference = sigma_M2L;

	FILE *out_1 = fopen("tests/test_calibration_k2_reference_points.dat", "w");
	for (int i = 0; i < gV.m + 1; i++)
	{
		fprintf(out_1, "%1.15e %1.15e %1.15e %1.15e\n", 
			sigma[i] / sigma_reference, Re_k2[i], -1.0 * Im_k2[i],
			sqrt(Re_k2[i] * Re_k2[i] + Im_k2[i] * Im_k2[i]));
	}
	fclose(out_1);

	/* Love number as a function of frequency */
	FILE *out_2 = fopen("tests/test_calibration_k2_plot.dat", "w");
	for (double sigma_loop = 1e-10; sigma_loop < 1e10; sigma_loop *= 1.01)
	{
		double tau_a_sigma, tau_b_sigma;

		double *C_and_D;

		C_and_D = convert_parameters_gV_summations_C_and_D(gV, sigma_loop);

		tau_a_sigma = (gV.eta / gV.alpha + gV.eta * C_and_D[0]) 
			/ (1.0 + sigma_loop * sigma_loop * gV.eta * C_and_D[1]);
		tau_b_sigma = (gV.eta / gV.alpha + gV.eta / gV.gamma + gV.eta * C_and_D[0]) 
			/ (1.0 + sigma_loop * sigma_loop * gV.eta * C_and_D[1]);

		free(C_and_D);

		double real_loop, imag_loop;

		calculate_k2(&real_loop, &imag_loop, sigma_loop, k0, 
			tau_b_sigma - tau_a_sigma, tau_b_sigma);
		
		fprintf(out_2, "%1.15e %1.15e %1.15e %1.15e\n", 
			sigma_loop / sigma_reference, real_loop, -1.0 * imag_loop,
			sqrt(real_loop * real_loop + imag_loop * imag_loop));
	}
	fclose(out_2);

	if (gV.m > 0)
	{
		free(gV.alpha_k);
		free(gV.eta_k);
	}

	return 0;
}

double
calibrate_Imk2(const double rate, const double dist, 
	const double m1, const double m2, const double I0, 
	const double R,	const double omega_z, const double G)
{
	double Imk2;

	double alpha = rate;
	double r = dist;
	double a = R;
	double Omega = omega_z;

	double M = m1 + m2;
	double m = (m1 * m2) / M;

	double r2 = pow(r, 2.0);
	double r3 = pow(r, 3.0);
	double r5 = pow(r, 5.0);

	double n = sqrt((G * M) / r3);

	double n2 = pow(n, 2.0);

	double n_p = -(3.0 / 2.0) * sqrt((G * M) / r5);
	double Omega_p = -(m / (2.0 * I0)) * sqrt((G * M) / r);

	double h_1 = m * n * n_p * r2;
	double h_2 = m * n2 * r;
	double h_3 = I0 * Omega * Omega_p;
	double h_4 = (G * m1 * m2) / r2;

	double h = h_1 + h_2 + h_3 + h_4;

	double N_2 = sqrt(5.0 / (4.0 * M_PI * 24.0));

	double a_hat = (G * m2) / (2.0 * N_2 * r3);

	double a5 = pow(a, 5.0);
	double a_hat2 = pow(a_hat, 2.0);

	double omega_SD = 2.0 * (Omega - n); // Semi-diurnal

	double beta = -(5.0 * a5 * a_hat2 * omega_SD) / (32.0 * M_PI * G);

	Imk2 = -(alpha * h) / beta;

	return Imk2;
}

int
calculate_tau_v_and_tau(double tau_v_pair[2], double tau_pair[2],
	const double nu, const double Imk2, const double dist,
	const double m1, const double m2, const double kf, 
	const double omega_z, const double G)
{
	double tau_local_minus;
	double tau_local_plus;
	double tau_v_local_minus;
	double tau_v_local_plus;

	double r = dist;
	double b = Imk2;
	double Omega = omega_z;

	double M = m1 + m2;

	double r3 = pow(r, 3.0);

	double n = sqrt((G * M) / r3);

	double omega_tilde = 2.0 * (Omega - n); // Semi-diurnal

	double kf2 = pow(kf, 2.0);
	double nu2 = pow(nu, 2.0);
	double b2 = pow(b, 2.0);

	tau_local_minus = (-kf * nu) / (2.0 * b * omega_tilde) 
		- (1.0 / (2.0 * omega_tilde)) * sqrt(((kf2 * nu2) / b2) - 4.0);
	tau_local_plus = (-kf * nu) / (2.0 * b * omega_tilde) 
		+ (1.0 / (2.0 * omega_tilde)) * sqrt(((kf2 * nu2) / b2) - 4.0);

	tau_v_local_minus = nu * tau_local_minus; 
	tau_v_local_plus = nu * tau_local_plus; 

	tau_v_pair[0] = tau_v_local_minus;
	tau_v_pair[1] = tau_v_local_plus;
	tau_pair[0] = tau_local_minus;
	tau_pair[1] = tau_local_plus;

	return 0;
}

int
calculate_tau_v_and_tau_from_Rek2(double *tau_v, 
	double *tau, const double Rek2, const double nu,
	const double kf, const double sigma)
{
	*tau = sqrt((kf-Rek2)/(sigma*sigma*(Rek2-kf*(1.0-nu))));
	*tau_v = nu * (*tau); 

	return 0;
}

int
calculate_k2(double *re, double *im, const double sigma, 
	const double kf, const double tau_v, const double tau)
{
	double tau_e = tau - tau_v;
	double sigma2 = sigma * sigma;
	double tau2 = tau * tau;

	*re =  kf * ((1.0 + sigma2 * tau_e * tau) / (1.0 + sigma2 * tau2));
	*im = -kf * ((sigma * tau_v) / (1.0 + sigma2 * tau2));

	return 0;
}

double*
convert_parameters_gV_summations_C_and_D(gvrheo gV,
										 double freq)
{
	double C = 0.0, D = 0.0;
	for (int i = 0; i < gV.m; i++)
	{
		double denominator 
			= gV.alpha_k[i] * gV.alpha_k[i]
			+ (freq * gV.eta_k[i]) * (freq * gV.eta_k[i]);
		C += gV.alpha_k[i] / denominator;
		D += gV.eta_k[i] / denominator;
	}

	double *C_and_D = (double *) malloc (2 * sizeof(double));
	C_and_D[0] = C;
	C_and_D[1] = D;

	return C_and_D;
}

int
convert_parameters_gV_f	(const gsl_vector * x,
						 void *params,
            			 gsl_vector * f)
{
	gvrheo *gV = (gvrheo *) params;

	gV->alpha 	= gsl_vector_get (x, 0);
	gV->eta 	= gsl_vector_get (x, 1);

	for (int i = 0; i < gV->m; i++)
	{
		gV->alpha_k[i] 	= gsl_vector_get (x, 2+i*2);
		gV->eta_k[i] 	= gsl_vector_get (x, 2+i*2+1);
	}

	for (int i = 0; i < gV->m + 1; i++)
	{
		double *C_and_D;
		
		C_and_D = convert_parameters_gV_summations_C_and_D(*gV, gV->sigma[i]);

		double component_f = 
			(1.0 + gV->sigma[i] * gV->sigma[i] * gV->eta * C_and_D[1]) * gV->tau_a[i]
				- gV->eta / gV->alpha - gV->eta * C_and_D[0];

		double component_g = 
			(1.0 + gV->sigma[i] * gV->sigma[i] * gV->eta * C_and_D[1]) * gV->tau_b[i]
				- gV->eta / gV->alpha - gV->eta / gV->gamma - gV->eta * C_and_D[0];

		gsl_vector_set (f, i, component_f);
		gsl_vector_set (f, i + (gV->m + 1), component_g);

		free(C_and_D);
	}

	return GSL_SUCCESS;
}

int
convert_parameters_gV_df(const gsl_vector * x, 
						 void *params,
               			 gsl_matrix * J)
{
	gvrheo *gV = (gvrheo *) params;

	for (int i = 0; i < gV->m + 1; i++)
	{
		double *C_and_D;
		
		C_and_D = convert_parameters_gV_summations_C_and_D(*gV, gV->sigma[i]);

		double component_del_f_del_alpha
			= gV->eta / (gV->alpha * gV->alpha);
		double component_del_f_del_eta
			= gV->sigma[i] * gV->sigma[i] * C_and_D[1] * gV->tau_a[i]
				- 1.0 / gV->alpha - C_and_D[0];
		double component_del_g_del_alpha
			= component_del_f_del_alpha;
		double component_del_g_del_eta
			= gV->sigma[i] * gV->sigma[i] * C_and_D[1] * gV->tau_b[i]
				- 1.0 / gV->alpha - 1.0 / gV->gamma - C_and_D[0];		

		free(C_and_D);

		gsl_matrix_set (J, i, 0, component_del_f_del_alpha);
		gsl_matrix_set (J, i, 1, component_del_f_del_eta);
		gsl_matrix_set (J, i + (gV->m + 1), 0, component_del_g_del_alpha);
		gsl_matrix_set (J, i + (gV->m + 1), 1, component_del_g_del_eta);

		for (int j = 0; j < gV->m; j++)
		{
			double denominator 
				= pow(gV->sigma[i] * gV->sigma[i] * gV->eta_k[j] * gV->eta_k[j] 
					+ gV->alpha_k[j] * gV->alpha_k[j], 2.0);
			double numerator_1 
				= gV->sigma[i] * gV->sigma[i] * gV->eta_k[j] * gV->eta_k[j] 
					- gV->alpha_k[j] * gV->alpha_k[j];
			double numerator_2 
				= - 2.0 * gV->sigma[i] * gV->sigma[i] * gV->alpha_k[j] * gV->eta_k[j];
			double numerator_3 
				= - 2.0 * gV->alpha_k[j] * gV->eta_k[j];
			double numerator_4 
				= -1.0 * numerator_1;

			double del_C_del_alpha_k	= numerator_1 / denominator;
			double del_C_del_eta_k 		= numerator_2 / denominator;
			double del_D_del_alpha_k 	= numerator_3 / denominator;
			double del_D_del_eta_k 		= numerator_4 / denominator;

			double component_del_f_del_alpha_k
				= gV->sigma[i] * gV->sigma[i] * gV->eta * gV->tau_a[i] 
					* del_D_del_alpha_k
					- gV->eta * del_C_del_alpha_k;

			double component_del_f_del_eta_k
				= gV->sigma[i] * gV->sigma[i] * gV->eta * gV->tau_a[i] 
					* del_D_del_eta_k
					- gV->eta * del_C_del_eta_k;

			double component_del_g_del_alpha_k
				= gV->sigma[i] * gV->sigma[i] * gV->eta * gV->tau_b[i] 
					* del_D_del_alpha_k
					- gV->eta * del_C_del_alpha_k;

			double component_del_g_del_eta_k
				= gV->sigma[i] * gV->sigma[i] * gV->eta * gV->tau_b[i]
					* del_D_del_eta_k
					- gV->eta * del_C_del_eta_k;

			gsl_matrix_set (J, i, 2 + j*2, component_del_f_del_alpha_k);
			gsl_matrix_set (J, i, 2 + j*2 + 1, component_del_f_del_eta_k);
			gsl_matrix_set (J, i + (gV->m + 1), 2 + j*2, component_del_g_del_alpha_k);
			gsl_matrix_set (J, i + (gV->m + 1), 2 + j*2 + 1, component_del_g_del_eta_k);
		}
	}

	return GSL_SUCCESS;
}

int
convert_parameters_gV_fdf 	(const gsl_vector * x, 
							 void *params,
                			 gsl_vector * f,
							 gsl_matrix * J)
{
	convert_parameters_gV_f (x, params, f);
	convert_parameters_gV_df(x, params, J);

	return GSL_SUCCESS;
}

int
print_state_f (size_t iter,
			   gsl_multiroot_fsolver * s)
{
	printf ("iter = %3lu ", iter);

	printf ("x = ");
	for (size_t i = 0; i < s->x->size; i++)
	{
		printf ("% .5e ", gsl_vector_get (s->x, i));
	}

	printf ("f(x) = ");
	for (size_t i = 0; i < s->f->size; i++)
	{
		printf ("% .5e ", gsl_vector_get (s->f, i));
	}

	printf ("\n");

	return 0;
}

int
print_state_fdf (size_t iter,
			 	 gsl_multiroot_fdfsolver * s)
{
	printf ("iter = %3lu ", iter);

	printf ("\n");

	printf ("x = ");
	for (size_t i = 0; i < s->x->size; i++)
	{
		printf ("% .5e ", gsl_vector_get (s->x, i));
	}

	printf ("\n");

	printf ("f(x) = ");
	for (size_t i = 0; i < s->f->size; i++)
	{
		printf ("% .5e ", gsl_vector_get (s->f, i));
	}

	printf ("\n");

	printf ("J(x) = \n");
	for (size_t i = 0; i < s->J->size1; i++)
	{
		for (size_t j = 0; j < s->J->size2; j++)
		{
			printf ("% .5e ", gsl_matrix_get (s->J, i, j));
		}
		printf ("\n");
	}

	printf ("\n");

	return 0;
}

int
convert_parameters_gV (gvrheo *gV)
{
	const int	 n = 2 * (gV->m + 1);
	const size_t n_size_t = (size_t) n;

	gsl_vector *x = gsl_vector_alloc (n_size_t);

	/* Maxwell solution */
	double alpha_init 	= gV->gamma * (gV->tau_b[0] - gV->tau_a[0]) / gV->tau_a[0];
	double eta_init		= gV->gamma * (gV->tau_b[0] - gV->tau_a[0]);

	if (gV->m == 0)
	{
		gsl_vector_set (x, 0, alpha_init);
		gsl_vector_set (x, 1, eta_init);
	}
	else
	{
		gsl_multiroot_function_fdf fdf = {&convert_parameters_gV_f,
										  &convert_parameters_gV_df,
										  &convert_parameters_gV_fdf,
										  n_size_t, gV};
		const gsl_multiroot_fdfsolver_type *T;
		gsl_multiroot_fdfsolver *s;
		T = gsl_multiroot_fdfsolver_gnewton;
		s = gsl_multiroot_fdfsolver_alloc (T, n_size_t);

		// gsl_multiroot_function f = {&convert_parameters_gV_f,
		// 							n_size_t, gV};
		// const gsl_multiroot_fsolver_type *T2;
		// gsl_multiroot_fsolver *s2;
		// T2 = gsl_multiroot_fsolver_hybrids;
		// s2 = gsl_multiroot_fsolver_alloc (T2, n_size_t);

		bool loop = false;
		int loop_counter = 0, loop_limit = 1000;

 		if (gV->m == 1)
		{
			gsl_vector_set (x, 0, alpha_init / 10000.0);
			gsl_vector_set (x, 1, eta_init * 1000.0); // overshooting for eta
			gsl_vector_set (x, 2, alpha_init / 1000.0);
			gsl_vector_set (x, 3, eta_init / 10000.0);
			/* trying for ps with a0 calibrated */
			// gsl_vector_set (x, 0, alpha_init / 100.0);
			// gsl_vector_set (x, 1, eta_init * 1000.0);
			// gsl_vector_set (x, 2, alpha_init / 10.0);
			// gsl_vector_set (x, 3, eta_init * 10.0);
		}
		else
		{
			/****** does not work (yet) *******/

			loop = true;
			/* turn off error handler so I can loop */
			gsl_set_error_handler_off();

			back:;

			const gsl_rng_type * T_rng;
			gsl_rng * r;

			gsl_rng_env_setup();

			T_rng = gsl_rng_default;
			r = gsl_rng_alloc (T_rng);

			double sum_min = 1e6;
			int interval_max = 22;
			int interval_min = 10;
			for (int k = 0; k < 1000000; k++)
			{
				gsl_vector *x_rng = gsl_vector_alloc (n_size_t);
				gsl_vector *f_rng = gsl_vector_alloc (n_size_t);

				gsl_rng_set (r, time(NULL));

				for (int i = 0; i < n; i++)
				{
					int exponent = (gsl_rng_get (r) % 
						(interval_max - interval_min + 1)) + interval_min;
					gsl_vector_set (x_rng, i, pow (10.0, exponent));
				}

				convert_parameters_gV_f(x_rng, gV, f_rng);

				double sum = 0.0;
				for (size_t i = 0; i < f_rng->size; i++)
				{
					sum += fabs(gsl_vector_get (f_rng, i));
				}
				if (sum < sum_min)
				{
					sum_min = sum;
					gsl_vector_memcpy (x, x_rng);
				}
				gsl_vector_free (x_rng);
				gsl_vector_free (f_rng);
			}
			gsl_rng_free (r);
		}

		gsl_multiroot_fdfsolver_set (s, &fdf, x);

		// gsl_multiroot_fsolver_set (s2, &f, x);
	
		int status_1, status_2;
		size_t iter = 0;

		print_state_fdf (iter, s);

		// print_state_f (iter, s2);

		do
		{
			iter++;

			status_1 = gsl_multiroot_fdfsolver_iterate (s);

			// status_1 = gsl_multiroot_fsolver_iterate (s2);

			print_state_fdf (iter, s);

			// print_state_f (iter, s2);

			int sinal_check = 0;
			for (size_t i = 0; i < s->x->size; i++)
			{
				if(gsl_vector_get (s->x, i) < 0.0)
				{
					sinal_check = -1;
				}
			}

			if ((status_1 || sinal_check) && loop)
			{
				loop_counter++;
				if (loop_counter == loop_limit)
				{
					break;
				}
				else
				{
					goto back;
				}
			}
				
			status_2 = gsl_multiroot_test_residual (s->f, 1e-10);

			// status_2 = gsl_multiroot_test_residual (s2->f, 1e-10);
		}
		while (status_2 == GSL_CONTINUE && iter < 1000);

		if (!loop)
		{
			printf ("status = %s\n", gsl_strerror (status_1));
		}

		gsl_vector_memcpy (x, s->x);

		// gsl_vector_memcpy (x, s2->x);

		gsl_multiroot_fdfsolver_free (s);

		// gsl_multiroot_fsolver_free (s2);
	}

	// update variables
	gV->alpha 	= gsl_vector_get (x, 0);
	gV->eta 	= gsl_vector_get (x, 1);
	for (int i = 0; i < gV->m; i++)
	{
		gV->alpha_k[i]	= gsl_vector_get (x,2+i*2);
		gV->eta_k[i] 	= gsl_vector_get (x,2+i*2+1);
	}

	gsl_vector_free (x);

	return 0;
}

int
calculate_tau_a_tau_b	(double *tau_a,
						 double *tau_b,
						 double k0,
						 double sigma,
						 double Re_k2,
						 double Im_k2)
{
	double tau_a_local = 0.0;
	double tau_b_local = 0.0;

	/* tau_b^3 + a tau_b^2 + b tau_b + c = 0 */
	double a = (k0 - Re_k2) / (sigma * Im_k2);
	double b = 1.0 / (sigma * sigma);
	double c = (k0 - Re_k2) / (sigma * sigma * sigma * Im_k2);

	/* from Wolfram Alpha */
	double d = pow((-2.0*a*a*a + 3.0*sqrt(3.0)*sqrt(4.0*a*a*a*c - a*a*b*b - 18.0*a*b*c 
		+ 4.0*b*b*b + 27.0*c*c) + 9.0*a*b - 27.0*c),(1.0/3.0));

	tau_b_local = d / (3.0 * pow(2.0,(1.0/3.0))) 
		- (pow(2.0,(1.0/3.0)) * (3.0*b-a*a)) / (3.0 * d)
		- a / 3.0;

	tau_a_local = (Re_k2 * (1.0 + sigma * sigma * tau_b_local * tau_b_local) - k0) 
		/ (sigma * sigma * tau_b_local * k0);

	*tau_a = tau_a_local;
	*tau_b = tau_b_local;

	return 0;
}
