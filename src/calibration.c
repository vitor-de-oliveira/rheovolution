/**
 * Code designed to use additional functions in tidal_theory lib
 * to calibrate the Love number k2 in various situations
 * and plot this paramter as a function of perturbing frequency
*/

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "../include/tidal_theory.h"

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

	// /* calibration */
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

	/***** Implementing gV ******/

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

	// M2(L) from Table 1 on Ragazzo 2017
	double sigma_M2L 	= (2.0 * M_PI) / (12.421 / (365.25 * 24.0));
	double Im_k2_M2L 	= -0.02496;
	double Re_k2_M2L 	= 0.2811;

	double tau_a_M2L, tau_b_M2L;
	
	calculate_tau_a_tau_b(&tau_a_M2L, &tau_b_M2L, k0, 
		sigma_M2L, Re_k2_M2L, Im_k2_M2L);

	// Chandler Wobble approx from Fig. 10 Ragazzo 2022
	double sigma_CW 	= (2.0 * M_PI) / (433.0 / 365.25);
	double Im_k2_CW 	= -0.003; // taken by eye from plot
	double Re_k2_CW 	= 0.358;
	
	double tau_a_CW, tau_b_CW;

	calculate_tau_a_tau_b(&tau_a_CW, &tau_b_CW, k0, 
		sigma_CW, Re_k2_CW, Im_k2_CW);

	int	   m = 1;
	double alpha, eta, alpha_k[m], eta_k[m];

	double gamma 	= 1.629106922962863e+09; // result for parameter_gamma() in data_processing.c for Earth
	double sigma[] 	= {sigma_M2L, sigma_CW};
	double tau_a[] 	= {tau_a_M2L, tau_a_CW};
	double tau_b[] 	= {tau_b_M2L, tau_b_CW};

	convert_parameters_gV(&alpha, &eta, alpha_k, eta_k,
		gamma, m, sigma, tau_a, tau_b);

	// printf ("alpha = %e\n", alpha);
	// printf ("eta = %e\n", eta);
	// for (int i = 0; i < m; i++)
	// {
	// 	printf ("alpha_%d = %e\n", i+1, alpha_k[i]);
	// 	printf ("eta_%d = %e\n", i+1, eta_k[i]);
	// }

	// Love number as a function of frequency
	for (double sigma_loop = 1e-10; sigma_loop < 1e10; sigma_loop *= 1.01)
	{
		// double sigma_loop = sigma_M2L;

		double *C_and_D;
		
		C_and_D = convert_parameters_gV_summations_C_and_D(m,
				sigma_loop, alpha_k, eta_k);

		double C = C_and_D[0];
		double D = C_and_D[1];

		// printf("%e\n", C);
		// printf("%e\n", D);

		free(C_and_D);

		double tau_a_sigma, tau_b_sigma;

		tau_a_sigma = (eta / alpha + eta * C) 
			/ (1.0 + sigma_loop * sigma_loop * eta * D);
		tau_b_sigma = (eta / alpha + eta / gamma + eta * C) 
			/ (1.0 + sigma_loop * sigma_loop * eta * D);

		double real_loop, imag_loop;

		calculate_k2(&real_loop, &imag_loop, sigma_loop, k0, 
			tau_b_sigma - tau_a_sigma, tau_b_sigma);
		printf("%1.15e %1.15e %1.15e %1.15e\n", 
			sigma_loop / sigma_M2L, real_loop, -1.0 * imag_loop,
			sqrt(real_loop * real_loop + imag_loop * imag_loop));

	}

	return 0;
}
