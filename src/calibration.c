/* Calibration of the Maxwell rheology parameters for the Earth 
 * based on the semi-major axis growth in the Earth--Moon system 
 * and plots of k2 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "tidal_theory.h"

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
	// printf("b = %1.5e\n", Imk2);
	// double nu = 0.680; // tau_v = nu * tau
	// double tau_v_pair[2];
	// double tau_pair[2];
	// calculate_tau_v_and_tau(tau_v_pair, tau_pair, 
	// 	nu, Imk2, a, m1, m2, kf, omega_z, G);
	// printf("tau_plus = %1.5e tau_v_plus = %1.5e = %1.5e s\n", tau_pair[1],
	// 	tau_v_pair[1], tau_v_pair[1] * (365.25 * 24.0 * 60.0 * 60.0));
	// printf("tau_minus = %1.5e tau_v_minus = %1.5e = %1.5e s\n", tau_pair[0],
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

	double Rek2 = 0.358;
	double kf = 0.933;
	double nu = 0.680;
	double sigma = (2.0 * M_PI) / (433.0 / 365.25);

	double tau, tau_v;
	calculate_tau_v_and_tau_from_Rek2(&tau_v, &tau,
		Rek2, nu, kf, sigma);
	printf("tau = %1.5e tau_v = %1.5e = %1.5e s\n", tau,
		tau_v , tau_v * (365.25 * 24.0 * 60.0 * 60.0));

	double real, imag;
	calculate_k2(&real, &imag, sigma, kf, tau_v, tau);
	printf("real = %1.5e imag = %1.15e\n", real, imag);

	return 0;
}
