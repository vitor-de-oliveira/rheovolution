/**
 * Library for reading an input file for Celestial Mechanics
 * written in a standard format and converting the variables
 * to the ones I need for Ragazzo and Ruiz' theory
 * 
 * Author: Vitor M. de Oliveira
 * Date: 19 june 2023
*/

#ifndef CON_H
#define CON_H

#define _GNU_SOURCE // getline
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h> // ssize_t

#include "algelin3d.h"
#include "celmec.h"

typedef struct CelestialBody{
    char    name[100];  // name
	double  mass;       // mass (sun mass)
	double	lod;	    // length of day (day)
	double 	obl;	    // obliquity (deg)
	double	psi;	    // longitude of the ascending node 
						// of the body equator or precession angle (deg)
	double	R;		    // radius (km)
	double	rg;		    // moment of inertia factor
	double	J2;		    // gravity field coefficient
	double	C22;        // gravity field coefficient
	double	lib;	    // angle between the ascending node and the
						// lowest principal moment of inertia (deg)
	double	kf;		    // fluid Love number
	double	Dt;		    // tidal lag (s)
	double	tau;	    // Maxwell relaxation time plus tidal lag (yr)
	double	a;		    // semi-major axis (AU)
	double	e;		    // orbit eccentricity
	double	I;		    // inclination (deg)
	double	M;		    // mean anomaly (deg)
	double	w;		    // argument of periapsis (deg)
	double	Omega;	    // longitude of the ascending node (deg)

	bool	centrifugal; // centrifugal force
	bool	tidal;		 // tidal force
} cltbdy;

int
count_columns(const char *s);

int
read_input(cltbdy **body, int *number_of_bodies, const char file[]);

int
convert_input	(double *m1, double *m2, double *I0, double *R,
				 double *kf, double omega[],
				 double *alpha, double *eta,
				 double tilde_x[], double tilde_x_dot[],
				 bool *centrifugal, bool *tidal,
				 const double G,
				 const char file[]);

#endif