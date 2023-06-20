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

typedef struct CelestialBody{
    char    name[100];  // name
	double  mass;       // mass (sun mass)
	double	lod;	    // length of day (day)
	double 	obl;	    // obliquity (deg)
	double	psi;	    // ? (deg)
	double	R;		    // radius (km)
	double	rg;		    // ?
	double	J2;		    //
	double	C22;        //
	double	lib;	    // ? (deg)
	double	k2;		    // Love number
	double	Dt;		    // ? (s)
	double	tau;	    // characteristic time (yr)
	double	a;		    // semi-major axis (AU)
	double	e;		    // orbital eccentricity
	double	I;		    // ? (deg)
	double	M;		    // ? (deg)
	double	w;		    // ? (deg)
	double	Omega;	    // ? (deg)
} cltbdy;

int
count_columns(const char *s);

int
read_input(cltbdy **body, int *number_of_bodies, const char file[]);

int
convert_input(const char file[]);

#endif