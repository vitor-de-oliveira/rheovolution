/**
 * Library for processing input and output data
**/

#ifndef DAT_H
#define DAT_H

#define _GNU_SOURCE // getline
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h> // ssize_t

#include "linear_algebra.h"
#include "celestial_mechanics.h"
#include "tidal_theory.h"

int
count_columns(const char *s);

/* struct for simulation info */

typedef struct SimulationInfo {

	/* simulation id */
	char	name[300];

	/* argument of main */
	char	main_input[300];

	/* files location */
	char	input_folder[300];
	char	output_folder[300];
	char	system_specs[300];
	char	integration_specs[300];
	char	dev_specs[300];

	/* general specs */
	double	G;							// gravitational parameter
	char	units[100];					// units used for simulation
	char	rheology_model[100];		// system file type
	int		number_of_bodies;			// number of bodies to be used
	bool	omega_correction;			// corrects initial angular velocity
	bool	keplerian_motion;			// restricts all orbits to keplerian
	bool	two_bodies_aprox;			// removes interaction between
										// orbiting bodies

	/* numerical specs given by user */
	double	t_init;						// initial time
	double	t_trans;					// transient time
	double	t_final;					// final time
	double	t_step;						// time step
	bool	t_step_received;			// true if user provided t_step

	/* numerical specs defined by the program */
	double 	t_step_init;				// initial time step
	double	t_step_min;					// minimum time step
	double	error_abs;					// absolute error
	double	error_rel;					// relative error

	/* output specs */
	double 	output_size;				// size of main output file
	int		data_skip;					// number of data points
										// to be skipped on printing

	/* auxiliary variables */
	int		counter;					// counter for data skipping
	double	t;							// simulation time
	int		omega_correction_counter;
	double	omega_correction_t_final;
	bool	write_to_file;

} siminf;

int
print_SimulationInfo(siminf simulation);

int
parse_input(siminf *simulation,
			const char input_file[]);

int
fill_in_bodies_data	(cltbdy	**bodies,
				 	 const siminf simulation);

/* output handling */

// calculates how much data to skip
// based on file size of main output
int
calculate_data_skip (siminf *simulation,
					 const cltbdy *bodies);

int
create_output_files	(const cltbdy *bodies,
			 		 const siminf simulation,
					 FILE *out[]);

int
write_output(const cltbdy *bodies,
			 const siminf simulation,
			 FILE *out[]);

int
close_output_files	(const siminf simulation,
					 FILE *out[]);

int
write_simulation_overview	(const int time_spent_in_seconds,
							 const siminf simulation);

// reads the output of the program
// and calculates orbital elements
int
output_to_orbit(cltbdy *bodies,
			 	const siminf simulation);

// reads the output of the program
// and calculates spin variables
int
output_to_spin	(cltbdy *bodies,
				 const siminf simulation);

// reads the output of the program
// as well as the ourput for orbit and
// spin calculations, and plot the
// results via Gnuplot
int
plot_output_comma_orbit_and_spin(const cltbdy *bodies,
								 const siminf simulation);

// returns shortest time scale for given bodies
double
find_shortest_time_scale(const cltbdy *bodies,
					 	 const siminf simulation);

// returns largest time scale for given bodies
double
find_largest_time_scale(const cltbdy *bodies,
					 	const siminf simulation);

#endif
