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
#include "dynamical_system.h"

int
count_columns(const char *s);

/* struct for simulation info */

typedef struct SimulationInfo {

	/* simlation id */
	char	name[100];

	/* argument of main */
	char	main_input[100];

	/* files location */
	char	input_folder[100];
	char	output_folder[100];
	char	system_specs[100];
	char	simulation_specs[100];
	char	dev_specs[100];

	/* general specs */
	double	G;							// gravitational parameter
	char	units[100];					// units used for simulation
	bool	write_to_file;
	bool	omega_correction;

	/* system specs */
	int		system_file_type;
	int		number_of_bodies;

	/* numerical specs */
	double	t_init;						// initial time
	double	t_trans;					// transient time
	double	t_final;					// final time
	double	t_step;						// time step
	double	eps_abs;					// absolute error
	double	eps_rel;					// relative error

	/* output specs */
	int		data_skip;					// number of data points
										// to be skipped on printing

	/* auxiliary variables */
	int		counter;					// counter for data skipping
	double	t;							// simulation time
	int		omega_correction_counter;
	double	omega_correction_t_final;

} siminf;

int
parse_input(siminf *simulation,
			const char input_file[]);

int
fill_in_bodies_data	(cltbdy	**bodies,
				 	 const siminf simulation);

/* output handling */

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

#endif
