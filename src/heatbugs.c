/*
 * This file is part of heatbugs_CPU.
 *
 * heatbugs_CPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * heatbugs_CPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with heatbugs_CPU. If not, see <http://www.gnu.org/licenses/>.
 * */



#define _GNU_SOURCE	/* this allow 'getopt(..)' function in <unistd.h> to  */
			/* compile under -std=c99 compiler option, because    */
			/* 'getopt' is POSIX, not c99.                        */

#include <stdio.h>	/* printf(...)	*/
#include <stdlib.h>	/* exit(...)	*/
#include <unistd.h>

#include <math.h>	/* fabs(...)	*/
#include <string.h>
#include <ctype.h>

#include "glib.h"	/* FALSE, TRUE, random's */

//#include "keyboard.h"
#include "heatbugs.h"



/** Default parameters. */

#define NUM_ITERATIONS		1000	/* 1000 iterations. (0 = non stop).   */
#define BUGS_NUMBER		100	/* 100 Bugs in the world.             */
#define WORLD_WIDTH		100
#define WORLD_HEIGTH		100

#define WORLD_DIFFUSION_RATE	0.90f	/* [0..1], % temperature to   */
					/* neighbour cells.           */

#define WORLD_EVAPORATION_RATE	0.01f	/* [0..1], % temperature loss */
					/* to 'ether'.                */

#define BUGS_RAND_MOVE_CHANCE	0.00f	/* [0..100], Chance a bug will move.  */
#define BUGS_TEMP_MIN_IDEAL	10	/* [0 .. 200] */
#define BUGS_TEMP_MAX_IDEAL	40	/* [0 .. 200] */
#define BUGS_HEAT_MIN_OUTPUT	5	/* [0 .. 100] */
#define BUGS_HEAT_MAX_OUTPUT	25	/* [0 .. 100] */

/* The file to send results. Directory must exist. */
#define OUTPUT_FILENAME		"../results/heatbugsCPU.csv"


/** Simulation constants. */
#define NUM_NEIGHBOURS 8


/** Used to drive what shall happen to the agent at each step. */
#define FIND_ANY_FREE		0x00ffffff
#define FIND_MAX_TEMPERATURE	0x00ffff00
#define FIND_MIN_TEMPERATURE	0x00ff00ff


/** Heatbugs related. */
#define OKI_DOKI	 0
#define NOT_DOKI	-1

#define RESET 0

#define A_BUG		0x00ff00dd
#define A_EMPTY_CELL	0x00000000


/** Use in swarm_map */
#define HAS_BUG( swarm_map ) ((swarm_map) != A_EMPTY_CELL)
#define HAS_NO_BUG( swarm_map ) ((swarm_map) == A_EMPTY_CELL)

#define NEW_BUG_IN( swarm_map ) swarm_map = A_BUG


/** Use in swarm. */

/* Every time bug change position, use SET_BUG_LOCAL(...) with new position. */
#define SET_BUG_LOCAL( swarm_locus, position ) swarm_locus = position
#define SET_BUG_IDEAL_TEMPERATURE( swarm_iTemp, iTemp ) swarm_iTemp = iTemp
#define SET_BUG_OUTPUT_HEAT( swarm_outHeat, outHeat ) swarm_outHeat = outHeat




/** Input data used for simulation. */
typedef struct parameters {
	/* Num Iterations to stop. (0 = non stop). */
	size_t numIterations;
	/* Number of bugs in the world. */
	size_t bugs_number;
	/* World width size. */
	size_t world_width;
	/* World height size. */
	size_t world_height;
	/* World's vector size = (world_height * world_width). */
	size_t world_size;
	/* [0..1], % temperature to adjacent cells. */
	float world_diffusion_rate;
	/* [0..1], % temperature's loss to 'ether'.  */
	float world_evaporation_rate;
	/* [0..100], Chance a bug will move. */
	float bugs_random_move_chance;
	/* [0 .. 200], bug's minimum prefered temperature. */
	unsigned int bugs_temperature_min_ideal;
	/* [0 .. 200], bug's maximum prefered temperature. */
	unsigned int bugs_temperature_max_ideal;
	/* [0 .. 100], min heat a bug leave in the world in each step. */
	unsigned int bugs_heat_min_output;
	/* [0 .. 100], max heat a bug leave in the world in each step. */
	unsigned int bugs_heat_max_output;
	/* Seed to be used as random generator initialization value. */
	unsigned int seed;
	/* File to send results. */
	char output_filename[256];
} Parameters_t;




/** The bug data type. */
typedef struct bug {
	size_t locus;
	unsigned int ideal_temperature;
	unsigned int output_heat;
} bug_t;


/** Simulation buffers. */
typedef struct hb_buffers {
	bug_t *swarm;			/* SIZE: BUGS_NUM			- Bug's position in the swarm_map. */
	unsigned int *swarm_map;	/* SIZE: WORLD_HEIGHT * WORLD_WIDTH	- Bug's map. Each cell points to NULL (has no bug), or to a bug_t. */
	float *heat_map[2];		/* SIZE: WORLD_HEIGHT * WORLD_WIDTH	- Temperature map (heat_map) & the buffer (heat_buffer). */
	float *unhappiness;		/* SIZE: NUM_BUGS			- The Unhappiness vector. */
} HBBuffers_t;


/** Buffers selector. In each step they swap value to indicate the correct 'heat_map' buffer to be sent to kernel. */
typedef struct hb_buffer_select {
	unsigned int main;
	unsigned int secd;
} hbBufferSelect_t;




const char version[] = "Heatbugs simulation for CPU (serial processing) v3.0 with Glib-2.0 randoms.";


/* To address bug's position neighbourhood. */
enum {SW = 0, S, SE, W, E, NW, N, NE};



#define HB_ERROR hb_error_quark()


static GQuark hb_error_quark( void ) {
	return g_quark_from_static_string( "hb-error-quark" );
}



/**
 * Sets the parameters passed as command line arguments.
 * If there are no parameters, default parameters are used.
 * Default parameter will be used for every omitted parameter.
 *
 * This function uses the GNU 'getopt' command line argument parser,
 * and requires the macro _GNU_SOURCE to be defined. As such, the
 * code using the 'getopt' function is not portable.
 * Consequence of using 'getopt' is that the previous result of 'argv'
 * parameter may change after 'getopt' is used, therefore 'argv'
 * should not be used again.
 * The function 'getopt' is marked as Thread Unsafe.
 *
 * @param[out]	params		- Parameters to be filled with default or
 *				  from command line.
 * @param[in]	argc		- Command line argument counter.
 * @param[in]	argv		- Command line arguments.
 * @param[out]	err		- GLib object for error reporting.
* */
void getSimulParameters( Parameters_t *const params, int argc,
					char *argv[], GError **err )
{
	FILE *uranddev = NULL;

	int c;	/* Parsed command line option. */

	/* The string 't:T:h:H:r:n:d:e:w:W:i:f:' is the parameter string to   */
	/* be checked by 'getopt' function.                                   */
	/* The ':' character means that a value is required after the         */
	/* parameter selector character (i.e. -t 50  or  -t50).               */
	const char matches[] = "t:T:h:H:r:n:d:e:w:W:i:s:f:";


	/* Default / hardcoded parameters. */
	params->numIterations = NUM_ITERATIONS;				/* i */
	params->bugs_number = BUGS_NUMBER;				/* n */
	params->world_width =  WORLD_WIDTH;				/* w */
	params->world_height = WORLD_HEIGTH;				/* W */
	params->world_diffusion_rate = WORLD_DIFFUSION_RATE;		/* d */
	params->world_evaporation_rate = WORLD_EVAPORATION_RATE;	/* e */

	params->bugs_random_move_chance = BUGS_RAND_MOVE_CHANCE;	/* r */
	params->bugs_temperature_min_ideal = BUGS_TEMP_MIN_IDEAL;	/* t */
	params->bugs_temperature_max_ideal = BUGS_TEMP_MAX_IDEAL;	/* T */
	params->bugs_heat_min_output = BUGS_HEAT_MIN_OUTPUT;		/* h */
	params->bugs_heat_max_output = BUGS_HEAT_MAX_OUTPUT;		/* H */

	strcpy( params->output_filename, OUTPUT_FILENAME );		/* f */


	/* Read initial seed from linux /dev/urandom */
	uranddev = fopen( "/dev/urandom", "r" );
	hb_if_err_create_goto( *err, HB_ERROR,
		uranddev == NULL,
		HB_UNABLE_OPEN_FILE, error_handler,
		"Could not open urandom device to get seed." );

	fread( &params->seed, sizeof( unsigned int ), 1, uranddev );
	fclose( uranddev );


	/* Parse command line arguments using GNU's getopt function. */

	while ( (c = getopt( argc, argv, matches )) != -1 )
	{
		switch (c)
		{
			case 't':
				params->bugs_temperature_min_ideal =
					atoi( optarg );
				break;
			case 'T':
				params->bugs_temperature_max_ideal =
					atoi( optarg );
				break;
			case 'h':
				params->bugs_heat_min_output =
					atoi( optarg );
				break;
			case 'H':
				params->bugs_heat_max_output =
					atoi( optarg );
				break;
			case 'r':
				params->bugs_random_move_chance =
					atof( optarg );
				break;
			case 'n':
				params->bugs_number =
					atoi( optarg );
				break;
			case 'd':
				params->world_diffusion_rate =
					atof( optarg );
				break;
			case 'e':
				params->world_evaporation_rate =
					atof( optarg );
				break;
			case 'w':
				params->world_width =
					atoi( optarg );
				break;
			case 'W':
				params->world_height =
					atoi( optarg );
				break;
			case 'i':
				params->numIterations =
					atoi( optarg );
				break;
			case 's':
				params->seed =
					atoi( optarg );
				break;
			case 'f':
				strcpy( params->output_filename, optarg );
				break;
			case '?':
				hb_if_err_create_goto( *err, HB_ERROR,
					(optopt != ':'
					&& strchr( matches, optopt ) != NULL),
					HB_PARAM_ARG_MISSING, error_handler,
					"Option required argument missing." );

				hb_if_err_create_goto( *err, HB_ERROR,
					(isprint( optopt )),
					HB_PARAM_OPTION_UNKNOWN, error_handler,
					"Unknown option." );

				hb_if_err_create_goto( *err, HB_ERROR,
					TRUE,
					HB_PARAM_CHAR_UNKNOWN, error_handler,
					"Unprintable character in command line." );
			default:
				hb_if_err_create_goto( *err, HB_ERROR,
					TRUE,
					HB_PARAM_PARSING, error_handler,
					"Weird error occurred while parsing parameter." );
		} /* end_switch */
	} /* end_while */

	/*
	   NOTE: Check here for extra arguments... see:
	   https://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html
	 */

	params->world_size = params->world_height * params->world_width;

	/* Check for bug's number related errors. */
	hb_if_err_create_goto( *err, HB_ERROR,
		params->bugs_number == 0,
		HB_BUGS_ZERO, error_handler,
		"There are no bugs." );

	hb_if_err_create_goto( *err, HB_ERROR,
		params->bugs_number >= params->world_size,
		HB_BUGS_OVERFLOW, error_handler,
		"Number of bugs exceed available world slots." );

	/* Check range related erros in bug's ideal temperature. */
	/* Checking order matters!                               */
	hb_if_err_create_goto( *err, HB_ERROR,
		params->bugs_temperature_min_ideal >
					params->bugs_temperature_max_ideal,
		HB_TEMPERATURE_OVERLAP, error_handler,
		"Bug's ideal temperature range overlaps." );

	hb_if_err_create_goto( *err, HB_ERROR,
		params->bugs_temperature_max_ideal >= 200,
		HB_TEMPERATURE_OUT_RANGE, error_handler,
		"Bug's max ideal temperature is out of range." );

	/* Check for range related error in bug's output heat. */
	/* Checking order matters!                             */
	hb_if_err_create_goto( *err, HB_ERROR,
		params->bugs_heat_min_output > params->bugs_heat_max_output,
		HB_OUTPUT_HEAT_OVERLAP, error_handler,
		"Bug's output heat range overlaps.");

	hb_if_err_create_goto( *err, HB_ERROR,
		params->bugs_heat_max_output >= 100,
		HB_OUTPUT_HEAT_OUT_RANGE, error_handler,
		"Bug's max output heat is out of range." );

	/* If numeber of bugs is 80% of the world space issue a warning. */
	if (params->bugs_number >= 0.8 * params->world_size)
		fprintf( stderr,
			"Warning: Bugs number near available world slots.\n" );


error_handler:
	/* If error handler is reached leave function imediately. */

	return;
}



/**
 * Create all the buffers for both, host and device.
 *
 * @param[out]	buff		- The structure with all host buffers to be
 *				  filled.
 * @param[in]	params		- The structure with all simulation parameters.
 *				  Used to compute 'buff' sizes.
 * @param[out]	err		- GLib object for error reporting.
 * */
void setupBuffers( HBBuffers_t *const buff, const Parameters_t *const params,
				GError **err )
{
	GError *err_setbuf = NULL;


	/** SWARM. */
	buff->swarm = (bug_t *) malloc( params->bugs_number * sizeof( bug_t ) );
	hb_if_err_create_goto( *err, HB_ERROR,
		buff->swarm == NULL,
		HB_MALLOC_FAILURE, error_handler,
		"Unable to allocate memory for swarm vector." );


	/** SWARM MAP. */
	buff->swarm_map = (unsigned int *) malloc( params->world_size * sizeof( unsigned int ) );
	hb_if_err_create_goto( *err, HB_ERROR,
		buff->swarm_map == NULL,
		HB_MALLOC_FAILURE, error_handler,
		"Unable to allocate memory for swarm map array." );


	/** HEAT MAP & Buffer. */
	buff->heat_map[ 0 ] = (float *) malloc( params->world_size * sizeof( float ) );
	hb_if_err_create_goto( *err, HB_ERROR,
		buff->heat_map[ 0 ] == NULL,
		HB_MALLOC_FAILURE, error_handler,
		"Unable to allocate memory for heat map array." );

	/* Buffer... */
	buff->heat_map[ 1 ] = (float *) malloc( params->world_size * sizeof( float ) );
	hb_if_err_create_goto( *err, HB_ERROR,
		buff->heat_map[ 1 ] == NULL,
		HB_MALLOC_FAILURE, error_handler,
		"Unable to allocate memory for heat map buffer array." );


	/** UNHAPPINESS */
	buff->unhappiness = (float *) malloc( params->bugs_number * sizeof( float ) );
	hb_if_err_create_goto( *err, HB_ERROR,
		buff->unhappiness == NULL,
		HB_MALLOC_FAILURE, error_handler,
		"Unable to allocate memory for unhappiness vector." );


error_handler:
	/* If error hadler is reached leave function imediately. */

	return;
}



/**
 * Initiate the world and create agents.
 *
 * @param[out]	buff	 	- The structure with all buffers to be filled,
 *			   	  and agents to be created.
 * @param[in]	params	 	- Provide the parameters for buffer's
 *				  initialization.
  * */
void initiate( HBBuffers_t *const buff, const Parameters_t *const params )
{
	size_t bug_locus;


	/* Set seed for GLibs's global random number generator. */
	g_random_set_seed( params->seed );


	/* Set vectors to zero. */
	memset( buff->swarm_map, RESET, params->world_size * sizeof( unsigned int ) );
	memset( buff->heat_map[ 0 ], RESET, params->world_size * sizeof( float ) );
	memset( buff->heat_map[ 1 ], RESET, params->world_size * sizeof( float ) );
	memset( buff->unhappiness, RESET, params->bugs_number * sizeof( float ) );


	/* Initiate swarm (that is, bug population) and swarm map. */
//	printf( "    > Generating bugs...\n" );

	/* Choose 'bugs_number' number of random world positions. */
	for (size_t bug_id = 0; bug_id < params->bugs_number; bug_id++)
	{
		/* Find a new free position. */
		do {
			bug_locus = (size_t) g_random_int_range( 0, params->world_size );	/* Interval [0..world_size[ as it should! */
		} while (HAS_BUG( buff->swarm_map[ bug_locus ] ));

		/* Free position found, create new bug in the swarm_map. */
		NEW_BUG_IN( buff->swarm_map[ bug_locus ] );


		/* Create the bug's data in the swarm. */

		SET_BUG_LOCAL( buff->swarm[ bug_id ].locus, bug_locus );

		SET_BUG_IDEAL_TEMPERATURE( buff->swarm[ bug_id ].ideal_temperature,
			g_random_int_range( params->bugs_temperature_min_ideal,
					params->bugs_temperature_max_ideal ) );

		SET_BUG_OUTPUT_HEAT( buff->swarm[ bug_id ].output_heat,
			g_random_int_range( params->bugs_heat_min_output ,
					params->bugs_heat_max_output ) );

		/* Update initial bug unhappiness as abs(ideal_temperature -
		 * temperature).
		 * Since world temperature is initial zero, it turns out that
		 * initial unhappiness  = ideal temperature.
		 * */

		buff->unhappiness[ bug_id ] = (float) buff->swarm[ bug_id ].ideal_temperature;
	}

	return;
}



/**
 * Initiate the world and create agents.
 *
 * @param[in]	heat_map	- The buffer with temperature data.
 * @param[out]	heat_buffer	- Buffer to be used as heat computation holder
 *				  (i.e. a double buffer).
 * @param[in]	params	 	- Provide buffer's and simulation parameters.
 * */
void comp_world_heat_v1( const float *const heat_map,
				float *const heat_buffer,
				const Parameters_t *const params )
{
	size_t ln, ls, ce, cw;	/* Line at North/South. Column at East/West. */
	size_t lc, cc;		/* Line at Center. Column at Center.	     */
	size_t pos;


/* Debug */
/*
for (size_t wpos = 0; wpos < params->world_size; wpos++)
	heat_map[wpos] = wpos + 1;
*/



	/*
	 * The vectors 'heat_map' and 'heat_buffer' grow from left to right
	 * (west to east) and then from bottom to top (south to north).
	 * Follow heat_buffer vector, compute Center Line and Center Column.
	 * Then compute the lines North/South of 'lc' and Columns East/West of
	 * 'cc'.
	 * Find heat contribution from all (lc, cc) neighbours present in
	 * 'heat_map'.
	 * */
	for (size_t wpos = 0; wpos < params->world_size; wpos++)
	{
		lc = (size_t) wpos / params->world_width;
		cc = (size_t) wpos - lc * params->world_width;

		ln = (lc + 1) % params->world_height;
		ls = (lc + params->world_height - 1) % params->world_height;
		ce = (cc + 1) % params->world_width;
		cw = (cc + params->world_width - 1) % params->world_width;


		/** Compute diffusion. */

		/* NW */
		pos = ln * params->world_width + cw;
		heat_buffer[ wpos ]  = heat_map[ pos ];	/* STORE the heat. */

		/* N */
		pos = ln * params->world_width + cc;
		heat_buffer[ wpos ] += heat_map[ pos ];	/* Accumulate the heat. */

		/* NE */
		pos = ln * params->world_width + ce;
		heat_buffer[ wpos ] += heat_map[ pos ];	/* Accumulate the heat. */

		/* W */
		pos = lc * params->world_width + cw;
		heat_buffer[ wpos ] += heat_map[ pos ];	/* Accumulate the heat. */

		/* E */
		pos = lc * params->world_width + ce;
		heat_buffer[ wpos ] += heat_map[ pos ];	/* Accumulate the heat. */

		/* SW */
		pos = ls * params->world_width + cw;
		heat_buffer[ wpos ] += heat_map[ pos ];	/* Accumulate the heat. */

		/* S */
		pos = ls * params->world_width + cc;
		heat_buffer[ wpos ] += heat_map[ pos ];	/* Accumulate the heat. */

		/* SE */
		pos = ls * params->world_width + ce;
		heat_buffer[ wpos ] += heat_map[ pos ];	/* Accumulate the heat. */

		/* Get the 8th of the diffusion percentage of all neighbour cells. */
		heat_buffer[ wpos ] = heat_buffer[ wpos ] * params->world_diffusion_rate / 8;

		/* Add cell's remaining heat. */
		heat_buffer[ wpos ] += heat_map[ wpos ] * (1 - params->world_diffusion_rate);


		/** Compute evaporation. */

		heat_buffer[ wpos ] = heat_buffer[ wpos ] * (1 - params->world_evaporation_rate);
	}

/* Debug */
/*
	size_t lin = params->world_height - 1;
	do {
		for (size_t col = 0; col < params->world_width; col++)
		{
			pos = lin * params->world_width + col;

			printf( "%05.2f ", heat_buffer[pos] );
		}
		printf( "\n" );
	} while (lin-- > 0);

	exit(0);
*/

	return;
}



void comp_world_heat_v2( const float *const heat_map,
				float *const heat_buffer,
				const Parameters_t *const params )
{
	size_t i, j, c;


/* Debug */
/*
for (size_t wpos = 0; wpos < params->world_size; wpos++)
	heat_map[wpos] = wpos + 1;
*/


	/** NORTH*/
		/* Compute north cells contribution. */

	i = 0;

	j = params->world_width;

	while (i < params->world_size - params->world_width)
		heat_buffer[i++] = heat_map[j++];

	j = 0;

	while (i < params->world_size)
		heat_buffer[i++] = heat_map[j++];


	/** SOUTH */
		/* Compute south cells contributions. */

	i = 0;

	j = params->world_size - params->world_width;

	while (i < params->world_width)
		heat_buffer[i++] += heat_map[j++];

	j = 0;

	while (i < params->world_size)
		heat_buffer[i++] += heat_map[j++];


	/** EAST */
		/* Compute east cells contributions. */

	i = 0;

	while (i < params->world_size)
	{
		j = i + 1;
		c = 0;

		while (c < params->world_width - 1)
		{
			heat_buffer[i++] += heat_map[j++];
			c++;
		}

		j = i + 1 - params->world_width;

		heat_buffer[i++] += heat_map[j];
	}


	/** WEST */
		/* Compute west cells contribution. */

	i = 0;

	while (i < params->world_size)
	{
		j = i + params->world_width - 1;

		heat_buffer[i++] += heat_map[j];

		j = i - 1;
		c = 0;

		while (c < params-> world_width - 1)
		{
			heat_buffer[i++] += heat_map[j++];
			c++;
		}
	}


	/** NORTHEAST */
		/* Compute northeast cells contributions. */

	i = 0;

	while (i < params->world_size - params->world_width)
	{
		j = i + params->world_width + 1;
		c = 0;

		while (c < params->world_width - 1)
		{
			heat_buffer[i++] += heat_map[j++];
			c++;
		}

		j = i + 1;

		heat_buffer[i++] += heat_map[j];
	}

	j = 1;
	c = 0;

	while (c < params->world_width - 1)
	{
		heat_buffer[i++] += heat_map[j++];
		c++;
	}

	j = 0;

	heat_buffer[i] += heat_map[j];


	/** NORTHWEST */
		/* Compute northwest cells contributions. */

	i = 0;

	while (i < params->world_size - params->world_width)
	{
		j = i + (params->world_width << 1) - 1;

		heat_buffer[i++] += heat_map[j];

		j = i + params->world_width - 1;
		c = 0;

		while (c < params->world_width - 1)
		{
			heat_buffer[i++] += heat_map[j++];
			c++;
		}
	}

	j = params->world_width - 1;

	heat_buffer[i++] += heat_map[j];

	j = 0;
	c = 0;

	while (c < params->world_width - 1)
	{
		heat_buffer[i++] += heat_map[j++];
		c++;
	}


	/** SOUTHEST */
		/* Compute southeast cells contributions. */

	i = 0;
	j = params->world_size - params->world_width + 1;
	c = 0;

	while (c < params->world_width - 1)
	{
		heat_buffer[i++] += heat_map[j++];
		c++;
	}

	j = params->world_size - params->world_width;

	heat_buffer[i++] += heat_map[j];

	while (i < params->world_size)
	{
		j = i - params->world_width + 1;
		c = 0;

		while (c < params->world_width - 1)
		{
			heat_buffer[i++] += heat_map[j++];
			c++;
		}

		j = i + 1 - (params->world_width << 1);

		heat_buffer[i++] += heat_map[j++];
	}


	/** SOUTHWEST */
		/* Compute southwest cells contribution. */

	i = 0;

	j = params->world_size - 1;

	heat_buffer[i++] += heat_map[j];

	j = params->world_size - params->world_width + 1;
	c = 0;

	while (c < params->world_width - 1)
	{
		heat_buffer[i++] += heat_map[j++];
		c++;
	}

	while (i < params->world_size)
	{
		j = i - 1;

		heat_buffer[i++] += heat_map[j];

		j = i - params->world_width - 1;
		c = 0;

		while (c < params->world_width - 1)
		{
			heat_buffer[i++] += heat_map[j++];
			c++;
		}
	}


	for (i = 0; i < params->world_size; i++)
	{
		/* Get the 8th of the diffusion percentage of all neighbour cells. */
		heat_buffer[ i ] = heat_buffer[ i ] * params->world_diffusion_rate / 8;

		/* Add cell's remaining heat. */
		heat_buffer[ i ] += heat_map[ i ] * (1 - params->world_diffusion_rate);


		/** Compute evaporation. */

		heat_buffer[ i ] = heat_buffer[ i ] * (1 - params->world_evaporation_rate);
	}


/* Debug */
/*
	size_t pos;
	size_t lin = params->world_height - 1;
	do {
		for (size_t col = 0; col < params->world_width; col++)
		{
			pos = lin * params->world_width + col;

			printf( "%05.2f ", heat_buffer[pos] );
		}
		printf( "\n" );
	} while (lin-- > 0);

	exit(0);
*/

	return;
}



unsigned int best_free_neighbour( const int todo, const float *const heat_map,
	const unsigned int *const swarm_map, const Parameters_t *const params,
	const size_t bug_locus)
{
	/* Agent position into the world / 2D position. */

	/* Central Row */
	const size_t rc = bug_locus / params->world_width;
	/* Contral Column */
	const size_t cc = bug_locus % params->world_width;

	/* Neighbouring rows and columns. */

	/* Row at North. */
	const size_t rn = (rc + 1) % params->world_height;
	/* Row at South. */
	const size_t rs = (rc + params->world_height - 1) % params->world_height;
	/* Column at East. */
	const size_t ce = (cc + 1) % params->world_width;
	/* Columns at West. */
	const size_t cw = (cc + params->world_width - 1) % params->world_width;

	struct {
		size_t pos;
		float heat;
	} best, neighbour[ NUM_NEIGHBOURS ];

	unsigned int NEIGHBOUR_IDX[ NUM_NEIGHBOURS ] = {SW, S, SE, W, E, NW, N, NE};


	/* Compute back the vector positions. Used on both, best */
	/* location and random location.                         */

	/* SW neighbour position. */
	neighbour[ SW ].pos = rs * params->world_width + cw;
	/* S neighbour position.  */
	neighbour[ S  ].pos = rs * params->world_width + cc;
	/* SE neighbour position. */
	neighbour[ SE ].pos = rs * params->world_width + ce;
	/* W neighbour position.  */
	neighbour[ W  ].pos = rc * params->world_width + cw;
	/* E neighbour position.  */
	neighbour[ E  ].pos = rc * params->world_width + ce;
	/* NW neighbour position. */
	neighbour[ NW ].pos = rn * params->world_width + cw;
	/* N neighbour position.  */
	neighbour[ N  ].pos = rn * params->world_width + cc;
	/* NE neighbour position. */
	neighbour[ NE ].pos = rn * params->world_width + ce;

	if (todo != FIND_ANY_FREE)
	{
		/* Fetch temperature of all neighbours. */

		/* SW neighbour temperature. */
		neighbour[ SW ].heat = heat_map[ neighbour[ SW ].pos ];
		/* S neighbour temperature. */
		neighbour[ S  ].heat = heat_map[ neighbour[ S  ].pos ];
		/* SE neighbour temperature. */
		neighbour[ SE ].heat = heat_map[ neighbour[ SE ].pos ];
		/* W neighbour temperature. */
		neighbour[ W  ].heat = heat_map[ neighbour[ W  ].pos ];
		/* E neighbour temperature. */
		neighbour[ E  ].heat = heat_map[ neighbour[ E  ].pos ];
		/* NW neighbour temperature. */
		neighbour[ NW ].heat = heat_map[ neighbour[ NW ].pos ];
		/* N neighbour temperature. */
		neighbour[ N  ].heat = heat_map[ neighbour[ N  ].pos ];
		/* NE neighbour temperature. */
		neighbour[ NE ].heat = heat_map[ neighbour[ NE ].pos ];

		/* Actual bug location is the best location, until otherwise. */
		best.pos = bug_locus;			/* Bug position. */
		best.heat = heat_map[ best.pos ];	/* Temperature at */
							/* bug position. */

		/* Loop unroll. */
		if (todo == FIND_MAX_TEMPERATURE)
		{
			if (neighbour[ SW ].heat > best.heat) best = neighbour[ SW ];
			if (neighbour[ S  ].heat > best.heat) best = neighbour[ S  ];
			if (neighbour[ SE ].heat > best.heat) best = neighbour[ SE ];
			if (neighbour[ W  ].heat > best.heat) best = neighbour[ W  ];
			if (neighbour[ E  ].heat > best.heat) best = neighbour[ E  ];
			if (neighbour[ NW ].heat > best.heat) best = neighbour[ NW ];
			if (neighbour[ N  ].heat > best.heat) best = neighbour[ N  ];
			if (neighbour[ NE ].heat > best.heat) best = neighbour[ NE ];
		}
		else	/* todo == FIND_MIN_TEMPERATURE */
		{
			if (neighbour[ SW ].heat < best.heat) best = neighbour[ SW ];
			if (neighbour[ S  ].heat < best.heat) best = neighbour[ S  ];
			if (neighbour[ SE ].heat < best.heat) best = neighbour[ SE ];
			if (neighbour[ W  ].heat < best.heat) best = neighbour[ W  ];
			if (neighbour[ E  ].heat < best.heat) best = neighbour[ E  ];
			if (neighbour[ NW ].heat < best.heat) best = neighbour[ NW ];
			if (neighbour[ N  ].heat < best.heat) best = neighbour[ N  ];
			if (neighbour[ NE ].heat < best.heat) best = neighbour[ NE ];
		}

		/*
		   Return if the bug is already in the best local or if the
		   best local is bug free.
		 * */
		if ((best.pos == bug_locus) || HAS_NO_BUG( swarm_map[ best.pos ] ))
			return best.pos;

	} /* end_if (todo != GOTO_ANY_FREE) */


	/*
	 * Here, there is a random moving chance, or best neighbour is not
	 * free. Try to find any available free place.
	 * */

	/*
	 * Fisher-Yates shuffle algorithm to shuffle an index vector so we can
	 * pick a random free neighbour by checking each index until find a
	 * first free.
	 * */
	for (size_t i = 0; i < NUM_NEIGHBOURS; i++)
	{
		size_t rnd_i = (size_t) g_random_int_range( i, NUM_NEIGHBOURS );

		if (rnd_i == i) continue;	/* Next shuffle. */

		/* Warning, this macro is using C99 extension. */
		SWAP( NEIGHBOUR_IDX[ i ], NEIGHBOUR_IDX[ rnd_i ] );
	}

	/* Loop unroll. */

	/* Find a first free neighbour. Index over the 8 neighbours. */
	best.pos = neighbour[ NEIGHBOUR_IDX[0] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	best.pos = neighbour[ NEIGHBOUR_IDX[1] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	best.pos = neighbour[ NEIGHBOUR_IDX[2] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	best.pos = neighbour[ NEIGHBOUR_IDX[3] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	best.pos = neighbour[ NEIGHBOUR_IDX[4] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	best.pos = neighbour[ NEIGHBOUR_IDX[5] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	best.pos = neighbour[ NEIGHBOUR_IDX[6] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	best.pos = neighbour[ NEIGHBOUR_IDX[7] ].pos;
	if (HAS_NO_BUG( swarm_map[ best.pos ] )) return best.pos;

	return bug_locus;	/* There is no free neighbour. */
}



void bug_step( bug_t *const swarm, unsigned int *const swarm_map,
			float *const heat_map, float *const unhappiness,
			const Parameters_t *const params )
{
	static size_t *ids = NULL;		/* Bugs id. */

	size_t bug_locus, bug_new_locus;
	int todo;


	/* One time vector initialization. Used to shuffle bugs. */
	if (!ids)
	{
		ids = (size_t *) malloc( params->bugs_number * sizeof( size_t ) );

		for (size_t idx = 0; idx < params->bugs_number; idx++)
			ids[ idx ] = idx;
	}

	/*
	 * Fisher-Yates shuffle algorithm.
	 * Use the vector to add randomness to the order bugs are selected
	 * for moving, preventing the same bug to always get the chance	to
	 * pick best location first.
	 * It should be more efficient to shuffle an integer bug indexer vector
	 * than shuffle the bugs vector, since every bug is a 3 integer
	 * structure.
	 * */

	/* Shuffle bugs indexer vector. */
	for (size_t idx = 0; idx < params->bugs_number; idx++)
	{
		/* The chance of j == i CANNOT be excluded because keeping the	*/
		/* value in the same position generates also a valid sequence.	*/
		size_t rnd_idx = (size_t) g_random_int_range( idx, params->bugs_number );

		if (rnd_idx == idx) continue;	/* Next shuffle.	*/

		/* Warning, this macro is using C99 extension. */
		SWAP( ids[ idx ], ids[ rnd_idx ] );
	}


	/* For each bug, indexed by bug_ids[ idx ]. */

	#define BUG ids[idx]

	for (size_t idx = 0; idx < params->bugs_number; idx++)
	{
		bug_locus = swarm[ BUG ].locus;

		/* Compute bug unhappiness, before trying to move. */
		unhappiness[ BUG ] =
			fabs( (float) swarm[ BUG ].ideal_temperature - heat_map[ bug_locus ] );

		/*
		 * Usually compare equality of floats is absurd. Netlogo
		 * wrapps the code with: if (unhappiness > 0) { do stuff... }
		 * I decided to unwrap, terminating as soon as possible using
		 * if (bug_unhappiness = 0.0f) {...}  and continue the
		 * remaining code as a fall out case for (bug_unhappiness >
		 * 0.0f). After all, this is semantically equivalent to
		 * netlogo version.
		 * */
		if (unhappiness[ BUG ] == 0.0f)
		{
			 /* Bug hasn't move, we don't need to update swarm. */
			 heat_map[ bug_locus ] += swarm[ BUG ].output_heat;
			 continue;	/* Next bug. */
		}

		/* Arriving here, means (unhappiness > 0.0f) */

		/*
			Find the best place for the bug to go.
			In order to implement netlogo approach, that is
			(in netlogo order), to compute:
			1) random-move-chance,
			2) best-patch for (temp <  ideal_temp) (when bug is in COLD),
			3) best-patch for (temp >= ideal_temp) (when bug is in HOT),
			we use the C conditional function twice.

			A variable is used to hold what to do, (1), (2) or (3).
			However the order must be reversed since (1), when
			happen, takes precedence over (2) or (3), whatever
			(2) XOR (3) is true or not.
		*/

		todo = (heat_map[ bug_locus ] < swarm[ BUG ].ideal_temperature)
				? FIND_MAX_TEMPERATURE : FIND_MIN_TEMPERATURE;

		todo = (g_random_double_range(0, 100) < params->bugs_random_move_chance)
				? FIND_ANY_FREE : todo;


		bug_new_locus = best_free_neighbour( todo, heat_map, swarm_map,
						params, bug_locus);

		/*
			Since the execution line is serial, 'bug_new_locus' is
			garantee to be: 1) the same as 'bug_locus', or 2) a new
			free location.
			It is not necessary to check again the free position,
			so we can update temperature at 'bug_new_locus',
			and only then check if the bug move or stay at the same
			'bug_locus' position.
		*/
		heat_map[ bug_new_locus ] += swarm[ BUG ].output_heat;


		/* If bug's current location is already the best one... */
		if (bug_new_locus == bug_locus)
		{
			/* Bug hasn't move, we don't need to update swarm's. */
			continue;	/* Next bug. */
		}

		/* Otherwise, move the bug to his new 'best' location. */

		/* Set new bug location in the "swarm". */
		swarm[ BUG ].locus = bug_new_locus; /* Update swarm. */

		/* Move the bug to the new location in the "swarm_map". */
		swarm_map[ bug_locus ] = A_EMPTY_CELL;
		swarm_map[ bug_new_locus ] = A_BUG;
	} /* END_for_each_bug */
} /* end bug_step(...) */



/**
 * Initiate the world and create agents.
 * */
void simulate( HBBuffers_t *const buff, const Parameters_t *const params,
			FILE *hbResultFile, GError **err )
{
	GError *err_simulate = NULL;

	float unhapp_average;
	hbBufferSelect_t bufsel;
	size_t iter_counter;


	/** Get unhappiness. */
	unhapp_average = average( buff->unhappiness, params->bugs_number );

	/* Output result to file. */
	fprintf( hbResultFile, "%.17g\n", unhapp_average );


	iter_counter = 0;


	bufsel.main = 0;    /* On first step, main buffer has index 0.	*/
	bufsel.secd = 1;    /* on first step, secondary buffer has index 1. */


	/*******************************/
	/**      SIMULATION LOOP      **/
	/*******************************/

	while ( (iter_counter < params-> numIterations)
		|| (params->numIterations == 0) )
	{
		/** Compute world heat, diffusion followed by evaporation. */
		/* Use 'bufsel' to point the current buffer. */
		comp_world_heat_v2( buff->heat_map[bufsel.main],
					buff->heat_map[bufsel.secd],
					params );

		/** Perform bug step. */
		/* Use 'bufsel' to point the correct buffer. */
		bug_step( buff->swarm, buff->swarm_map,
				buff->heat_map[bufsel.secd], buff->unhappiness,
				params );

		/** Get unhappiness. */
		unhapp_average = average( buff->unhappiness, params->bugs_number );

		/* Output result to file. */
		fprintf( hbResultFile, "%.17g\n", unhapp_average );

		/** Prepare next iteration. */

		iter_counter++;

		/* Warning, this macro is using C99 extension. */
		SWAP( bufsel.main, bufsel.secd );
	}
}




int main( int argc, char *argv[] )
{
	FILE *hbResultFile = NULL;
	GError *err_main = NULL;	/* Error reporting object, from Glib. */

	Parameters_t params;		/* Simulation parameters. */

	HBBuffers_t buff = { NULL, NULL, { NULL, NULL }, NULL };	/* Buffers used for simulation. */



	getSimulParameters( &params, argc, argv, &err_main );
	hb_if_err_goto( err_main, error_handler );

	setupBuffers( &buff, &params, &err_main );
	hb_if_err_goto( err_main, error_handler );


	/* Open output file for results. */
	hbResultFile = fopen(params.output_filename, "w+");	/* Open file overwrite. */
	hb_if_err_create_goto( err_main, HB_ERROR,
		hbResultFile == NULL, HB_UNABLE_OPEN_FILE, error_handler,
		"Could not open output file." );


	/* Initiate. */
	initiate( &buff, &params );

	/* Simulate */
	simulate( &buff, &params, hbResultFile, &err_main );


//	printf( "End...\n\n" );

	/* TODO: Profiling. */


	goto clean_all;


error_handler:

	/* Handle error. */
	fprintf( stderr, "Error: %s\n", err_main->message );
	g_error_free( err_main );


clean_all:

	fclose( hbResultFile );

	if (buff.unhappiness) free( buff.unhappiness );
	if (buff.heat_map[ 1 ])	free( buff.heat_map[ 1 ] );
	if (buff.heat_map[ 0 ]) free( buff.heat_map[ 0 ] );
	if (buff.swarm_map) free( buff.swarm_map );
	if (buff.swarm) free( buff.swarm );

	// if (err_main) g_error_free( err_main );


	return OKI_DOKI;
}
