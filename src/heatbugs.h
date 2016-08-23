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

#ifndef __HEATBUGS_CPU_H_
#define __HEATBUGS_CPU_H_


/* Comment next to use floats. Default is doubles. */
#define USE_PRECISION_MATH




#ifdef USE_PRECISION_MATH
#define real_t double
#else
#define real_t float
#endif


#define VSIZE( v ) ( sizeof( v ) / sizeof( v[0] ) )


/* Input data used for simulation. */
typedef struct simulation {
	real_t bugs_min_ideal_temperature;
	real_t bugs_max_ideal_temperature;
	real_t bugs_min_output_heat;
	real_t bugs_max_output_heat;
	real_t bugs_random_move_chance;		/* Chance a bug will move [0..100].       */
	real_t wrl_diffusion_rate;		/* % temperature to adjacent cells [0..1] */
	real_t wrl_evaporation_rate;		/* % temperature's loss to 'ether' [0..1] */
	unsigned int world_height;
	unsigned int world_width;
	unsigned int number_of_bugs;		/* The number of bugs in the world.       */
	unsigned int numIterations;		/* Iterations to stop. (0 = no stop).     */
} simulatio_t;

/* Coordinates used for location. */
typedef struct locus {
	unsigned int lin;
	unsigned int col;
} locus_t;

/* Related with 2D dimensions. */
typedef struct dimensio {
	unsigned int height;
	unsigned int width;
} dimensio_t;

/* Bug specific data, different for each bug. */
/* Initialized from simulation parameters.    */
typedef struct bug {
	real_t ideal_temperature;
	real_t output_heat;
	real_t unhappiness;
	locus_t local;
} bug_t;


typedef struct populatio {
	unsigned int number_of_bugs;
	bug_t *swarm;				/* Bugs vector.	The swarm...	*/
} populatio_t;


typedef struct mundus {
	dimensio_t dimension;			/* World Dimensions.	*/
	real_t **temperature_map;		/* Temperature Array.	*/
	char **swarm_map;			/* Bugs location Indicator Array.	*/
} mundus_t;


typedef enum neighbours {SW, S, SE, W, E, NW, N, NE} neighbours_t;

typedef enum movselect {RANDOM_CHANCE, MIN_TEMPERATURE, MAX_TEMPERATURE} movselect_t;


/* Create float(FArray) and char(CArray) arrays. */
real_t **newFArray( unsigned int height, unsigned int width );
char **newCArray( unsigned int height, unsigned int width );

/* Set float array to zero. */
void resetFArray( real_t **array, unsigned int height, unsigned int width );
void resetCArray( char **array, unsigned int height, unsigned int width );

real_t mean_unhappiness( populatio_t * const population );

void setup( mundus_t * const world, populatio_t * const population, simulatio_t * const simul );

locus_t bestFreeNeighbour( movselect_t const moveselect, mundus_t const * const world, locus_t const * const local );

void world_diffuse( mundus_t * const world, const real_t diffusion_rate );
void world_evaporate( mundus_t * const world, const real_t evaporation_rate );
void bug_step( mundus_t * const world, populatio_t * const population, const real_t random_move_chance );

#endif
