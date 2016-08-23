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


#include <stdio.h>			/* printf(...)	*/
#include <stdlib.h>			/* exit(...)	*/
#include <math.h>			/* fabs(...)	*/
#include <string.h>

#include "glib.h"
//#include "keyboard.h"
#include "heatbugs.h"


#define OKI_DOKI 0
#define NOT_DOKI -1

#define BUG '@'
#define EMPTY '-'



const char version[] = "Heatbugs simulation for CPU (serial processing) v2.8 with Glib-2.0 randoms.";


/**
 * Return an array of floats|doubles.
 * */
real_t **newFArray( unsigned int height, unsigned int width )
{
	real_t **array;

	array = (real_t **) malloc( height * sizeof(real_t *) );

	for (unsigned int i = 0; i < height; i++) {
		array[i] = (real_t *) malloc( width * sizeof(real_t) );
		memset( array[i], 0, width * sizeof(real_t) );
	}

	return array;
}

/**
 * Return an array of chars.
 * */
char **newCArray( unsigned int height, unsigned int width )
{
	char **array;

	array = (char **) malloc( height * sizeof(char *) );

	for (unsigned int i = 0; i < height; i++) {
		array[i] = (char *) malloc( width * sizeof(char) );
		memset( array[i], EMPTY, width * sizeof(char) );
	}

	return array;
}

/*
 * Reset float|double array to zero.
 * */
void resetFArray( real_t **array, unsigned int height, unsigned int width )
{
	for (unsigned int i = 0; i < height; i++)
		memset( array[i], 0, width * sizeof(real_t) );
}

/*
 * Reset char array to space.
 * */
void resetCArray( char **array, unsigned int height, unsigned int width )
{
	for (unsigned int i = 0; i < height; i++)
		memset( array[i], EMPTY, width * sizeof(char) );
}


/*
 * Compute bug's unhapiness average.
 * */
real_t mean_unhappiness( populatio_t * const population )
{
	real_t total_unhappiness = 0.0;


	for (unsigned int bug = 0; bug < population->number_of_bugs; bug++)
	{
		total_unhappiness += population->swarm[bug].unhappiness;
	}

	return (total_unhappiness / population->number_of_bugs);
}


void setup( mundus_t * const world, populatio_t * const population, simulatio_t * const simul )
{
	unsigned int l, c;

	/* Issue warning for number of bugs near | above available world space.	*/

	unsigned int worldSize = simul->world_height * simul->world_width;

	if (simul->number_of_bugs == 0) {
		fprintf(stderr, "Error: There are no bugs.\n\n");
		exit(NOT_DOKI);
	}

	if (simul->number_of_bugs >= worldSize) {
		fprintf(stderr, "Error: Bugs number exceed available world slots.\n\n");
		exit(NOT_DOKI);
	}

	if (simul->number_of_bugs >= 0.8 * worldSize) {
		fprintf(stderr, "Warning: Bugs number close available world slots.\n");
		fprintf(stderr, "         Generation may be slow due to collision.\n\n");
	}

	world->dimension.height = simul->world_height;
	world->dimension.width  = simul->world_width;

	/* Temperature map. [Initial slots state is zero] */
	world->temperature_map = newFArray( world->dimension.height, world->dimension.width );

	/* Bug's location indication in the world. */
	world->swarm_map = newCArray( world->dimension.height, world->dimension.width );


	/* Initialize bug population. */

	population->number_of_bugs = simul->number_of_bugs;

	/* Bugs vector. */
	population->swarm = (bug_t *) malloc( population->number_of_bugs * sizeof(bug_t) );

	/* Issue information "Generating bugs population." */
	printf("    > Generating bugs...\n");

	/* Choose a 'number_of_bugs' number of random world positions,	*/
	/* allow bugs to be initialized in them...			*/
	for (unsigned int bug = 0; bug < population->number_of_bugs; bug++)
	{
		do {
			l = (unsigned int) g_random_int_range( 0, world->dimension.height );	/* Interval [0..height[	as it should! */
			c = (unsigned int) g_random_int_range( 0, world->dimension.width );	/* Interval [0..width[ as it should!  */
		} while (world->swarm_map[l][c] == BUG);

		world->swarm_map[l][c] = BUG;	/* Set the bug location in the world. */

		/* Initialize the bug location in the swarm... */
		population->swarm[bug].local.lin = l;
		population->swarm[bug].local.col = c;

		/* Initialize other bug related stuff. */

		population->swarm[bug].ideal_temperature = (real_t) g_random_int_range( simul->bugs_min_ideal_temperature, simul->bugs_max_ideal_temperature );

		population->swarm[bug].output_heat = (real_t) g_random_int_range( simul->bugs_min_output_heat, simul->bugs_max_output_heat );

		population->swarm[bug].unhappiness = (real_t) fabs( population->swarm[bug].ideal_temperature - world->temperature_map[l][c] );
	} /* end_for_loop */
}


/* Check "Pointers and references" in the link:								*/
/* https://en.wikipedia.org/wiki/Const_%28computer_programming%29			*/

/* Return in 'target' a neighbour of 'local':								*/
/* - The best min/max neighbour temperature location if	bug free! 			*/
/* - Any available / bug free neighbour location, as NetLogo does.			*/
/* - 'local' if this is already the best temperature location or if there	*/
/*   isn't any bug free neighbour.											*/
locus_t bestFreeNeighbour( movselect_t const moveselect, mundus_t const * const world, locus_t const * const local )
{
	unsigned int north, south, west, east;	/* From current local position.	*/
	real_t bestTemperature;
	locus_t best_local, tmp_local;

	static unsigned int neighbours[] = {SW, S, SE, W, E, NW, N, NE};

	south = (local->lin + world->dimension.height - 1) % world->dimension.height;
	north = (local->lin + world->dimension.height + 1) % world->dimension.height;
	west = (local->col + world->dimension.width - 1) % world->dimension.width;
	east = (local->col + world->dimension.width + 1) % world->dimension.width;

	/* Find best temperature, check all neighbours. */
	if (moveselect != RANDOM_CHANCE)
	{
		/* Consider the 'local' temperature first. Try to get a location	*/
		/* with better temperature and if there isn't, then a simple bug	*/
		/* free location.													*/
		bestTemperature = world->temperature_map[ local->lin ][ local->col ];
		best_local = *local;

		/* Check all 'local' neighbour's temperature. */
		for (neighbours_t i = SW; i <= NE; i++)
		{
			switch (i)
			{
				case SW: case S: case SE:
					tmp_local.lin = south;
					break;
				case NW: case N: case NE:
					tmp_local.lin = north;
					break;
				default:
					tmp_local.lin = local->lin;
					break;
			} /* end_switch */

			switch (i)
			{
				case SW: case W: case NW:
					tmp_local.col = west;
					break;
				case SE: case E: case NE:
					tmp_local.col = east;
					break;
				default:
					tmp_local.col = local->col;
					break;
			} /* end_switch */

			if (moveselect == MAX_TEMPERATURE)
			{
				if (world->temperature_map[ tmp_local.lin ][ tmp_local.col ] > bestTemperature)
				{
					bestTemperature = world->temperature_map[ tmp_local.lin ][ tmp_local.col ];
					best_local = tmp_local;
				}
			}
			else /* moveselect is MIN_TEMPERATURE */
			{
				if (world->temperature_map[ tmp_local.lin ][ tmp_local.col ] < bestTemperature)
				{
					bestTemperature = world->temperature_map[ tmp_local.lin ][ tmp_local.col ];
					best_local = tmp_local;
				}
			}
		} /* end_for all neighbours */

		/* If 'best_loc' == 'local' already has the best temperature so, return...	*/
		if (best_local.lin == local->lin  &&  best_local.col == local->col) return *local;

		/* 'best_loc' is not 'local' (has changed)... and...	*/
		/* If 'target' location is bug free, return as well.	*/
		if (world->swarm_map[ best_local.lin ][ best_local.col ] != BUG) return best_local;
	}

	/* Otherwise try to find any bug free neighbour or do a RANDOM_MOVE_CHANCE.	*/

	/* Shuffle neighbour's vector so we can pick a random free neighbour. */
	for (unsigned int i = 0; i < (unsigned int) VSIZE(neighbours); i++)
	{
		unsigned int rnd = (unsigned int) g_random_int_range( i, VSIZE(neighbours) );

		if (rnd == i) continue;

		unsigned int tmp = neighbours[i];
		neighbours[i] = neighbours[rnd];
		neighbours[rnd] = tmp;
	}

	/* Try to find a bug free neighbour. */
	for (unsigned int i = 0; i < (unsigned int) VSIZE(neighbours); i++)
	{
		switch (neighbours[i])
		{
			case SW: case S: case SE:
				tmp_local.lin = south;
				break;
			case NW: case N: case NE:
				tmp_local.lin = north;
				break;
			default:
				tmp_local.lin = local->lin;
				break;
		} /* end_switch */

		switch (neighbours[i])
		{
			case SW: case W: case NW:
				tmp_local.col = west;
				break;
			case SE: case E: case NE:
				tmp_local.col = east;
				break;
			default:
				tmp_local.col = local->col;
				break;
		} /* end_switch */

		if (world->swarm_map[ tmp_local.lin ][ tmp_local.col ] != BUG)
		{
			return tmp_local;	/* A bug free neighbour. */
		}
	}

	return *local;	/* Cannot move to any location, return current location. */
}


/* Note that bug world is a zero based array in first quadrant,	so that		*/
/* (0,0) coordinate is at leftmost bottom position.							*/
/* The bugs world is circular, means it folds and join at the edges...		*/
void world_diffuse( mundus_t * const world, const real_t diffusion_rate )
{
	static real_t **buffer = NULL;

	real_t diffusion_leftover;
	unsigned int north, south, west, east;	/* From current (i,j) position. */
	unsigned int oldBuffer, newBuffer;		/* Buffer indexers.				*/

	/* A three line buffer used to hold diffusion calculations until it's   */
	/* safe to transfer back to array. Allocated just once.					*/
	if (!buffer) {
		buffer = newFArray( 3, world->dimension.width );
	}

	diffusion_leftover = 1.0 - diffusion_rate;	/* That remains in the cell. */

	/* Handle first special case, (line 0 (bottom) of temperature array). */
	north = 1;
	south = world->dimension.height - 1;

	for (unsigned int j = 0; j < world->dimension.width; j++)
	{
		west = (j + world->dimension.width - 1) % world->dimension.width;
		east = (j + 1) % world->dimension.width;

		/* Sum all neighbours heat to compute contribution. */

		/* Patch 1: North */
		buffer[0][j]  = world->temperature_map[ north ][ j ];
		/* Patch 2: South */
		buffer[0][j] += world->temperature_map[ south ][ j ];
		/* Patch 3: West */
		buffer[0][j] += world->temperature_map[ 0 ][ west ];
		/* Patch 4: East */
		buffer[0][j] += world->temperature_map[ 0 ][ east ];

		/* Patch 5: Northwest */
		buffer[0][j] += world->temperature_map[ north ][ west ];
		/* Patch 6: Northeast */
		buffer[0][j] += world->temperature_map[ north ][ east ];
		/* Patch 7: Southwest */
		buffer[0][j] += world->temperature_map[ south ][ west ];
		/* Patch 8: Southeast */
		buffer[0][j] += world->temperature_map[ south ][ east ];

		/* Each neighbour heat contribution is 1/8 of diffusion rate. */
		buffer[0][j] = buffer[0][j] * diffusion_rate / 8.0;

		/* Add the remaining heat left over. */
		buffer[0][j] = buffer[0][j]	+ world->temperature_map[0][j] * diffusion_leftover;
	}

	/* Handle remaining cases except the last one. */
	newBuffer = 1;		/* Required initializations so buffers are used for  */
	oldBuffer = 2;		/* each cycle in the right order, and then swapped.  */

	for (unsigned int i = 1; i < world->dimension.height - 1; i++)
	{
		north = i + 1;
		south = i - 1;

		for (unsigned int j = 0; j < world->dimension.width; j++)
		{
			west = (j + world->dimension.width - 1) % world->dimension.width;
			east = (j + 1) % world->dimension.width;

			/* Sum all neighbours heat to compute contribution. */

			/* Patch 1: North */
			buffer[newBuffer][j]  = world->temperature_map[ north ][ j ];
			/* Patch 2: South */
			buffer[newBuffer][j] += world->temperature_map[ south ][ j ];
			/* Patch 3: West */
			buffer[newBuffer][j] += world->temperature_map[ i ][ west ];
			/* Patch 4: East */
			buffer[newBuffer][j] += world->temperature_map[ i ][ east ];

			/* Patch 5: Northwest */
			buffer[newBuffer][j] += world->temperature_map[ north ][ west ];
			/* Patch 6: Northeast */
			buffer[newBuffer][j] += world->temperature_map[ north ][ east ];
			/* Patch 7: Southwest */
			buffer[newBuffer][j] += world->temperature_map[ south ][ west ];
			/* Patch 8: Southeast */
			buffer[newBuffer][j] += world->temperature_map[ south ][ east ];

			/* Each neighbour heat contribution is 1/8 of diffusion rate. */
			buffer[newBuffer][j] = buffer[newBuffer][j] * diffusion_rate / 8.0;

			/* Add the remaining heat left over. */
			buffer[newBuffer][j] = buffer[newBuffer][j]	+ world->temperature_map[i][j] * diffusion_leftover;
		}

		/* Replace previous buffer. */
		if (i >= 2) {
			real_t *tmp;

			tmp = world->temperature_map[south];
			world->temperature_map[south] = buffer[oldBuffer];
			buffer[oldBuffer] = tmp;
		}

		/* Swap buffer's indices. */
		do {
			unsigned int tmp = oldBuffer;
			oldBuffer = newBuffer;
			newBuffer = tmp;
		} while( 0 );
	} /* end_for i */

	/* Handle last case (last line). */
	/* This last case, prevents modular calculation for all lines. */
	north = 0;
	south = world->dimension.height - 2;

	/* Note: var newBuffer is being used here.	*/
	/* It comes correct from the last loop.		*/

	for ( unsigned int j = 0; j < world->dimension.width; j++ )
	{
		west = (j + world->dimension.width - 1) % world->dimension.width;
		east = (j + 1) % world->dimension.width;

		/* Patch 1: North */
		buffer[newBuffer][j]  = world->temperature_map[ north ][ j ];
		/* Patch 2: South */
		buffer[newBuffer][j] += world->temperature_map[ south ][ j ];
		/* Patch 3: West */
		buffer[newBuffer][j] +=
				world->temperature_map[ world->dimension.height - 1 ][ west ];
		/* Patch 4: East */
		buffer[newBuffer][j] +=
				world->temperature_map[ world->dimension.height - 1 ][ east ];

		/* Patch 5: Northwest */
		buffer[newBuffer][j] += world->temperature_map[ north ][ west ];
		/* Patch 6: Northeast */
		buffer[newBuffer][j] += world->temperature_map[ north ][ east ];
		/* Patch 7: Southwest */
		buffer[newBuffer][j] += world->temperature_map[ south ][ west ];
		/* Patch 8: Southeast */
		buffer[newBuffer][j] += world->temperature_map[ south ][ east ];

		/* Each neighbour heat contribution is 1/8 of diffusion rate. */
		buffer[newBuffer][j] = buffer[newBuffer][j] * diffusion_rate / 8.0;

		/* Add the remaining heat left over. */
		buffer[newBuffer][j] = buffer[newBuffer][j]	+ world->temperature_map[ world->dimension.height - 1 ][j] * diffusion_leftover;
	}

	/* Now update the three buffer entries (last entries) in the array.	*/
	/* Dont forget, it is a swap with the buffer.                       */
	do {
		real_t *tmp;

		tmp = world->temperature_map[0];
		world->temperature_map[0] = buffer[0];
		buffer[0] = tmp;

		tmp = world->temperature_map[world->dimension.height - 2];
		world->temperature_map[world->dimension.height - 2] = buffer[oldBuffer];
		buffer[oldBuffer] = tmp;

		tmp = world->temperature_map[world->dimension.height - 1];
		world->temperature_map[world->dimension.height - 1] = buffer[newBuffer];
		buffer[newBuffer] = tmp;
	} while(0);

}	/* Diffusion done... */


void world_evaporate( mundus_t * const world, const real_t evaporation_rate )
{
	real_t evaporation_leftover;

	evaporation_leftover = 1.0 - evaporation_rate;

	for (unsigned int i = 0; i < world->dimension.height; i++)
	{
		for (unsigned int j = 0; j < world->dimension.width; j++)
		{
			world->temperature_map[i][j] *= evaporation_leftover;
		}
	}
}	/* Evaporation done... */


void bug_step( mundus_t * const world, populatio_t * const population, const real_t random_move_chance )
{
	static unsigned int *bug_ids = NULL;		/* To add randomness in the		*/
											/* order bugs are moved.		*/

	locus_t bug_local;		/* Bug current location. */
	locus_t bug_target;		/* Bug target location.	 */

	/* One time vector initialization. */
	if (!bug_ids)
	{
		bug_ids = (unsigned int *) malloc( population->number_of_bugs * sizeof(unsigned int) );

		for (unsigned int id = 0; id < population->number_of_bugs; id++)
		{
			bug_ids[id] = id;
		}
	}

	/* Fisher-Yates shuffle algorithm. */
	/* Use the vector to add randomness to the order bugs are selected		*/
	/* for moving, preventing the same bug to always get the chance	to 		*/
	/* pick best location first.											*/

	/* Shuffle bugs. */
	for (unsigned int idx = 0; idx < population->number_of_bugs; idx++)
	{
		/* The chance of j == i CANNOT be excluded because keeping the		*/
		/* value in the same position generates also a valid sequence.		*/
		unsigned int rnd_idx = (unsigned int) g_random_int_range( idx, population->number_of_bugs );

		if (rnd_idx == idx) continue;	/* Next shuffle.	*/

		unsigned int tmp_bug_id = bug_ids[idx];
		bug_ids[idx] = bug_ids[rnd_idx];
		bug_ids[rnd_idx] = tmp_bug_id;
	}

	/* For each bug, indexed by idx... */
	for (unsigned int idx = 0; idx < population->number_of_bugs; idx++)
	{
		/* Get the bug current location. */
		bug_local = population->swarm[ bug_ids[idx] ].local;

		/* The bug unhapiness is the magnitude (absolute value) of the		*/
		/* difference between bug ideal temperature and the temperature at	*/
		/* bug current location.						*/
		/* NetLogo computes unhapiness before bugs moves away from current	*/
		/* position for each step, so first two values should be always the	*/
		/* same; unhapiness from setup time and unhapiness of first step!	*/
		population->swarm[ bug_ids[idx] ].unhappiness =
			(real_t) fabs( population->swarm[ bug_ids[idx] ].ideal_temperature - world->temperature_map[ bug_local.lin ][ bug_local.col ] );

		if (population->swarm[ bug_ids[idx] ].unhappiness > 0.0)
		{
			/** REPORT_PATCH: Find new target location / patch, for the bug. */

			if ((real_t) g_random_double_range(0, 100) < random_move_chance)
			{
				/* Find a free neighbour. */
				bug_target = bestFreeNeighbour( RANDOM_CHANCE, world, &bug_local );
			}
			else if (world->temperature_map[ bug_local.lin ][ bug_local.col ] < population->swarm[ bug_ids[idx] ].ideal_temperature)
			{
				/* Find the hottest location next to me; also take my current location into consideration. */
				bug_target = bestFreeNeighbour( MAX_TEMPERATURE, world, &bug_local );
			}
			else
			{
				/* Find the coolest location next to me; also take my current location into consideration. */
				bug_target = bestFreeNeighbour( MIN_TEMPERATURE, world, &bug_local );
			}

			/** BUG_MOVE: Move the bug to target, if it is not already in the target / current location. */

			if (bug_target.lin != bug_local.lin || bug_target.col != bug_local.col)
			{
				world->swarm_map[ bug_local.lin ][ bug_local.col ] = EMPTY;
				world->swarm_map[ bug_target.lin ][ bug_target.col ] = BUG;

				population->swarm[ bug_ids[idx] ].local = bug_target;

				/* The bug moved, so set the new bug local. */
				bug_local = bug_target;
			}
		} /* END if random move chance */


		/* UPDATE: Recompute the temperature in the new bug location. Here, 'bug_local' may be the old bug position or may be the new bug position. */
		world->temperature_map[ bug_local.lin ][ bug_local.col ] += population->swarm[ bug_ids[idx] ].output_heat;
	} /* END for (each_bug) */
} /* Bug_step done... */


/* Debug... */

void dataInput( mundus_t * const world )
{
	for (unsigned int j = 0; j < world->dimension.height; j++)
	{
		for (unsigned int i = 0; i < world->dimension.width; i++)
		{
			printf("[%d, %d] :> ", j, i);

			#ifdef USE_PRECISION_MATH
				scanf("%lf", &(world->temperature_map[j][i]));
			#else
				scanf("%f", &(world->temperature_map[j][i]));
			#endif
		}
	}
}

void displayFArray( mundus_t const * const world )
{
	for (unsigned int i = world->dimension.height; i > 0 ; i--)
	{
		for (unsigned int j = 0; j < world->dimension.width; j++)
		{
			printf("%f ", world->temperature_map[i - 1][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void displayCArray( mundus_t const * const world )
{
		for (unsigned int i = world->dimension.height; i > 0; i--)
		{
			for (unsigned int j = 0; j < world->dimension.width; j++)
			{
				printf("%c", world->swarm_map[i - 1][j]);
			}
			printf("\n");
		}
		printf("\n");
}


int main( int argc, char *argv[] )
{
	FILE *hbResultFile;

	simulatio_t simul;
	populatio_t population;
	mundus_t world;

	real_t average_unhappiness;


	/* TODO: Get simulation parameters. */
	/* setparam(); */
	simul.bugs_min_ideal_temperature = 20.0;	/* 10.0 */
	simul.bugs_max_ideal_temperature = 30.0;	/* 40.0 */
	simul.bugs_min_output_heat =  15.0;			/*  5.0 */
	simul.bugs_max_output_heat = 25.0;			/* 25.0 */
	simul.wrl_evaporation_rate = 0.01;			/*  1%. */
	simul.wrl_diffusion_rate  =  0.09;			/* 90%. */
	simul.bugs_random_move_chance = 0.0;		/* 10%. Valid:[0 .. 100] */

	simul.world_height = 5;
	simul.world_width =  5;

	simul.number_of_bugs = 15;					/* 100 Bugs in the world. */
	simul.numIterations =  10;				/* 0 = NonStop. */





	hbResultFile = fopen("../results/heatbugsC.csv", "w+");	/* Open file overwrite. */

	if (hbResultFile == NULL) {
		fprintf(stderr, "Error: Could not open output file.\n\n");
		/* TODO: flush stdin/stdout. */
		exit(-1);
	}

	printf("Using globally initiated Glib2.xx Random Generator...\n");

	/* Setup the environment and the bugs. */
	printf("Setup...\n");
	setup( &world, &population, &simul );

	/* Run simulation. */
	printf("Simulating...\n    > output to file.\n");

	// displayFArray( &world );
	// dataInput( &world );
	// displayFArray( &world );

	// displayCArray( &world );

	unsigned int iter_counter = 0;

	average_unhappiness = mean_unhappiness( &population );
	fprintf( hbResultFile, "%d, %.17g\n", iter_counter, average_unhappiness );

	/* If numIteration == 0 -> run non-stop. */
	while ( (simul.numIterations == 0) || (iter_counter < simul.numIterations) )
	{
		world_diffuse( &world, simul.wrl_diffusion_rate );
		world_evaporate( &world, simul.wrl_evaporation_rate );
		bug_step( &world, &population, simul.bugs_random_move_chance );

		/* dumpdata( &bugWorld ); */
		// displayCArray( &world );

		iter_counter++;

		average_unhappiness = mean_unhappiness( &population );
		fprintf( hbResultFile, "%d, %.17g\n", iter_counter, average_unhappiness );
	}

	// displayFArray( &world );
	// displayCArray( &world );
	printf("    > end...\n\n");

	fclose( hbResultFile );

	return OKI_DOKI;
}
