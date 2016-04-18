void bug_step( mundus_t * const world, populatio_t * const population, const real_t random_move_chance )
{
	static unsigned int *bug_ids = NULL;		/* To add randomness in the		*/
											/* order bugs are moved.		*/

	bug_t bug;			/* Current bug, as seen by the processor. */
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
		/* Get the current bug from the swarm.  */
		bug = population->swarm[ bug_ids[idx] ];

		/* The bug unhapiness is the magnitude (absolute value) of the		*/
		/* difference between bug ideal temperature and the temperature at	*/
		/* bug current location.											*/
		/* NetLogo computes unhapiness before bugs moves away from current	*/
		/* position for each step, so first two values should be always the	*/
		/* same; unhapiness from setup time and unhapiness of first step!	*/
		bug.unhappiness =
			(real_t) fabs( bug.ideal_temperature - world->temperature_map[ bug.local.lin ][ bug.local.col ] );

		if (bug.unhappiness > 0.0)
		{
			/** REPORT_PATCH: Find new target location / patch, for the bug. */

			if ((real_t) g_random_double_range(0, 100) < random_move_chance)
			{
				/* Find a free neighbour. */
				bug_target = bestFreeNeighbour( RANDOM_CHANCE, world, &bug.local );
			}
			else if (world->temperature_map[ bug.local.lin ][ bug.local.col ] < bug.ideal_temperature)
			{
				/* Find the hottest location next to me; also take my current location into consideration. */
				bug_target = bestFreeNeighbour( MAX_TEMPERATURE, world, &bug.local );
			}
			else
			{
				/* Find the coolest location next to me; also take my current location into consideration. */
				bug_target = bestFreeNeighbour( MIN_TEMPERATURE, world, &bug.local );
			}

			/** BUG_MOVE: Move the bug to target, if it is not already in the target / current location. */

			if (bug_target.lin != bug.local.lin || bug_target.col != bug.local.col)
			{
				world->swarm_map[ bug.local.lin ][ bug.local.col ] = EMPTY;
				world->swarm_map[ bug_target.lin ][ bug_target.col ] = BUG;

				/* The bug moved, so set the new bug local. */
				bug.local = bug_target;

				/* Update the bug in the swarm, this updates both, bug location and bug unhappiness. */
				population->swarm[ bug_ids[idx] = bug;

			}
		} /* END if random move chance */


		/* UPDATE: Recompute the temperature for the new bug location. Here, 'bug.local' may be the old bug position or may be the new bug position. */
		world->temperature_map[ bug.local.lin ][ bug.local.col ] += bug.output_heat;
	} /* END for (each_bug) */
} /* Bug_step done... */

