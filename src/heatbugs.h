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


#include  <stdio.h>


/**
* Error reporting macros from cf4ocl OpenCL library by Nuno Fachada, using
* GLib 2.0.
* Check:
*	https://fakenmc.github.io/cf4ocl
*	https://fakenmc.github.io/cf4ocl/docs/latest/ccl__common_8h_source.html
*/



#define NDEBUG




/* Define NDEBUG for specific debug. Does not compile in windows. */

#ifdef NDEBUG
	#define CCL_STRD G_STRFUNC
#else
	#define CCL_STRD G_STRLOC
#endif


/* Using C99 variadic macros. */

#define hb_if_err_create_goto( err, quark, error_condition, error_code, label, msg, ... ) \
	if (error_condition) { \
		g_debug( CCL_STRD ); \
		g_set_error( &(err), (quark), (error_code), (msg), ##__VA_ARGS__ ); \
		goto label; \
	}

#define hb_if_err_goto( err, label ) \
	if ((err) != NULL) { \
		g_debug( CCL_STRD ); \
		goto label; \
	}

#define hb_if_err_propagate_goto( err_dest, err_src, label ) \
	if ((err_src) != NULL) { \
		g_debug( CCL_STRD ); \
		g_propagate_error( err_dest, err_src ); \
		goto label; \
	}


/* Swap to values of indicated 'type'. Warning! C99 extension 'typeof'. */
#define SWAP( a, b ) { __typeof__(a) t = a; a = b; b = t; }


/** Heatbugs own error codes. **/
enum hb_error_codes {
	/** Successfull operation. */
	HB_SUCCESS = 0,
	/** Invalid parameters. */
	HB_INVALID_PARAMETER = -1,
	/** A command line option with missing argument. */
	HB_PARAM_ARG_MISSING = -2,
	/** Unknown option in the command line. */
	HB_PARAM_OPTION_UNKNOWN = -3,
	/** Unknown option characters in command line. */
	HB_PARAM_CHAR_UNKNOWN = -4,
	/** Weird error occurred while parsing parameter. */
	HB_PARAM_PARSING = -5,
	/** Number of bugs is zero. */
	HB_BUGS_ZERO = -6,
	/** Bugs exceed world slots. */
	HB_BUGS_OVERFLOW = -7,
	/** Bug's ideal temperature range overlaps. */
	HB_TEMPERATURE_OVERLAP = -8,
	/** Bug's max ideal temperature exceeds range. */
	HB_TEMPERATURE_OUT_RANGE = -9,
	/** Bug's output heat range overlap. */
	HB_OUTPUT_HEAT_OVERLAP = -10,
	/** Bug's max output heat exceeds range. */
	HB_OUTPUT_HEAT_OUT_RANGE = -11,
	/** Unable to open file. */
	HB_UNABLE_OPEN_FILE = -12,
	/** Memory alocation failed. */
	HB_MALLOC_FAILURE = -13
};



float average( const float *const vector, const size_t vsize )
{
	float sum = 0.0;

	for (size_t idx = 0; idx < vsize; idx++)
		sum += vector[ idx ];

	return sum / vsize;
}


#endif
