/*
 * Problem Set 4 Utility Library
 * TMA4280
 * Håkon Åmdal
 */

/* system headers */
#include <sys/time.h>
#include <stdlib.h>

/* local headers */
#include "util.h"

/*
 * Returns current wall time.
 * Function copied and modified from Arne Morten Kvarving's common library.
 */
double
wall_time()
{
	#ifdef HAVE_MPI
  	return MPI_Wtime();
	#else
  	struct timeval tmpTime;
  	gettimeofday(&tmpTime,NULL);
  	return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
	#endif
}

