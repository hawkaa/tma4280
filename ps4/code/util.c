/*
 * Problem Set 4 Utility Library
 * TMA4280
 * Håkon Åmdal
 */

/* system headers */
#include <sys/time.h>
#include <stdlib.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* local headers */
#include "util.h"

/*
 * Returns current wall time.
 * Function copied a from Arne Morten Kvarving's common library.
 */
double
wall_time()
{
	#ifdef HAVE_MPI
  	return MPI_Wtime();
	#elif defined(HAVE_OPENMP)
  	return omp_get_wtime();
	#else
	struct timeval tmpTime;
	gettimeofday(&tmpTime,NULL);
	return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
	#endif
}

/*
 * Power function for integers.
 * Created because the pow in stdlib did not accept integer constants.
 */
int
intpow(const int b, int e)
{
	int i, r;
	r = 1;
	for(i = 0; i < e; ++i)
		r = r*b;
	return r;
}


