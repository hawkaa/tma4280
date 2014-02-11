/*
 * Problem set 4
 * TMA4280
 * Håkon Åmdal
 */

/* system headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* local headers */
#include "util.h"

/*
 * Constants for determining the output, as specified in the task.
 */
static const int K_MIN = 3;
static const int K_MAX = 14;


/*
 * Fills a vector of length n with the inverse squared of the index (as 
 * specified in the task).
 */
static void
fill_invsquare_vector(int n, double *v)
{
	int i;
	#pragma omp parallel for schedule(static)
	for(i = 0; i < n; ++i) {
		v[i] = 1.0 / (i * i + 2 * i + 1);
	}
}

/*
 * Calculates the sum of an array of doubles with the length n.
 */
static double
sum_vector(int n, double *v)
{
	double sum;
	int i;
	sum = 0.0;
	#pragma omp parallel for schedule(static) reduction(+:sum)
	for(i = 0; i < n; ++i) 
		sum += v[i];	
	return sum;
}

/*
 * Power function for integers.
 */
static int
intpow(int b, int e)
{
	int r;
	r = b;
	while(e > 1){
		r = r * b;
		--e;
	}
	return r;
}

/*
 * Program entry function.
 */
int
main(int argc, char** argv)
{
	#ifdef HAVE_OPENMP
	printf("OpenMP suppoort enabled\n");
	#endif
	double pi_reference, sum_reference, sum;
	int n, k;

	/* timestamps */
	double t1, t2;

	/* start timer */
	t1 = wall_time(); 

	/* reference values */
	pi_reference = 4.0 * atan(1.0);
	sum_reference = pi_reference * pi_reference / 6.0;
	
	/* allocating and filling vector with values */
	n = intpow(2, K_MAX);
	double *v = malloc(sizeof (double) * n);
	fill_invsquare_vector(n, v);

	/* summing and printing results */
	for(k = K_MIN; k <= K_MAX; ++k) {
		sum = sum_vector(intpow(2,k), v);
		printf("%i:\t %e\n", k, sum-sum_reference);
	}

	/* stop timer */
	t2 = wall_time();

	/* printing timing information */
	printf("Time: %e\n", t2-t1);

	/* preferred over return statements to end application */
	exit(0);
}
