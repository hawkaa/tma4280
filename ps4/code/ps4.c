/*
 * Problem set 4
 * TMA4280
 * Håkon Åmdal
 */

/* system headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* local headers */
#include "util.h"

/*
 * Constants for determining the output, as specified in the task.
 */
static const int K_MIN = 3;
static const int K_MAX = 14;

/*
 * Timestamps
 */
static double t1, t2;

/*
 * MPI variables
 */ 
static int rank, size;
static MPI_Status status;
static int tag = 100;

/*
 * Fills a vector of length n with the inverse squared of the index (as 
 * specified in the task).
 */
static void
fill_partial_vector(const int n, double *vector, const int num_parts, const int offset)
{
	int i, j;
	#pragma omp parallel for schedule(static)
	for(i = 0; i < n; ++i) {
		j = i*num_parts + offset;
		vector[i] = 1.0 / (j * j + 2 * j + 1);
	}
}

/*
 * Calculates the sum of an array of doubles with the length n.
 */
static double
sum_vector(int n, const double *v)
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
intpow(const int b, int e)
{
	int i, r;
	r = 1;
	for(i = 0; i < e; ++i)
		r = r*b;
	return r;
}

static void
initialize(int argc, char** argv)
{
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
	rank = 0;
	size = 1;
	#endif
	
	if (rank == 0) {
		/* this could be relevant debug info */	
		#ifdef HAVE_OPENMP
		printf("OpenMP support enabled.\n");
		#endif
	
		#ifdef HAVE_MPI
		printf("Open MPI support enabled.\n");
		#endif
	}
	
	/* start timer */
	t1 = wall_time();

}

static void
finalize()
{
	
	/* stop timer */
	t2 = wall_time();
	
	if(rank == 0)
		printf("Total time: %e\n", t2-t1);

	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif

	exit(0);
}

/*
 * Program entry function.
 */
int
main(int argc, char** argv)
{
	/* variables for references and sum */
	double pi_reference, sum_reference, partial_sum, sum;

	/* vector to sum */
	double *v;
		
	/* various */
	int n, k, i, num_elements, partial_num_elements;

	initialize(argc, argv);

	/* reference values */
	pi_reference = 4.0l * atan(1.0);
	sum_reference = pi_reference * pi_reference / 6.0;
	
	/* allocate room for vector */
	n = intpow(2, K_MAX) / size;
	v  = malloc(sizeof (double) * n);

	/*
	 * Giving everyone a partial vector
	 */ 
	if(rank == 0) {
		int i;
		for(i = 1; i < size; ++i){
			fill_partial_vector(n, v, size, i);
			MPI_Send(v, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD); 
		}
		fill_partial_vector(n, v, size, rank);
	} else {
		MPI_Recv(v, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);	
	}

	/* now everyone has their partial vector stored in v */

	/* iterating over k's */
	for(k = K_MIN; k <= K_MAX; ++k) {
		
		/* how many elements are needed */
		num_elements = intpow(2, k);
		partial_num_elements = num_elements / size;
		
		/* summing vector */
		partial_sum = sum_vector(partial_num_elements, v);

		/* reducing into p0 */
		MPI_Reduce (&partial_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if(rank == 0) {
			/* reporting if p0 */
			printf("%i:\t %e\n", k, sum-sum_reference);
		}
	}
	
	/* finalizing */
	finalize();
}
