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
#ifdef HAVE_MPI
static MPI_Status status;
static int tag = 100; 		/* arbitary value */
#endif

/*
 * Fills a partial vector of length n with the inverse squared of the index
 * (as specified in the problem set). The function will only fill parts of the
 * vector, dependent on the number of parts the main vector is split into
 * (world size) and offset (processor rank). 
 */
static void
fill_partial_vector(const int n, double *vector, const int num_parts, const int offset)
{
	int i, j;
	#pragma omp parallel for schedule(static) private(j)
	for(i = 0; i < n; ++i) {
		j = (i * num_parts + offset) + 1;
		vector[i] = 1.0 / (j * j);
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
 * Initialize the program by setting MPI variables and starting timers.
 */
static void
initialize(int argc, char** argv)
{
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
	/* if no mpi, current processor is 0 and size of world is 1 */
	rank = 0;
	size = 1;
	#endif
	
	if (rank == 0) {
		/* printing relevant info of modules enabled */	
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

/*
 * Ending the program by printing timing results and finalizing MPI.
 */
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
 * Program main function.
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
	
	/* initializing */
	initialize(argc, argv);

	/* reference values */
	pi_reference = 4.0l * atan(1.0);
	sum_reference = pi_reference * pi_reference / 6.0;
	
	/* allocate room for (partial) vector */
	n = intpow(2, K_MAX) / size;
	v  = malloc(sizeof (double) * n);

	/*
	 * Giving everyone a partial vector, a share of the main vector to
	 * sum.
	 */ 
	if(rank == 0) {

		int i;

		/* one for every other processor, if mpi enabled */
		#ifdef HAVE_MPI
		for(i = 1; i < size; ++i){
			fill_partial_vector(n, v, size, i);
			MPI_Send(v, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD); 
		}
		#endif

		/* one for itself */
		fill_partial_vector(n, v, size, rank);

	} else {
		#ifdef HAVE_MPI
		MPI_Recv(v, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		#endif
	}

	/* now everyone has their part of the main vector stored in v */

	/* iterating over k's, as specified in the problem set */
	for(k = K_MIN; k <= K_MAX; ++k) {
		
		/* how many elements are needed? */
		num_elements = intpow(2, k);
		partial_num_elements = num_elements / size;
		
		/* summing the vector */
		partial_sum = sum_vector(partial_num_elements, v);

		/* reducing into p0 */
		#ifdef HAVE_MPI
		MPI_Reduce (&partial_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		#else
		sum = partial_sum;
		#endif

		if(rank == 0) {
			/* reporting if p0 */
			printf("%i:\t %e\n", k, sum_reference - sum);
		}
	}

	/* freeing up the vector */
	free(v);
	
	/* finalizing */
	finalize();
}
