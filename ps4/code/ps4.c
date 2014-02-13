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

static void
distribute_vector(int n, double *v, int world_size)
{
	int rank, k, i;
	double *buf;

	k = n / world_size;
	buf = malloc(sizeof (double) * k);
	for(rank = 1; rank < world_size; ++rank) {
		for(i = 0; i < k; ++i) {
			buf[i] = v[i*world_size + rank];
		}
		//MPI_Send(buf, k, MPI_DOUBLE, rank, 100, MPI_COMM_WORLD);
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
		
	double pi_reference, sum_reference, sum;
	int rank, size;
	int n, k;
	double *v;
	/* timestamps */
	double t1, t2;
	
	/* mpi initialization */
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
	rank = 0;
	size = 1;
	#endif

	if (rank == 0) {
		
		#ifdef HAVE_OPENMP
		printf("OpenMP support enabled.\n");
		#endif
	
		#ifdef HAVE_MPI
		printf("Open MPI support enabled.\n");
		#endif
	}
	
	/* start timer */
	t1 = wall_time();

	MPI_Status status;

	if(rank == 0) {
		n = intpow(2, K_MAX);
		v = malloc(sizeof (double) * n);
		fill_invsquare_vector(n,v);
		distribute_vector(n,v,size-1);
	} else {
		MPI_Recv(v, 13, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, status);
	}



	/*double d = 2.0;
	double res = 0.1;
	
	MPI_Reduce(&d, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank == 0)
		printf("%f\n", res);
	*/
	
	/* stop timer */
	t2 = wall_time();
	
	if(rank == 0)
		printf("Total time: %e\n", t2-t1);

	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif

	exit(0);
	

	/* start timer *
	t1 = wall_time(); 

	/* reference values *
	pi_reference = 4.0l * atan(1.0);
	sum_reference = pi_reference * pi_reference / 6.0;
	
	/* allocating and filling vector with values 
	n = intpow(2, K_MAX);
	double *v = malloc(sizeof (double) * n);
	fill_invsquare_vector(n, v);
	

	/* summing and printing results *
	/*
	for(k = K_MIN; k <= K_MAX; ++k) {
		sum = sum_vector(intpow(2,k), v);
		printf("%i:\t %e\n", k, sum-sum_reference);
	}*/

	/* stop timer *
	t2 = wall_time();

	/* printing timing information *
	printf("Time: %e\n", t2-t1);

	/* preferred over return statements to end application *
	MPI_Finalize();
	exit(0);*/
}
