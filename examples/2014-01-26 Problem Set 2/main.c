#include <stdlib.h>
#include <stdio.h>
#include "p2math.h"

int main(int argc, char** argv) {
	if (argc < 3) {
		printf("Not enough arguments\n");
		return 1;
	}
	
	int n;
	double scalar;
	
	n = atoi(argv[1]);
	scalar = strtod(argv[2], NULL);

	double *A, *a, *b, *x, *y;
	double alpha;
	A = malloc(sizeof(double)*n*n);
	a = malloc(sizeof(double)*n);
	b = malloc(sizeof(double)*n);
	x = malloc(sizeof(double)*n);
	y = malloc(sizeof(double)*n);

	
	fill_up_vector_double(n, a);
	fill_up_vector_double(n, b);
	fill_up_matrix_double(n, n, A);
	double_vector_add_scalar_multiply(n, scalar, a, b, x);
	add_multiply_matrix(n, a, A, b, y);

	alpha = dot_product(n, x,y);
//	print_vector(n, x);
//	print_vector(n, y);
	printf("%f\n", alpha);

}
