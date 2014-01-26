#include "p2math.h"
#include <stdlib.h>
#include <stdio.h>

void double_vector_add_scalar_multiply(int n_length, double scalar, double* a, double* b, double* c)
{
	int i;
	for(i=0; i<n_length; ++i) {
		c[i] = a[i] + scalar*b[i];
	}
}

void fill_up_matrix_double(int m, int n, double* A)
{
	int i,j;
	for(i=0;i<m;++i) {
		for(j=0;j<n;++j) {
			A[n*i+j] = i*j;
		}
	}

}



void print_vector(int n_length, double* a)
{
	int i;
	for(i=0;i<n_length;++i) {
		printf("%4.4f ", a[i]);
	}
	printf("\n");
}

void print_matrix(int m, int n, double* A)
{
	int i,j;
	for(i=0;i<m;++i) {
		for(j=0;j<n;++j) {
			printf("%f ", A[n*i+j]);
		}
		printf("\n");
	}
}

void fill_up_vector_double(int n_length, double* a){
	int i;
	for(i=0;i<n_length;++i) {
		a[i] = i;
	}
}


void add_multiply_matrix(int n, double* a, double* A, double*b, double* c)
{
	int i,j;
	double d;
	for(i=0;i<n;++i){
		d = 0.0;
		for(j=0;j<n;++j){
			d += A[i*n + j]*b[j];
		}
		c[i] = a[i]+ d;
	}
}

double dot_product(int n, double* a, double* b)
{
	int i;
	double d = 0.0;
	for(i=0;i<n;++i) {
		d = d + a[i]*b[i];
	}
	return d;
}
