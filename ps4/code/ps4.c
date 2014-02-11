
#include <stdio.h>
#include <stdlib.h>

#include "ps4_util.h"

static int n = 5;

int main(int argc, char** argv)
{
	
	double *v = (double*)malloc(sizeof(double)*n);
	
	fill_vector(n, v);
	printf("hehe%f\n", v[2]);
	return 0;
}
