#include "ps4_util.h"

void fill_vector(int n, double* v)
{
	/*
	 * This function fills a vector with elements
	 */
	int i;
	for(i=0; i<n; ++i)
		v[i] = 1/i*i;
}
