#ifndef P2MATH_H_
#define P2MATH_H_

void fill_up_vector_double(int n_length, double* a); 
void fill_up_matrix_double(int m, int n, double* A);
void print_vector(int n_length, double* a);
void print_matrix(int m, int n, double* A);
void double_vector_add_scalar_multiply(int n_length, double scalar, double* a, double* b, double* c);
void add_multiply_matrix(int n, double* a, double* A, double*b, double* c);
double dot_product(int n, double* a, double* b);

#endif
