#include<complex.h>
#include<math.h>
#include<stdio.h>

#define PI 3.14159265

void cprintlist(complex *c, int);
void printlist(float *x, int);
//void cprint_abs(complex c);
void cprint(complex c);
void cprintlist_tofile(FILE *file, complex *x, int N);
int isPowerOfTwo (unsigned int x);
unsigned int bit_reversed(unsigned int, int);
int bitflip(int p, int d);
int get_exponent(int d, int k_global, int N);

void seq_fft(complex *x, complex *y, int N);
void seq_fft_proc(complex *y, int rank, int size, int n, int N, complex *twiddle_factors);
//void seq_fft_proc(complex *y, int rank, int size, int n, int N);
void bit_reverse_list(complex *x, int N);
