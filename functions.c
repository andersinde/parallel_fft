#include "functions.h"
#include<stdio.h>
#include <assert.h>

void cprint(complex c) {
    printf("%f%+fi\n", crealf(c), cimagf(c));
}
void cprintlist(complex *x, int N) {
    printf("{\n");
    for (int i=0; i<N; i++)
	cprint(x[i]);
    printf("}\n");
    printf("\n");
}
void cprintlist_tofile(FILE *file, complex *x, int N) {
    fprintf(file,"{\n");
    for (int i=0; i<N; i++)
	fprintf(file, "%f%+fi\n", crealf(x[i]), cimagf(x[i]));
	//cprint(x[i]);
    fprintf(file,"}\n");
}
void printlist(float *x, int N) {
    printf("{\n");
    for (int i=0; i<N; i++) \
	printf("%f\n", x[i]);
    printf("}\n");
    printf("\n");
}
int isPowerOfTwo (unsigned int x) {
  return ((x != 0) && !(x & (x - 1)));
}

// Function to reverse bits of num
unsigned int bit_reversed(unsigned int num, int NO_OF_BITS) {
    unsigned int reverse_num = 0;
    for (int i = 0; i < NO_OF_BITS; i++)
        if ((num & (1 << i)))
            reverse_num |= 1 << ((NO_OF_BITS - 1) - i);

    return reverse_num;
}

//  returns p with the the d'th bit flipped
int bitflip(int p, int d) {
    int mask = pow(2, d);
    return p ^ mask;
}

// returns exponent to w_N as function of d, k and N
int get_exponent(int d, int k, int N) {
    return ((int)pow(2,d)) * bit_reversed(k,(int) log2(N)) %N;
    // = 2**d * bit_reversed(k, log(N)) mod N
}

void bit_reverse_list(complex *x, int N) {
    complex temp;
    int D = (int)log2(N);
    for (int i=1; i<(int)N/2; i++) {
	temp = x[i];
	x[i] = x[bit_reversed(i,D)];
	x[bit_reversed(i,D)] = temp;
    }
}

void seq_fft(complex *x, complex *y, int N) {
    int D = (int) log2(N);
    complex w_N = cexp(-2i*PI/N);
    complex wk;

    // initialize result list y
    for (int i=0; i<N; i++)
       	y[i] = x[i];

    // aux list to have constant y per FFT step
    complex y_aux[N];
    for (int i=0; i<N; i++)
	y_aux[i] = x[i];


    // every step of FFT algorithm, log(N) steps
    for (int d=D-1; d>=0; d--) {
	for (int k = 0; k<N; k++) {

	    // use bitflip to aquire butterfly "neighbour"
	    int q = bitflip(k,d);

	    // calculate wiggle factor
	    //int exponent = ((int)pow(2,d)) *bit_reversed(k,(int) log2(N)) %N;
	    int exponent = get_exponent(d, k, N);
	    wk = cpow(w_N, exponent);

	    // figure out which term should be multiplied by w
	    //(num >> n) & 1); returns the n'th bit of num
	    int is_odd = (k >> d) & 1;

	    if (is_odd)
		y[k] = wk * y[k] + y_aux[q]; // y_aux is constant throughout k loop
	    else
		y[k] = y[k] + wk * y_aux[q];
	}

	// update y_aux for next round
	for (int i=0; i<N; i++)
	    y_aux[i] = y[i];
    }
    bit_reverse_list(y,N);
}

// perform seq fft on section of a bigger array starting at index rank
void seq_fft_proc(complex *y, int rank, int size, int n, int N, complex *twiddle_factors) {

    int D = (int) log2(n);
    complex w_N = cexp(-2i*PI/N);
    complex wk;
    int q, k_global, exponent, is_odd;

    // aux list to have constant y per FFT step
    complex y_aux[n];

    // every step of FFT algorithm, log(n) steps
    for (int d=D-1; d>=0; d--) {

	// auxiliary, to y, list that remains constant through k loop
	for (int i=0; i<n; i++)
	    y_aux[i] = y[i];

	for (int k = 0; k<n; k++) {

	    // use bitflip to aquire butterfly "neighbour"
	    q = bitflip(k,d);
	    assert(q>=0);
	    assert(q<n);

	    // get twiddle factor, NOTE, need to use global k since this
	    // is a part of a larger list
	    k_global = rank*n + k;
	    exponent = get_exponent(d, k_global, N);
	    int sign = exponent < N/2 ? 1 : -1;
	    //wk = twiddle_factors[exponent];
	    wk = sign*twiddle_factors[exponent%(N/2)];

	    // figure out which term should be multiplied by w
	    //(num >> n) & 1); returns the n'th bit of num
	    is_odd = (k >> d) & 1;

	    if (is_odd)
		y[k] = wk * y[k] + y_aux[q];
	    else
		y[k] = y[k] + wk * y_aux[q];
	}
    }
}


