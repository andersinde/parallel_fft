#include <string.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include "functions.h"

#define PRINT 0

int main(int argc, char **argv) {

    int rank, size, tag, rc;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
	MPI_Wtime(); // start timer

    int D = atoi(argv[1]); // NOTE: input to program
    int N = pow(2,D);
    int n = N/size; // local size
    int k_global, exponent, q_rank, is_odd;
    int print_list = N<64;

    assert(isPowerOfTwo(N) && isPowerOfTwo(size));

    complex x[N]; // store big array on rank 0
    complex *res;
    complex y[n];
    complex y_recv[n];

    if (rank == 0) { // only allocate memory for big arrays on master
	res = (complex *)malloc(N*sizeof(complex));
	//x = (complex *)malloc(N*sizeof(complex));
    }
 
    srandom(123);

    ///*
    for (int i=0; i<N; i++)
	x[i] = ((float)random())/(RAND_MAX);


    //NOTE: temporary code below
    for (int i=0; i<n; i++)
	y[i] = x[rank*n + i]; // initialize result list y
    //*/

    const complex w_N = cexp(-2i*PI/N);
    complex wk;

    complex twiddle_factors[N/2]; // is it necessary to init whole array on every processor?

    for (int i=0; i<N/2; i++)
	twiddle_factors[i] = cpow(w_N, i);

    /*
	ALT1
    exponent = get_exponent(d, k_global, N);
    int sign = exponent < N/2 ? 1 : -1;
    wk = sign*twiddle_factors[exponent%(N/2)];

    ALT2
    exponent = get_exponent(d, k_global, N);
    wk = cpow(w_N, exponent);
    int exponent = get_exponent(d, k, N);
    wk = cpow(w_N, exponent);
     */

    //char filename[8+2];
    //sprintf(filename, "log/%d.txt", rank);
    //FILE *file  = fopen(filename, "w");
    //if (file == NULL) printf("Error. Could not open file\n"); 


    /*
     * FFT algorithm, log(N) steps
     *
     * k        = local element in processor
     * k_global = global element
     * q 	= global element to communicate with
     * q_rank	= other processor to communicate with
     * q_local  = local element in the processor to communicate with
     *
     * N	= global array size
     * n	= number of elements per processor: N/size
     */

    // MAIN LOOP
    for (int d=D-1; d>=0; d--) {
	// use bitflip to aquire butterfly "neighbour" q, i.e. the processor to exhange with
	q_rank = bitflip(rank, d - D + log2(size)); // rank of processor to communicate with
	is_odd = (rank*n >> d) & 1; // determines how to add up the numbers


	// SEQUENTIAL STEP
	if (rank == q_rank) { 
	    // communication steps are over: perform fft on processor (sequential)
	    seq_fft_proc(y, rank, size, n, N, twiddle_factors);
	    break;
	}

	// COMMUNICATION STEP
	if (is_odd) {
	    MPI_Send(&y,      2*n, MPI_COMPLEX, q_rank, 0, MPI_COMM_WORLD);
	    MPI_Recv(&y_recv, 2*n, MPI_COMPLEX, q_rank, 1, MPI_COMM_WORLD, &status);

	    // add received list to current
	    for (int k = 0; k<n; k++) {
		k_global = rank*n + k;
		exponent = get_exponent(d, k_global, N); // todo: lookup table
		int sign = exponent < N/2 ? 1 : -1;
		//wk = twiddle_factors[exponent];
		wk = sign*twiddle_factors[exponent%(N/2)];
		y[k] = wk*y[k] +    y_recv[k];
	    }
	} else {
	    MPI_Recv(&y_recv, 2*n, MPI_COMPLEX, q_rank, 0, MPI_COMM_WORLD, &status);
	    MPI_Send(&y,      2*n, MPI_COMPLEX, q_rank, 1, MPI_COMM_WORLD);

	    // add received list to current
	    for (int k = 0; k<n; k++) {
		k_global = rank*n + k;
		exponent = get_exponent(d, k_global, N); // todo: lookup table
		int sign = exponent < N/2 ? 1 : -1;
		//wk = twiddle_factors[exponent];
		wk = sign*twiddle_factors[exponent%(N/2)];
		y[k] =    y[k] + wk*y_recv[k];
	    }
	}
    }

    // at this point all processor contain their part of the final result in y
    // send y to proc one and store it in res

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(y, 2*n, MPI_COMPLEX, res, 2*n, MPI_COMPLEX, 0, MPI_COMM_WORLD);

    // algorithm finished
    if (rank == 0) {
	double t =  MPI_Wtime(); 
	char filename[50];
	sprintf(filename, "results/%d.txt", N);

	// verification

	bit_reverse_list(res,N);

	// perform seq fft on original input to compare results
	seq_fft(x, x, N);
	//seq_fft_proc(x, 0, 1, n, N, twiddle_factors);
	//bit_reverse_list(x,N);

	float SS_res = 0;
	for (int i=0; i<N; i++)
	    SS_res += pow(crealf(x[i]-res[i]),2) + pow(cimagf(x[i]-res[i]),2);

	FILE *fp; // append execution time to file
	fp = fopen(filename, "a");
	fprintf(fp, "%d %d %f %f\n", N, size, t, SS_res);
	fclose(fp);

	printf("N = %d, P = %d, Time = %f, Error = %f\n", N,size,t,SS_res);
    }
    MPI_Finalize();

    return 0;
}

