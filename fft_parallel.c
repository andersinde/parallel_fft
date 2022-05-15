#include <string.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include "functions.h"

#define PRINT 0
//#define D 20
//#define N 1048576


/*  compilation on dardel:
 *
 * cc -o program fft_parallel.c functions.c -lm -w -O
 *
 */

int main(int argc, char **argv) {

    int rank, size, tag, rc;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int D = atoi(argv[1]); // NOTE: input to program
    int N = pow(2,D);
    int n = N/size; // local size
    int k_global, exponent, q_rank, is_odd;
    int print_list = N<64;

    assert(isPowerOfTwo(N) && isPowerOfTwo(size));

    static complex *x;
    static complex *res;
    complex y[n];
    complex y_recv[n];
    double t0;

    if (rank == 0) { // only allocate memory for big arrays on master
	res = (complex *)malloc(N*sizeof(complex));
	x = (complex *)malloc(N*sizeof(complex));
    }
 
    srandom(123);

    for (int i=0; i<n; i++)
	y[i] = ((float)random())/(RAND_MAX);
    
    if (rank == 0)
	t0 = MPI_Wtime(); // start timer

    const complex w_N = cexp(-2i*PI/N);
    complex wk;


    static complex *twiddle_factors;

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
		//exponent = get_exponent(d, k_global, N); // todo: lookup table
		//int sign = exponent < N/2 ? 1 : -1;
		//wk = sign*twiddle_factors[exponent%(N/2)];
		exponent = get_exponent(d, k_global, N);
		wk = cpow(w_N, exponent);
		y[k] = wk*y[k] +    y_recv[k];
	    }
	} else {
	    MPI_Recv(&y_recv, 2*n, MPI_COMPLEX, q_rank, 0, MPI_COMM_WORLD, &status);
	    MPI_Send(&y,      2*n, MPI_COMPLEX, q_rank, 1, MPI_COMM_WORLD);

	    // add received list to current
	    for (int k = 0; k<n; k++) {
		k_global = rank*n + k;
		//exponent = get_exponent(d, k_global, N); // todo: lookup table
		//int sign = exponent < N/2 ? 1 : -1;
		//wk = sign*twiddle_factors[exponent%(N/2)];
		exponent = get_exponent(d, k_global, N);
		wk = cpow(w_N, exponent);
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
	double time =  MPI_Wtime() - t0; 
	char filename[50];
	sprintf(filename, "results/%dD.txt", D);

	FILE *fp; // append execution time to file
	fp = fopen(filename, "a");
	fprintf(fp, "%f %d %d\n", time, size, N);
	fclose(fp);

	printf("t = %f, P = %d, N = %d\n", time,size,N);
    }
    MPI_Finalize();
    return 0;
}

