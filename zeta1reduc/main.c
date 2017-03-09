#include "zeta1.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char** argv) {
	clock_t start_t, diff_t;
	start_t = clock();

	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc != 2) {
		if (rank == 0) {
			printf("Input required, please set number of iterations, n.\n");
		}
		MPI_Finalize();
		return 1;
	}
	int n = atoi(argv[1]);
	
	double sol = zeta1(n);
	if (rank == 0) {
		printf("\nSolution using Riemann zeta function with n = %d is %.16e\n", n, sol);
		diff_t = (clock() - start_t);
		double total_t = (double) diff_t/CLOCKS_PER_SEC;
		printf("Walltime: %f\n\n", total_t);
	}
	
	MPI_Finalize();
	return 0;
}
