#include "mach1.h"

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
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

	double sol = mach1(n);
	if (rank == 0) {
		printf("Solution using Machin formula with n = %d is %.16e\n\n", n, sol);
	}

	MPI_Finalize();
	return 0;
}
