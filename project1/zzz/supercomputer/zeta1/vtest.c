#include "zeta1.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char** argv) {
	clock_t start_t, diff_t;
	start_t = MPI_Wtime();

	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	FILE *f;

	int ret = 0;
	double exp_val = M_PI;
	int k = 24;
	double error[24]; 

	if (rank == 0) {
		f = fopen("vtest_result.txt", "w");
	}

	for (int i = 1; i <= k; ++i) {
		double sol;
		for (int b=0; b < 1000; ++b)
		{
			sol = zeta1(pow(2, i));
		}
		error[i-1] = fabs(exp_val - sol);
		if (rank == 0) {
			fprintf(f, "n = 2^%d, error = %.16e ", i, error[i-1]);
			diff_t = (MPI_Wtime() - start_t);
			double total_t = (double) diff_t/CLOCKS_PER_SEC;
			fprintf(f, " t = %.16e\n", total_t);
		}
	}

	if (rank == 0) {
		fclose(f);
	}
	MPI_Finalize();

	return 0;
}
