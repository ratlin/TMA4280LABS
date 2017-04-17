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
	MPI_Init(&argc, &argv);
	double s_time = MPI_Wtime();
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
		s_time = MPI_Wtime();
		double sol = zeta1(pow(2, i));
		error[i-1] = fabs(exp_val - sol);
		if (rank == 0) {
			fprintf(f, "n = 2^%d, error = %.16e ", i, error[i-1]);
			double t_time = MPI_Wtime() - s_time;
			fprintf(f, " t = %.16e\n", t_time);
		}
	}

	if (rank == 0) {
		fclose(f);
	}
	MPI_Finalize();

	return 0;
}
