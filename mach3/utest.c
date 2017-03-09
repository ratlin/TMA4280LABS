#include "mach3.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int ret = 0;
	int n = 3;

	double exp_val = 3.1416210293250346;
	double comp_val = mach3(n);

	if (rank == 0) {
		printf("\nResult expected value = %.16e; computed value = %.16e\n\n", exp_val, comp_val);

	double rel_err = fabs(1.0 - comp_val/exp_val);

	printf("Relative error = %.16e\n\n", rel_err);

	if (rel_err > 1.0e-15) {
		ret = 1;
	}

	printf("Unit test %s.\n", (ret == 0 ? "SUCCESS" : "FAILURE"));
	}
	
	MPI_Finalize();
	return 0;
}