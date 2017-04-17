#include "recursive_doubling_sum.h"

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double sum_of_ranks = recursive_doubling_sum(rank);
	if (rank == 0) {
		printf("blabla %e\n", sum_of_ranks);
	}
	MPI_Finalize();
	return 0;
}