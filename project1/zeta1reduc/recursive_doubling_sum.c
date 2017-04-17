#include "recursive_doubling_sum.h"
#include <mpi.h>

double recursive_doubling_sum(double value) {
	int rank, nprocs;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	for (int i = 1; i < nprocs;) {
		double received_value;
		if ((rank/i)%2 == 0 && rank+i < nprocs) {
			MPI_Send(&value, 1, MPI_DOUBLE, rank+i, 100, MPI_COMM_WORLD);
			MPI_Recv(&received_value, 1, MPI_DOUBLE, rank+i, 100, MPI_COMM_WORLD, &status);
		}
		else {
			MPI_Recv(&received_value, 1, MPI_DOUBLE, rank-i, 100, MPI_COMM_WORLD, &status);
			MPI_Send(&value, 1, MPI_DOUBLE, rank-i, 100, MPI_COMM_WORLD);
		}
		i = i<<1; // i is bit-shifted one step up (0001 to 0010)
		value += received_value;
	} 
}