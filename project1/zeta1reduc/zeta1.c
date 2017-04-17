#include "zeta1.h"
#include "recursive_doubling_sum.h"

#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double zeta1(int n) {
	int nprocs, rank, tag;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	tag = 100;
	int m = n/nprocs;
	int r = n%nprocs;
	int my_size = m + (rank<r);
	double* vec;

	if (rank == 0) {
		vec = malloc(n*sizeof(double));
		
		for (int i = 1; i <= n; i++) {
			double d = i;
			vec[i-1] = 1.0/(d*d);
		}

		double * vi = vec + my_size;

		for (int i = 1; i < nprocs; i++) {
			if(i < r){
				MPI_Send(vi, m+1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
				vi += m+1;
			}
			else {
				MPI_Send(vi, m, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
				vi += m;
			}
		}	
	}
	else {
		vec = malloc(my_size*sizeof(double));
		MPI_Recv(vec, my_size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	}

	double partial_sum = 0.0;
	for (int i = my_size - 1; i >= 0; i--) {
		partial_sum += vec[i];
	}

	free(vec);

	double sum = recursive_doubling_sum(partial_sum);

	if (rank == 0) {
		return sqrt(6*sum);
	}
	else {
		return 0;
	}
}