/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 * Parallelized by Pauline Jensen, April 2017
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
// INPUT TO TRANSPOSE MUST CHANGE
void transpose(real **bt, real **b, size_t m);
void mpi_transpose(real **bt, real **b, int *send_displ, int *displ, int *send_count, int m, int num_col, real *send_buf, real *rec_buf, int rank, int nprocs);
real rhs(real x, real y);
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int nprocs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int num_threads = omp_get_max_threads();
	double global_u_max = 0;

    if (argc < 2) {
        printf("Usage:\n");
        printf("  poisson n\n\n");
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

    double start_time;
    if (rank == 0)
    {
    	start_time = MPI_Wtime();
    }

    // The number of grid points in each direction is n+1
    // The number of degrees of freedom in each direction is n-1
    int n = atoi(argv[1]);
    int m = n - 1;
    int leftover_col = m%nprocs;
    int nn = 4 * n;
    real h = 1.0 / n;

    int *count = (int *) malloc(nprocs*sizeof(int));
    int *displ = (int *) malloc((nprocs+1)*sizeof(int));
    displ[nprocs] = m;
    displ[0] = 0;

    for (int i = 0; i < nprocs; i++)
    {
    	count[i] = m/nprocs;
    	if (leftover_col != 0)
    	{
    		count[i]++;
    		leftover_col--;
    	}
    	if (i < nprocs-1)
    	{
    		displ[i+1] = displ[i]+count[i];
    	}
    }

    int num_col = count[rank];


    // Grid points
    real *grid = mk_1D_array(n+1, false);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    // The diagonal of the eigenvalue matrix of T
    real *diag = mk_1D_array(m, false);
    real *send_buf = mk_1D_array(num_col*m, false);
    real *rec_buf = mk_1D_array(num_col*m, false);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    // Initialize the right hand side data
    real **b = mk_2D_array(num_col, m, false);
    real **bt = mk_2D_array(num_col, m, false);
    real **z = mk_2D_array(num_threads, nn, false);
    int *send_count = (int *) malloc((nprocs+1)*sizeof(int));
    int *send_displ = (int *) malloc((nprocs+1)*sizeof(int));
    send_displ[0] = 0;
    for (int i = 0; i < nprocs; i++)
    {
    	send_count[i] = count[i]*count[rank];
    	send_displ[i] = displ[i]*count[rank];
    }
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_col; i++) {
        for (size_t j = 0; j < m; j++) {
            b[i][j] = h * h * rhs(grid[i+displ[rank]], grid[j]);
        }
    }

    // Calculate Btilde^T = S^-1 * (S * B)^T
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_col; i++) {
        fst_(b[i], &n, z[omp_get_thread_num()], &nn);
    }
    //transpose(bt, b, m);
    mpi_transpose(bt, b, send_displ, displ, send_count, m, num_col, send_buf, rec_buf, rank, nprocs);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_col; i++) {
        fstinv_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }

    // Solve Lambda * Xtilde = Btilde
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_col; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[i+displ[rank]] + diag[j]);
        }
    }

    // Calculate X = S^-1 * (S * Xtilde^T)
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_col; i++) {
        fst_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }
    //transpose(b, bt, m);
    mpi_transpose(b, bt, send_displ, displ, send_count, m, num_col, send_buf, rec_buf, rank, nprocs);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_col; i++) {
        fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
    }

    // Calculate maximal value of solution
    double u_max = 0.0;
    double current_max;
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_col; i++) {
        for (size_t j = 0; j < m; j++) {
        	current_max = b[i][j] - sin(PI*(grid[displ[rank]+1]))*sin(2.0*PI*grid[j]);
        	if (current_max > u_max)
        	{
        		u_max = current_max;
        	}
        }
    }

    MPI_Reduce(&u_max, &global_u_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
    	printf("u_max = %e\n", global_u_max);
    }

    MPI_Finalize();
    return 0;
}

real rhs(real x, real y) {
	// CHANGE BETWEEN: 	return 5*PI*PI*sin(PI*x)*sin(2.0*PI*y);
    //					return 2 * (y - y*y + x - x*x);
    return 5*PI*PI*sin(PI*x)*sin(2.0*PI*y);
}


// POSSIBLY MAKE ANOTHER RHS WITH RETURN VALUE: sin(PI*x)*sin(2.0*PI*y)

void transpose(real **bt, real **b, size_t m)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
}

void mpi_transpose(real **bt, real **b, int *send_displ, int *displ, int *send_count, int m, int num_col, real *send_buf, real *rec_buf, int rank, int nprocs)
{
	int counter = 0;

	for (int i = 0; i < nprocs; i++)
	{
		for (int j = 0; j < num_col; j++)
		{
			for (int k = displ[i]; k < displ[i+1]; k++)
			{
				send_buf[counter] = b[j][k];
				counter++;
			}
		}
	}

	MPI_Alltoallv(send_buf, send_count, send_displ, MPI_DOUBLE, rec_buf, send_count, send_displ, MPI_DOUBLE, MPI_COMM_WORLD);

	counter = 0;

	for (int i = 0; i < nprocs; i++)
	{
		for (int j = displ[i]; j < displ[i+1]; j++)
		{
			for (int k = 0; k < num_col; k++)
			{
				bt[k][j] = rec_buf[counter];
				counter++;
			}
		}
	}
}

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    real **ret = (real **)malloc(n1 * sizeof(real *));

    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }

    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}
