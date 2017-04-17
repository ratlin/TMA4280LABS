/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 * Parallelized by Pauline H. Jensen, April 2017
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
real rhs(real x, real y);
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int nprocs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int numthreads = omp_get_max_threads();

	// NOT SURE, MAYBE ONLY RANK 0 SHOULD PRINT THIS?
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
    int nn = 4 * n;
    real h = 1.0 / n;

    // SPLIT MATRIX INTO COLUMNS BASED ON nprocs
    // lokal antall kolonner og lokal displacement
    // ta hensyn til kolonner til overs

    // Grid points
    real *grid = mk_1D_array(n+1, false);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    // The diagonal of the eigenvalue matrix of T
    real *diag = mk_1D_array(m, false);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    // Initialize the right hand side data
    // HERE JUST THE PART THE PROCESS OWNS SHOULD BE INITIALIZED
    real **b = mk_2D_array(m, m, false);
    real **bt = mk_2D_array(m, m, false);
    // INITIALIZING THE z FOR THE SINE TRANSFORM
    real *z = mk_1D_array(nn, false);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            b[i][j] = h * h * rhs(grid[i], grid[j]);
        }
    }

    // Calculate Btilde^T = S^-1 * (S * B)^T
    // CHANGE SOME THINGS HERE TO ACCOMODATE THE PARTS
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        fst_(b[i], &n, z, &nn);
    }
    transpose(bt, b, m);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        fstinv_(bt[i], &n, z, &nn);
    }

    // Solve Lambda * Xtilde = Btilde
    // CHANGE SOME RANGES AND STUFF
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[i] + diag[j]);
        }
    }

    // Calculate X = S^-1 * (S * Xtilde^T)
    // CHANGE SOME STUFF
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z, &nn);
    }
    transpose(b, bt, m);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z, &nn);
    }

    // Calculate maximal value of solution
    // CHANGE RANGES AND PUT IN SOME MORE INFO, LIKE ERROR MAX
    double u_max = 0.0;
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            u_max = u_max > b[i][j] ? u_max : b[i][j];
        }
    }

    // SOME MPI_Reduce TO FIND ACTUAL MAXIMUM AND ERROR MAX
    // CODE HERE


    // FIX SO ONLY RANK 0 SHOULD PRINT MAX AND ERROR MAX AND TIMING
    printf("u_max = %e\n", u_max);

    MPI_Finalize();
    return 0;
}

real rhs(real x, real y) {
	// CHANGE TO RETURN 5*PI*PI*sin(PI*x)*sin(2.0*PI*y)
    return 2 * (y - y*y + x - x*x);
}

// POSSIBLY MAKE ANOTHER RHS WITH RETURN VALUE: sin(PI*x)*sin(2.0*PI*y)

// CHANGE SO IT HAS CORRECT INPUT
void transpose(real **bt, real **b, size_t m)
{
	// CREATE SOME VECTORS TO SEND AND RECEIVE INFO

	// CREATE SOME INTS TO KEEP TRACK OF SEND/RECEIVE COUNT AND DISPLACEMENT
	// MAKE A FOR LOOP TO SET THESE INTS CORRECTLY, DO NOT USE PRAGMA ON THE LOOP
	// HMMMM....

	// FIXING THIS IS THE HARD PART!!!!!!!! !!! !
	// FIRST MAKE DOUBLE FOR TO PUT CONTENTS INTO SEND VECTOR
	// THEN USE MPI_Alltoallv() TO COMMUNICATE SEND VECTORS
	// THEN A TRIPLE(?) FOR TO PUT THE RECEIVED DATA INTO THE MATRIX
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
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
