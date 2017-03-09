#include "zeta2.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
	double s_time = omp_get_wtime();

	int ret = 0;
	double exp_val = M_PI;
	int k = 24;
	double error[24]; 

	FILE *f = fopen("vtest_result.txt", "w");

	for (int i = 1; i <= k; ++i) {
		s_time = omp_get_wtime();
		double sol = zeta2(pow(2, i));
		error[i-1] = fabs(exp_val - sol);
		fprintf(f, "n = 2^%d, error = %.16e", i, error[i-1]);
		double t_time = omp_get_wtime() - s_time;
		fprintf(f, " t = %.16e\n", t_time);	
	}

	fclose(f);

	return 0;
}
