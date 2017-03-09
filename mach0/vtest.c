#include "mach0.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
	int ret = 0;
	double exp_val = M_PI;
	int k = 24;
	double error[24];

	FILE *f = fopen("vtest_result.txt", "w");

	for (int i = 1; i <= k; ++i) {
		double sol = mach0(pow(2, i));
		error[i-1] = fabs(exp_val - sol);
		fprintf(f, "n = 2^%d, error = %.16e\n", i, error[i-1]);
	}

	fclose(f);

	return 0;
}