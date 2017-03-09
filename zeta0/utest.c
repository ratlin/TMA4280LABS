#include "zeta0.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
	int ret = 0;
	int n = 3;

	double exp_val = 2.8577380332470415;
	double comp_val = zeta0(n);

	printf("\nResult expected value = %.16e; computed value = %.16e\n\n", exp_val, comp_val);

	double rel_err = fabs(1.0 - comp_val/exp_val);

	printf("Relative error = %.16e\n\n", rel_err);

	if (rel_err > 1.0e-15) {
		ret = 1;
	}

	printf("Unit test %s.\n\n", (ret == 0 ? "SUCCESS" : "FAILURE"));
	return 0;
}
