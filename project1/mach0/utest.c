#include "mach0.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
	int ret = 0;
	int n = 3;

	double exp_val = 3.1416210293250346; // change exp_val
	double comp_val = mach0(n);

	printf("\nResult expected value = %.16e; computed value = %.16e\n\n", exp_val, comp_val);

	double rel_err = fabs(1.0 - comp_val/exp_val);

	printf("Relative error = %.16e\n\n", rel_err);

	if (rel_err > 1.0e-15) { // change limit for success
		ret = 1;
	}

	printf("Unit test %s.\n", (ret == 0 ? "SUCCESS" : "FAILURE"));
	return 0;
}