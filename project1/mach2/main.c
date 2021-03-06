#include "mach2.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

	if (argc != 2) {
		printf("Input required, please set number of iterations, n.\n");
		return 1;
	}

	int n = atoi(argv[1]);
	double sol = mach2(n);
	printf("Solution using Machin formula with n = %d is %.16e\n\n", n, sol);
	return 0;
}
