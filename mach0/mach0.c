#include "mach0.h"

#include <math.h>

double arctan(double x, int n) {
	double sum = 0.0;

	for (int i = 1; i < n+1; i++) {
		sum += pow(-1, i-1)*pow(x, 2*i-1)/(2*i-1);
	}
	return sum;
}

double mach0(int n) {
	return (4*4*arctan(1.0/5.0, n)-4*arctan(1.0/239.0, n));
}