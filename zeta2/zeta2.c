#include "zeta2.h"

#include <math.h>


double zeta2(int n) {
	double sum = 0.0;

	#pragma omp parallel for schedule(static) reduction(+:sum)
	for (int i = 1; i < n+1; i++) {
		double d = i;
		sum += 1.0/(d*d);
	}

	return sqrt(6*sum);
}
