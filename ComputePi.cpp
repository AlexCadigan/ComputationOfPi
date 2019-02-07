#include <cmath>
#include <iostream>
#include "ComputePi.h"

int main() {
	double N = 100000000;

	sequential(N);
	// sequentialBroken(N);
	parallel(N);
}

void sequential(double N) {
	double x, y, Pi;
	for (double i = 1; i <= N; i ++) {
		x = (1 / N) * (i - 0.5);
		y = sqrt(1 - pow(x, 2));
		Pi += 4 * (y / N);
	}

	printf("Sequential:\t%f\n", Pi);
}

void sequentialBroken(double N) {
	double x, y, Pi;
	for (double i = 1; i <= N; i ++) {
		x = (1 / N) * (i - 0.5);
		y = sqrt(1 - pow(x, 2));
		Pi += 4 * (y / N);
	}

	printf("Sequential:\t%f\n", Pi);
}

void parallel(double N) {
	double x2, y2, Pi2;
	for (double i = 1; i <= N; i ++) {
		x2 = (1 / N) * (i - 0.5);
		y2 = sqrt(1 - pow(x2, 2));
		Pi2 += 4 * (y2 / N);
	}

	printf("Parallel:\t%f\n", Pi2);
}
