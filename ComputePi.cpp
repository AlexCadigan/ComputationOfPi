/*
Jennifer Cho + Alex Cadigan
2/6/2019
COMP-481 Final Project
*/

#include <cmath>
#include <iostream>
#include "ComputePi.h"

int main() {
	double N = 10;

	sequential(N);
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

void parallel(double N) {
	double x, y, Pi = 0;
	for (double i = 1; i <= N; i ++) {
		x = (1 / N) * (i - 0.5);
		y = sqrt(1 - pow(x, 2));
		Pi += 4 * (y / N);
	}

	printf("SequentialBroken:\t%f\n", Pi);
}
