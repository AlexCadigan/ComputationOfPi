/*
Jennifer Cho + Alex Cadigan
2/6/2019
COMP-481 Final Project
*/

#include <cmath>
#include <iostream>
#include <omp.h>
#include "ComputePi.h"

int main() {
	double N = 1000;

	double start = omp_get_wtime();
	sequential(N);
	double timeS = omp_get_wtime() - start;
	start = omp_get_wtime();
	parallel(N);
	double timeP = omp_get_wtime() - start;
	printf("Sequential Time:\t%f\n", timeS);
	printf("Parallel Time:\t%f\n", timeP);
}

void sequential(double N) {
	double x, y, Pi = 0;
	for (double i = 1; i <= N; i ++) {
		x = (1 / N) * (i - 0.5);
		y = sqrt(1 - pow(x, 2));
		Pi += 4 * (y / N);
	}

	printf("Sequential:\t%f\n", Pi);
}

void parallel(double N) {
	double Pi = 0;
	#pragma omp parallel 
	{
		double x, y;
		#pragma omp for reduction(+:Pi)
		for (int i = 1; i <= (int) N; i ++) {
			x = (1 / N) * (i - 0.5);
			y = sqrt(1 - pow(x, 2));
			Pi += 4 * (y / N);
		}
	}

	printf("Parallel:\t%f\n", Pi);
}
