/*
Jennifer Cho + Alex Cadigan
2/6/2019
COMP-481 Final Project
*/

#include <cmath>
#include <iostream>
#include <omp.h>
#include "ComputePi.h"

using namespace std;

int main() {
	// Gets user input:
	int numTrials, subInt;
	cout << "Enter the number of trials to run:\t";
	cin >> numTrials;
	cout << "Enter the number of sub-intervals on interval [0, 1] to use:\t";
	cin >> subInt;

	// Stores the average times
	double start, seqAve = 0, intAve;

	// Runs sequential simulations
	for (int trialNum = 0; trialNum < numTrials; trialNum ++) {
		start = omp_get_wtime();
		seq(subInt);
		seqAve += omp_get_wtime() - start;
	}
	printf("Sequential Time:\t%f\n", seqAve / numTrials);

	// Runs through the different thread amounts
	for (int numThreads = 2; numThreads <= 64; numThreads *= 2) {
		intAve = 0;
		// Runs parallel simulations 
		for (int trialNum = 0; trialNum < numTrials; trialNum ++) {
			start = omp_get_wtime();
			paraInt(subInt, numThreads);
			intAve += omp_get_wtime() - start;
		}
		printf("Number of threads:\t%d\tInterval Time:\t%f\n", numThreads, intAve / numTrials);
	}
}

void seq(double N) {
	double x, y, Pi = 0;
	for (double i = 1; i <= N; i ++) {
		x = (1 / N) * (i - 0.5);
		y = sqrt(1 - pow(x, 2));
		Pi += 4 * (y / N);
	}
}

void paraInt(double N, int numThreads) {
	double Pi = 0;
	#pragma omp parallel num_threads(numThreads)
	{
		double x, y;
		#pragma omp for reduction(+:Pi) schedule(auto)
		for (int i = 1; i <= (int) N; i ++) {
			x = (1 / N) * (i - 0.5);
			y = sqrt(1 - pow(x, 2));
			Pi += 4 * (y / N);
		}
	}
}

//Page 68 - calculation of pi through integration - figure 3.15

// PAge 69 - figure 3.17 - Monte Carlo
