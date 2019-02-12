/*
Jennifer Cho + Alex Cadigan
2/6/2019
COMP-481 Final Project
*/

#include <cmath>
#include <iostream>
#include <omp.h>

using namespace std;

void seq(double N);
void paraInterval(double N, int numThreads);
void paraIntegral(int intervals, int numThreads);
void paraMonteCarlo(int num_shots, int numThreads);
double rnd(unsigned int * seed);

int main() {
	// Gets user input:
	int numTrials, numCalc;
	cout << "Enter the number of trials to run:\t";
	cin >> numTrials;
	cout << "Enter the number of calculation iterations to use:\t";
	cin >> numCalc;

	// Stores the average times
	double start, seqAve = 0, intervalAve, integralAve, monteCarloAve;

	// Runs sequential simulations
	for (int trialNum = 0; trialNum < numTrials; trialNum ++) {
		start = omp_get_wtime();
		seq(numCalc);
		seqAve += omp_get_wtime() - start;
	}
	printf("Sequential Time:\t%f\n", seqAve / numTrials);

	// Runs through the different thread amounts
	for (int numThreads = 2; numThreads <= 64; numThreads *= 2) {
		intervalAve = 0, integralAve = 0, monteCarloAve = 0;
		// Runs parallel simulations 
		for (int trialNum = 0; trialNum < numTrials; trialNum ++) {
			start = omp_get_wtime();
			paraInterval(numCalc, numThreads);
			intervalAve += omp_get_wtime() - start;
			start = omp_get_wtime();
			paraIntegral(numCalc, numThreads);
			integralAve += omp_get_wtime() - start;
			start = omp_get_wtime();
			paraMonteCarlo(numCalc, numThreads);
			monteCarloAve += omp_get_wtime() - start;
		}
		printf("Number of threads:\t%d\tInterval Time:\t%f\tIntegral Time:\t%f\tMonte Carlo Time:\t%f\n", numThreads, intervalAve / numTrials, integralAve / numTrials, monteCarloAve / numTrials);
	}
}

/*
Runs sequential algorithm
*/
void seq(double N) {
	double x, y, Pi = 0;
	for (double i = 1; i <= N; i ++) {
		x = (1 / N) * (i - 0.5);
		y = sqrt(1 - pow(x, 2));
		Pi += 4 * (y / N);
	}
}

/*
Runs parallel interval algorithm
*/
void paraInterval(double N, int numThreads) {
	double Pi = 0;
	#pragma omp parallel num_threads(numThreads)
	{
		double x, y;
		#pragma omp for reduction(+:Pi)
		for (int i = 1; i <= (int) N; i ++) {
			x = (1 / N) * (i - 0.5);
			y = sqrt(1 - pow(x, 2));
			Pi += 4 * (y / N);
		}
	}
}

/*
Runs parallel integral algorithm
*/
void paraIntegral(int intervals, int numThreads) {
	double integral = 0;
	double dx = 1 / intervals;
	#pragma omp parallel for num_threads(numThreads) reduction(+:integral) 
	for (int i = 0; i < intervals; i ++) {
		double x = i * dx;
		double fx = sqrt(1 - x * x);
		integral = integral + fx * dx;
	}
	double Pi = 4 * integral;
}

/*
Runs parallel Monte Carlo algorithm
*/
void paraMonteCarlo(int num_shots, int numThreads) {
	unsigned int seeds[numThreads];
	for (int thread = 0; thread < numThreads; thread ++) {
		seeds[thread] = thread;
	}
	int num_hits = 0;
	#pragma omp parallel for num_threads(numThreads) reduction(+:num_hits)
	for (int shot = 0; shot < num_shots; shot ++) {
		int thread = omp_get_thread_num();
		double x = rnd(& seeds[thread]);
		double y = rnd(& seeds[thread]);
		if (x * x + y * y <= 1) {
			num_hits = num_hits + 1;
		}
	}
	double Pi = 4 * (double) num_hits / (double) num_shots;
}

double rnd(unsigned int * seed) {
	* seed = (1140671485 * (* seed) + 12820163) % (1 << 24);
	return ((double) ( * seed)) / (1 << 24);
}
