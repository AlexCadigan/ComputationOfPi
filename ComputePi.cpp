/*
Jennifer Cho + Alex Cadigan
2/6/2019
COMP-481 Final Project
*/

#include <cmath>
#include <iostream>
#include <omp.h>

using namespace std;

double seq(double N);
double paraInterval(int N, int numThreads);
double paraIntegral(int intervals, int numThreads);
double paraMonteCarlo(int num_shots, int numThreads);

int main() {
	// Gets user input:
	int numTrials, numCalc;
	cout << "Enter the number of trials to run:\t";
	cin >> numTrials;
	cout << "Enter the number of calculation iterations to use:\t";
	cin >> numCalc;

	// Stores the average times
	double start, seqAve = 0, intervalAve, integralAve, monteCarloAve, PiSeq, PiInterval, PiIntegral, PiMonte;

	// Runs sequential simulations
	for (int trialNum = 0; trialNum < numTrials; trialNum ++) {
		start = omp_get_wtime();
		// PiSeq = seq(numCalc);
		seqAve += omp_get_wtime() - start;
	}
	printf("T_S:\t%f\tE:\t%f%%\n", seqAve / numTrials, (abs(M_PI - PiSeq) / M_PI) * 100);
	// Runs through the different thread amounts
	for (int numThreads = 2; numThreads <= 1024; numThreads *= 2) {
		intervalAve = 0, integralAve = 0, monteCarloAve = 0;
		// Runs parallel simulations 
		for (int trialNum = 0; trialNum < numTrials; trialNum ++) {
			start = omp_get_wtime();
			PiInterval = paraInterval(numCalc, numThreads);
			intervalAve += omp_get_wtime() - start;
			start = omp_get_wtime();
			PiIntegral = paraIntegral(numCalc, numThreads);
			printf("Pi Integral:\t%f\n", PiIntegral);
			integralAve += omp_get_wtime() - start;
			start = omp_get_wtime();
			PiMonte = paraMonteCarlo(numCalc, numThreads);
			monteCarloAve += omp_get_wtime() - start;
		}
		printf("Th:\t%d\tT_Intv:\t%f\tE:\t%f%%\tT_Intg:\t%f\tE:\t%f%%\tT_M:\t%f\tE:\t%f%%\n", numThreads, intervalAve / numTrials, (labs(M_PI - PiInterval) / M_PI) * 100, integralAve / numTrials, (labs(M_PI - PiIntegral) / M_PI) * 100, monteCarloAve / numTrials, (labs(M_PI - PiMonte) / M_PI) * 100);
	}
}

/*
Runs sequential algorithm
*/
double seq(double N) {
	double x, y, Pi = 0;
	for (double i = 1; i <= N; i ++) {
		x = (1 / N) * (i - 0.5);
		y = sqrt(1 - pow(x, 2));
		Pi += 4 * (y / N);
	}
	printf("Pi Seq:\t%f\n", Pi);
	return Pi;
}

/*
Runs parallel interval algorithm
*/
double paraInterval(int N, int numThreads) {
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
	printf("Pi Interval:\t%f\n", Pi);
	return Pi;
}

/*
Runs parallel integral algorithm
*/
double paraIntegral(int intervals, int numThreads) {
	double integral = 0;
	double dx = 1 / (double) intervals;
	#pragma omp parallel for num_threads(numThreads) reduction(+:integral) 
	for (int i = 0; i < intervals; i ++) {
		double x = i * dx;
		double fx = sqrt(1 - x * x);
		integral = integral + fx * dx;
	}
	printf("Pi Integral:\t%f\n", 4 * integral);
	return 4 * integral;
}

/*
Runs parallel Monte Carlo algorithm
*/
double paraMonteCarlo(int num_shots, int numThreads) {
	int num_hits = 0;
	#pragma omp parallel num_threads(numThreads)
	{	
		unsigned int seed = omp_get_thread_num();
		#pragma omp for reduction(+:num_hits)
		for (int shot = 0; shot < num_shots; shot ++) {
			double x = (double) rand_r(& seed) / (double) RAND_MAX;
			double y = (double) rand_r(& seed) / (double) RAND_MAX;
			if (x * x + y * y <= 1) {
				num_hits = num_hits + 1;
			}
		}
	}
	printf("Pi Monte:\t%f\n", 4 * (double) num_hits / (double) num_shots);
	return 4 * (double) num_hits / (double) num_shots;
}
