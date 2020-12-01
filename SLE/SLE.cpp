#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "headers.h"
using namespace std;

int main() {
	int lolo;
	cout << "Enter 1 - enter your data, 2 - random data\n";
	cin >> lolo;
	int size;
	vector<vector<double>> A;
	vector<double> b;
	double eps;
	if (lolo == 1) {
		cout << "Enter system order:\n";
		cin >> size;
		cout << "Enter the matrix coefficients:\n";
		A = inputMatrixA(size);
		cout << "Enter the vector of free terms (vector b):\n";
		b = inputB(size);
		cout << "Enter precision (epsilon):\n";
		cin >> eps;
	}
	else {
		cout << "Enter system order:\n";
		cin >> size;
		A = randomMatrixA(size);
		//printMatrixA(A);
		b = randomB(size);
		//printAnsX(b);
		eps = 0.01;
	}
	cout << "*********************************************************************************************\n";
	cout<<"Cramer's method:\n";
	double T1start = clock();
	vector<double> x1 = Cramer(A,b);
	double T1end = clock();
	printAnsX(x1);
	T1end-=T1start;
	cout << fixed << setprecision(20) << "Runtime:  " << 1000.0 * T1end / CLOCKS_PER_SEC << "ms" << '\n';

	cout << "*********************************************************************************************\n";
	cout<<"Gauss method:\n";
	double T2start = clock();
	vector<double> x2 = Gauss(A,b);
	double T2end = clock();
	printAnsX(x2);
	T2end-=T2start;
	cout << fixed << setprecision(20) << "Runtime:  " << 1000.0 * T2end / CLOCKS_PER_SEC << "ms" << '\n';

	cout << "*********************************************************************************************\n";
	cout << "Method of Simple Iteration:\n";
	int k1 = 0;
	double T3start = clock();
	vector<double> x3 = SimpleIteration(A, b,eps,k1);
	double T3end = clock();
	cout << "Number of iterations:  " << k1 << '\n';
	printAnsX(x3);
	T3end -= T3start;
	cout << fixed << setprecision(20) << "Runtime:  " << 1000.0 * T3end / CLOCKS_PER_SEC << "ms" << '\n';

	cout << "*********************************************************************************************\n";
	cout << "Seidel's method:\n";
	int k2 = 0;
	double T4start = clock();
	vector<double> x4 = Seidel(A, b, eps, k2);
	double T4end = clock();
	cout << "Number of iterations:  " << k2 << '\n';
	printAnsX(x4);
	T4end -= T4start;
	cout << fixed << setprecision(20) << "Runtime:  " << 1000.0 *T4end / CLOCKS_PER_SEC <<"ms"<< '\n';

	cout << "*********************************************************************************************\n";
	cout << "Relaxation method:\n";
	int k3 = 0;
	double T5start = clock();
	vector<double> x5 = Relaxation(A, b, eps, k3);
	double T5end = clock();
	cout << "Number of iterations:  " << k3 << '\n';
	printAnsX(x5);
	T5end -= T5start;
	cout << fixed << setprecision(20) << "Runtime:  " << 1000.0 * T5end / CLOCKS_PER_SEC << "ms" << '\n';

	cout << "*********************************************************************************************\n";
	cout << "Jacobi method:\n";
	int k4 = 0;
	double T6start = clock();
	vector<double> x6 = Jacobi(A, b, eps, k4);
	double T6end = clock();
	cout << "Number of iterations:  " << k4 << '\n';
	printAnsX(x6);
	T6end -= T6start;
	cout << fixed << setprecision(20) << "Runtime:  " << 1000.0 * T5end / CLOCKS_PER_SEC << "ms" << '\n';

	return 0;
}