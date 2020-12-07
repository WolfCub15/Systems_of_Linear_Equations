#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
#include <iomanip>


using namespace std;

void printMatrixA(const vector<vector<double>>& A);
void printAnsX(const vector<double>& x);
vector<vector<double>> inputMatrixA(int n);
vector<double> inputB(int n);
vector<vector<double>> randomMatrixA(int n);
vector<double> randomB(int n);
double determinant(vector<vector<double>> A);
vector<double> cramer(const vector<vector<double>>& A, const vector<double>& b, double eps, int& iterations);
vector<double> gauss(vector<vector<double>> A, vector<double> b, double eps, int& iterations);
double converge(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b, double eps);
vector<double> simpleIteration(const vector<vector<double>>& A, const vector<double>& b, double eps, int& count);
vector<double> seidel(const vector<vector<double>>& A, const vector<double>& b, double eps, int& count);
vector<double> relaxation(const vector<vector<double>>& A, const vector<double>& b, double eps, int& count);
vector<double> jacobi(const vector<vector<double>>& A, const vector<double>& b, double eps, int& count);
