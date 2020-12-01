#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
#include <iomanip>


using namespace std;

void printMatrixA(vector<vector<double>> A) ;
vector<vector<double>> inputMatrixA(int n) ;
void printAnsX(vector<double> x);
vector<double> inputB(int n) ;
vector<vector<double>> randomMatrixA(int n) ;
vector<double> randomB(int n) ;
double Determinant(vector<vector<double>> A) ;
vector<double> Cramer(vector<vector<double>> A, vector<double> b) ;
vector<double> Gauss(vector<vector<double>> A, vector<double> b) ;
int converge(vector<vector<double>> A, vector<double> x, vector<double> b, double eps) ;
vector<double> SimpleIteration(vector<vector<double>> A, vector<double> b, double eps,int &count) ;
vector<double> Seidel(vector<vector<double>> A, vector<double> b, double eps, int &count) ;
vector<double> Relaxation(vector<vector<double>> A, vector<double> b, double eps,int &count) ;
vector<double> Jacobi(vector<vector<double>> A, vector<double> b, double eps, int& count);


