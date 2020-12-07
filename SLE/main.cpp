#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "methods.h"
using namespace std;

#define solveTask(algoName, func, A, b, eps)                                                                    \
    cout << "*********************************************************************************************\n";  \
    cout << algoName << ":\n";                                                                                  \
    iterations = 0;                                                                                             \
    tStart = clock();                                                                                           \
    x = func(A, b, eps, iterations);                                                                            \
    tEnd = clock();                                                                                             \
    cout << "Number of iterations:  " << iterations << '\n';                                                    \
    printAnsX(x);                                                                                               \
    tEnd -= tStart;                                                                                             \
    cout << fixed << setprecision(20) << "Mean squared error: " << converge(A, x, b, eps) << '\n';              \
    cout << fixed << setprecision(20) << "Runtime:            " << tEnd / CLOCKS_PER_SEC << "s\n";

int main() {
    int lolo;
    cout << "Enter 1 - enter your data, 2 - random data\n";
    cin >> lolo;
    int size;
    vector<vector<double>> A;
    vector<double> b;
    double eps;
    if (lolo == 1) {
        cout << "Enter size of system :\n";
        cin >> size;
        cout << "Enter the matrix coefficients:\n";
        A = inputMatrixA(size);
        cout << "Enter the vector of free terms (vector b):\n";
        b = inputB(size);
        cout << "Enter precision (epsilon):\n";
        cin >> eps;
    }
    else {
        cout << "Enter size of system :\n";
        cin >> size;
        A = randomMatrixA(size);
        // printMatrixA(A);
        b = randomB(size);
        // printAnsX(b);
        eps = 0.01;
    }
    int iterations;
    vector<double> x;
    double tStart, tEnd;
    solveTask("Cramer's method",            cramer,             A, b, eps);
    solveTask("Gauss method",               gauss,              A, b, eps);
    solveTask("Method of Simple Iteration", simpleIteration,    A, b, eps);
    solveTask("Seidel's method",            seidel,             A, b, eps);
    solveTask("Relaxation method",          relaxation,         A, b, eps);
    solveTask("Jacobi method",              jacobi,             A, b, eps);
    return 0;
}
