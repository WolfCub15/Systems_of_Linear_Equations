#include "methods.h"

void printMatrixA(const vector<vector<double>> &A) {
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) {
            cout << A[i][j] << ' ';
        }
        cout << '\n';
    }
}

void printAnsX(const vector<double>& x){
    for (int i = 0; i < x.size(); ++i) {
        cout << fixed << setprecision(10) << "x" << i << "  :  " << x[i] << '\n';
    }
    cout << '\n';
}

vector<vector<double>> inputMatrixA(int n) {
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> A[i][j];
        }
    }
    return A;
}

vector<double> inputB(int n) {
    vector<double> v(n);
    for (int i = 0; i < n; ++i) {
        cin >> v[i];
    }
    return v;
}

vector<vector<double>> randomMatrixA(int n) {
    mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double val = gen() % 100;
            if (i == j) A[i][j] = val * 99 * n + 1000;
            else A[i][j] = val;
        }
    }
    return A;
}

vector<double> randomB(int n) {
    mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    vector<double> vec(n);
    for (int i = 0; i < n; ++i) {
        vec[i] = gen() % 100 + 99 * n;
    }
    return vec;
}

double determinant(vector<vector<double>> A) {
    double det = 1;
    const double EPS = 1e-9;
    for (int i = 0; i < A.size(); ++i) {
        int k = i;
        for (int j = i + 1; j < A.size(); ++j) {
            if (abs(A[j][i]) > abs(A[k][i]))
                k = j;
        }
        if (abs(A[k][i]) < EPS) {
            det = 0;
            break;
        }
        swap(A[i], A[k]);
        if (i != k) {
            det = -det;
        }
        det *= A[i][i];
        for (int j = i + 1; j < A.size(); ++j) {
            A[i][j] /= A[i][i];
        }
        for (int j = 0; j < A.size(); ++j) {
            if (j != i && abs(A[j][i]) > EPS) {
                for (int k = i + 1; k < A.size(); ++k) {
                    A[j][k] -= A[i][k] * A[j][i];
                }
            }
        }
    }
    return det;
}

vector<double> cramer(const vector<vector<double>>& A, const vector<double>& b, double eps, int& iterations) {
    vector<double> x;
    x.reserve(b.size());
    double det = determinant(A);
    vector<double> ds;
    ds.reserve(b.size());
    if (fabs(det) > 1e-9) {
        for (int i = 0; i < (int)A.size(); ++i) {
            vector<vector<double>> tmp = A;
            for (int j = 0; j < tmp.size(); ++j) 
                tmp[j][i] = b[j];
            double lolo = determinant(tmp);
            ds.push_back(lolo);
        }
    }
    else cout << "The system has infinitely many solutions or is incompatible\n";
    for (double di : ds) {
        double ans = di / det;
        x.push_back(ans);
    }
    return x;
}

vector<double> gauss(vector<vector<double>> A, vector<double> b, double eps, int& iterations) {
    int n = A.size();
    vector<double> x(n, 0);
    for (int q = 0; q < n; ++q) {
        for (int i = q + 1; i < n; ++i) {
            double k = A[i][q] / A[q][q];
            for (int j = q; j < n; ++j) {
                A[i][j] -= k * A[q][j];
            }
            b[i] -= k * b[q];
        }

    }
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}

double converge(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b, double eps) {
    double ans = 0;
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        ans += pow(b[i]-sum, 2);
    }
    return sqrt(ans);
}

vector<double> simpleIteration(const vector<vector<double>>& A, const vector<double>& b, double eps,int &count) {
    int n = A.size();
    vector<double> x(n,0);
    vector<double> x0;
    int max_count = 100;
    count = 0;
    do {
        for (int i = 0; i < n; ++i) {
            x0 = x;
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A[i][j] * x0[j];
                }
            }
            x[i] = (b[i]-sum) / A[i][i];
        }
        count++;
        max_count--;
    } while(converge(A, x, b, eps) >= eps && max_count);
    
    return x;
}

vector<double> seidel(const vector<vector<double>>& A, const vector<double>& b, double eps, int &count) {
    int n = A.size();
    vector<double> x(n,0),pr(n,0);
    count =0;
    int max_count = 100;
    do {
        for (int i = 0; i < n; ++i) {
            x[i] = 0;
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        count++;
        max_count--;
    } while (converge(A, x, b, eps) >= eps && max_count);
    return x;
}

vector<double> relaxation(const vector<vector<double>>& A, const vector<double>& b, double eps,int &count) {
    int n = A.size();
    vector<double> x(n, 0), xn(n, 0);
    double w;
    int max_count = 100;
    mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    w = 100 + gen() % 100;
    w /= 100;
    do {
        for (int i = 0; i < n; ++i) {
            x[i] = b[i];
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                if(i!=j) x[i]-= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
            x[i] = w * x[i] + xn[i] - xn[i] * w;
            xn[i] = x[i];
        }
        count++;
        max_count--;
    } while (converge(A, x, b, eps) >= eps && max_count);
    return x;
}

vector<double> jacobi(const vector<vector<double>>& A, const vector<double>& b, double eps, int& count) {
    int n = A.size();
    vector<double> x(n, 0);
    vector<double> x0 = b;
    int max_count = 100;
    count = 0;
    do {
        for (int i = 0; i < n; ++i) {
            x[i] = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    x[i] -= A[i][j] * x0[j];
                }
            }
            x[i] /= A[i][i];
            for (int i = 0; i < n; ++i) {
                x0[i] = x[i];
            }
        }
        count++;
        max_count--;
    } while (converge(A, x, b, eps) >= eps && max_count);
    return x;
}
