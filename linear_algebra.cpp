#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

const double EPS = 1e-9;

// Utility: Print matrix
void printMatrix(const Matrix& A) {
    for (const auto& row : A) {
        for (double val : row)
            cout << val << "\t";
        cout << "\n";
    }
}

// Matrix Addition
Matrix add(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = A[0].size();
    Matrix C(n, Vector(m));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

// Matrix Subtraction
Matrix subtract(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = A[0].size();
    Matrix C(n, Vector(m));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

// Matrix Multiplication
Matrix multiply(const Matrix& A, const Matrix& B) {
    int n = A.size(), p = A[0].size(), m = B[0].size();
    Matrix C(n, Vector(m, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < p; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// Transpose
Matrix transpose(const Matrix& A) {
    int n = A.size(), m = A[0].size();
    Matrix T(m, Vector(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            T[j][i] = A[i][j];
    return T;
}

// Identity Matrix
Matrix identity(int n) {
    Matrix I(n, Vector(n, 0));
    for (int i = 0; i < n; ++i)
        I[i][i] = 1;
    return I;
}

// Recursive Determinant
double determinant(const Matrix& A) {
    int n = A.size();
    if (n == 1) return A[0][0];

    double det = 0;
    for (int p = 0; p < n; ++p) {
        Matrix sub(n - 1, Vector(n - 1));
        for (int i = 1; i < n; ++i)
            for (int j = 0, col = 0; j < n; ++j)
                if (j != p)
                    sub[i - 1][col++] = A[i][j];

        det += (p % 2 == 0 ? 1 : -1) * A[0][p] * determinant(sub);
    }
    return det;
}

// Matrix Inversion using Gaussian Elimination
Matrix inverse(Matrix A) {
    int n = A.size();
    Matrix I = identity(n);

    for (int i = 0; i < n; ++i) {
        // Find pivot
        int pivot = i;
        for (int j = i + 1; j < n; ++j)
            if (fabs(A[j][i]) > fabs(A[pivot][i]))
                pivot = j;

        if (fabs(A[pivot][i]) < EPS)
            throw runtime_error("Matrix is singular");

        swap(A[i], A[pivot]);
        swap(I[i], I[pivot]);

        double div = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= div;
            I[i][j] /= div;
        }

        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = 0; k < n; ++k) {
                    A[j][k] -= factor * A[i][k];
                    I[j][k] -= factor * I[i][k];
                }
            }
        }
    }

    return I;
}

// Power Iteration (dominant eigenvalue)
pair<double, Vector> powerIteration(const Matrix& A, int maxIter = 1000) {
    int n = A.size();
    Vector b(n, 1.0), b_next(n);
    double lambda = 0.0;

    for (int iter = 0; iter < maxIter; ++iter) {
        // Multiply A * b
        for (int i = 0; i < n; ++i) {
            b_next[i] = 0;
            for (int j = 0; j < n; ++j)
                b_next[i] += A[i][j] * b[j];
        }

        // Normalize
        double norm = 0;
        for (double val : b_next) norm += val * val;
        norm = sqrt(norm);

        for (int i = 0; i < n; ++i)
            b[i] = b_next[i] / norm;

        // Rayleigh quotient
        Vector Ab(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                Ab[i] += A[i][j] * b[j];

        lambda = 0;
        for (int i = 0; i < n; ++i)
            lambda += b[i] * Ab[i];
    }

    return {lambda, b};
}

// Example usage
int main() {
    Matrix A = {{4, 2}, {3, 1}};
    Matrix B = {{1, 0}, {0, 1}};

    cout << "Matrix A:\n"; printMatrix(A);
    cout << "\nMatrix B:\n"; printMatrix(B);

    cout << "\nA + B:\n"; printMatrix(add(A, B));
    cout << "\nA * B:\n"; printMatrix(multiply(A, B));
    cout << "\nTranspose(A):\n"; printMatrix(transpose(A));
    cout << "\nDeterminant of A: " << determinant(A) << endl;

    try {
        cout << "\nInverse of A:\n"; printMatrix(inverse(A));
    } catch (runtime_error& e) {
        cout << e.what() << endl;
    }

    auto [lambda, eigenvec] = powerIteration(A);
    cout << "\nDominant Eigenvalue (Power Iteration): " << lambda << "\nEigenvector: ";
    for (double v : eigenvec) cout << v << " ";
    cout << endl;

    return 0;
}
