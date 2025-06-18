# Matrix Operations Library

A comprehensive C++ library for basic matrix operations including arithmetic, linear algebra computations, and eigenvalue analysis.

## Features

- **Basic Matrix Operations**
  - Addition and subtraction
  - Matrix multiplication
  - Matrix transpose
  - Identity matrix generation

- **Advanced Operations**
  - Determinant calculation (recursive method)
  - Matrix inversion using Gaussian elimination
  - Dominant eigenvalue computation using Power Iteration

- **Utilities**
  - Matrix printing for debugging
  - Error handling for singular matrices
  - Configurable precision with epsilon tolerance

## Requirements

- C++11 or later
- Standard libraries: `<iostream>`, `<vector>`, `<cmath>`, `<limits>`

## Usage

### Basic Setup

```cpp
#include "matrix_lib.cpp"  // or compile and link separately

using Matrix = vector<vector<double>>;
using Vector = vector<double>;
```

### Creating Matrices

```cpp
// Create a 2x2 matrix
Matrix A = {{4, 2}, 
            {3, 1}};

// Create identity matrix
Matrix I = identity(3);  // 3x3 identity matrix
```

### Basic Operations

```cpp
Matrix A = {{4, 2}, {3, 1}};
Matrix B = {{1, 0}, {0, 1}};

// Addition
Matrix sum = add(A, B);

// Subtraction  
Matrix diff = subtract(A, B);

// Multiplication
Matrix product = multiply(A, B);

// Transpose
Matrix At = transpose(A);
```

### Advanced Operations

```cpp
// Determinant
double det = determinant(A);
cout << "Determinant: " << det << endl;

// Matrix Inversion
try {
    Matrix inv_A = inverse(A);
    cout << "Inverse found successfully" << endl;
} catch (runtime_error& e) {
    cout << "Error: " << e.what() << endl;
}

// Eigenvalue Analysis
auto [eigenvalue, eigenvector] = powerIteration(A);
cout << "Dominant eigenvalue: " << eigenvalue << endl;
```

## Function Reference

### Basic Operations

| Function | Description | Parameters | Returns |
|----------|-------------|------------|---------|
| `add(A, B)` | Matrix addition | Two matrices of same dimensions | Result matrix |
| `subtract(A, B)` | Matrix subtraction | Two matrices of same dimensions | Result matrix |
| `multiply(A, B)` | Matrix multiplication | Compatible matrices (A cols = B rows) | Result matrix |
| `transpose(A)` | Matrix transpose | Any matrix | Transposed matrix |
| `identity(n)` | Identity matrix | Size n | n×n identity matrix |

### Advanced Operations

| Function | Description | Parameters | Returns |
|----------|-------------|------------|---------|
| `determinant(A)` | Calculate determinant | Square matrix | Determinant value |
| `inverse(A)` | Matrix inversion | Non-singular square matrix | Inverse matrix |
| `powerIteration(A, maxIter)` | Find dominant eigenvalue | Square matrix, max iterations (default: 1000) | `pair<eigenvalue, eigenvector>` |

### Utilities

| Function | Description | Parameters | Returns |
|----------|-------------|------------|---------|
| `printMatrix(A)` | Print matrix to console | Any matrix | void |

## Error Handling

The library includes error handling for common issues:

- **Singular Matrix**: The `inverse()` function throws `runtime_error` if the matrix is singular (determinant ≈ 0)
- **Dimension Mismatch**: Ensure matrices have compatible dimensions for operations
- **Numerical Precision**: Uses `EPS = 1e-9` for floating-point comparisons

## Example Output

```
Matrix A:
4    2    
3    1    

Matrix B:
1    0    
0    1    

A + B:
5    2    
3    2    

A * B:
4    2    
3    1    

Transpose(A):
4    3    
2    1    

Determinant of A: -2

Inverse of A:
-0.5    1    
1.5    -2    

Dominant Eigenvalue (Power Iteration): 5.37228
Eigenvector: 0.816497 0.577350
```

## Compilation

```bash
g++ -std=c++11 -o matrix_operations matrix_lib.cpp
./matrix_operations
```

## Limitations

- **Memory**: Uses `vector<vector<double>>` which may not be optimal for very large sparse matrices
- **Performance**: Recursive determinant calculation has O(n!) complexity
- **Eigenvalues**: Power iteration only finds the dominant eigenvalue
- **Precision**: Limited by double-precision floating-point arithmetic

## Potential Improvements

- Add LU decomposition for better performance
- Implement QR algorithm for all eigenvalues
- Add support for complex matrices
- Optimize memory layout for cache efficiency
- Add matrix decomposition methods (SVD, Cholesky)

## License

This code is provided as-is for educational and research purposes.
