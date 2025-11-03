#include <iostream>
#include <stdexcept>

using namespace std;

long long gcd(long long a, long long b) {
    while (b != 0) {
        long long t = b;
        b = a % b;
        a = t;
    }
    return a < 0 ? -a : a;
}

// Rational number class
class Rational {
private:
    long long num;  // numerator
    long long den;  // denominator

    void normalize() {
        if (den == 0) throw runtime_error("Division by zero in Rational");
        if (den < 0) { num = -num; den = -den; }
        long long g = gcd(abs(num), abs(den));
        num /= g;
        den /= g;
    }

public:
    Rational(long long n = 0, long long d = 1) : num(n), den(d) {
        normalize();
    }

    Rational operator+(const Rational& r) const {
        return Rational(num * r.den + r.num * den, den * r.den);
    }

    Rational operator-(const Rational& r) const {
        return Rational(num * r.den - r.num * den, den * r.den);
    }

    Rational operator*(const Rational& r) const {
        return Rational(num * r.num, den * r.den);
    }

    Rational operator/(const Rational& r) const {
        return Rational(num * r.den, den * r.num);
    }

    friend ostream& operator<<(ostream& os, const Rational& r) {
        if (r.den == 1) os << r.num;
        else os << r.num << "/" << r.den;
        return os;
    }
};

// Matrix class
class Matrix {
private:
    size_t rows, cols;
    Rational** data;

public:
    // Constructors
    Matrix(size_t r = 0, size_t c = 0) : rows(r), cols(c) {
        data = new Rational * [rows];
        for (size_t i = 0; i < rows; ++i)
            data[i] = new Rational[cols];
    }

    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols) {
        data = new Rational * [rows];
        for (size_t i = 0; i < rows; ++i) {
            data[i] = new Rational[cols];
            for (size_t j = 0; j < cols; ++j)
                data[i][j] = other.data[i][j];
        }
    }

    ~Matrix() {
        for (size_t i = 0; i < rows; ++i) delete[] data[i];
        delete[] data;
    }

    Rational* operator[](size_t i) { return data[i]; }
    const Rational* operator[](size_t i) const { return data[i]; }

    size_t rowCount() const { return rows; }
    size_t colCount() const { return cols; }

    // Input matrix elements
    void input() {
        cout << "Enter matrix elements (" << rows << "x" << cols << "):\n";
        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < cols; ++j) {
                long long a, b;
                cout << "Element [" << i << "][" << j << "] numerator: ";
                cin >> a;
                cout << "Element [" << i << "][" << j << "] denominator (1 if integer): ";
                cin >> b;
                data[i][j] = Rational(a, b);
            }
    }

    // Print matrix
    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j)
                cout << data[i][j] << "\t";
            cout << "\n";
        }
    }

    // Basic matrix operations
    static Matrix add(const Matrix& A, const Matrix& B) {
        if (A.rows != B.rows || A.cols != B.cols)
            throw runtime_error("Matrix dimensions do not match for addition.");
        Matrix C(A.rows, A.cols);
        for (size_t i = 0; i < A.rows; ++i)
            for (size_t j = 0; j < A.cols; ++j)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    static Matrix subtract(const Matrix& A, const Matrix& B) {
        if (A.rows != B.rows || A.cols != B.cols)
            throw runtime_error("Matrix dimensions do not match for subtraction.");
        Matrix C(A.rows, A.cols);
        for (size_t i = 0; i < A.rows; ++i)
            for (size_t j = 0; j < A.cols; ++j)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    static Matrix classicalMultiply(const Matrix& A, const Matrix& B) {
        if (A.colCount() != B.rowCount())
            throw runtime_error("Matrix dimensions do not match for multiplication.");
        Matrix C(A.rowCount(), B.colCount());
        for (size_t i = 0; i < A.rowCount(); ++i)
            for (size_t j = 0; j < B.colCount(); ++j) {
                C[i][j] = Rational(0);
                for (size_t k = 0; k < A.colCount(); ++k)
                    C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        return C;
    }
};

// Strassen's algorithm (for square matrices only)
Matrix strassenMultiply(const Matrix& A, const Matrix& B) {
    size_t n = A.rowCount();
    if (A.rowCount() != A.colCount() || B.rowCount() != B.colCount() || A.rowCount() != B.rowCount()) {
        throw runtime_error("Strassen's algorithm supports only square matrices of the same size.");
    }

    if (n <= 2) return Matrix::classicalMultiply(A, B);

    size_t k = n / 2;
    Matrix A11(k, k), A12(k, k), A21(k, k), A22(k, k);
    Matrix B11(k, k), B12(k, k), B21(k, k), B22(k, k);

    for (size_t i = 0; i < k; ++i)
        for (size_t j = 0; j < k; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];
            A22[i][j] = A[i + k][j + k];
            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];
            B22[i][j] = B[i + k][j + k];
        }

    Matrix M1 = strassenMultiply(Matrix::add(A11, A22), Matrix::add(B11, B22));
    Matrix M2 = strassenMultiply(Matrix::add(A21, A22), B11);
    Matrix M3 = strassenMultiply(A11, Matrix::subtract(B12, B22));
    Matrix M4 = strassenMultiply(A22, Matrix::subtract(B21, B11));
    Matrix M5 = strassenMultiply(Matrix::add(A11, A12), B22);
    Matrix M6 = strassenMultiply(Matrix::subtract(A21, A11), Matrix::add(B11, B12));
    Matrix M7 = strassenMultiply(Matrix::subtract(A12, A22), Matrix::add(B21, B22));

    Matrix C11 = Matrix::add(Matrix::subtract(Matrix::add(M1, M4), M5), M7);
    Matrix C12 = Matrix::add(M3, M5);
    Matrix C21 = Matrix::add(M2, M4);
    Matrix C22 = Matrix::add(Matrix::subtract(Matrix::add(M1, M3), M2), M6);

    Matrix C(n, n);
    for (size_t i = 0; i < k; ++i)
        for (size_t j = 0; j < k; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + k] = C12[i][j];
            C[i + k][j] = C21[i][j];
            C[i + k][j + k] = C22[i][j];
        }

    return C;
}

// Main function
int main() {
    cout << "Strassen Algorithm for Rational Matrices (T**)\n\n";

    size_t rowsA, colsA, rowsB, colsB;
    cout << "Enter dimensions for Matrix A (rows cols): ";
    cin >> rowsA >> colsA;
    cout << "Enter dimensions for Matrix B (rows cols): ";
    cin >> rowsB >> colsB;

    Matrix A(rowsA, colsA);
    Matrix B(rowsB, colsB);

    cout << "\nInput Matrix A\n";
    A.input();

    cout << "\nInput Matrix B\n";
    B.input();

    cout << "\nMatrix A:\n";
    A.print();
    cout << "\nMatrix B:\n";
    B.print();

    Matrix C = Matrix::classicalMultiply(A, B);

    cout << "\nResult (Classical Multiplication):\n";
    C.print();

    if (rowsA == colsA && rowsA == rowsB && rowsB == colsB && (rowsA & (rowsA - 1)) == 0) {
        cout << "\nResult (Strassen Multiplication for square matrices):\n";
        Matrix S = strassenMultiply(A, B);
        S.print();
    } else {
        cout << "\nStrassen multiplication skipped (requires square matrices of power of 2 size).\n";
    }

    // Future development:
    // 1. Implement inverse matrix using Gauss-Jordan method
    // 2. Implement inverse matrix using LU decomposition
    // 3. Implement linear regression using given data and labels

    return 0;
}
