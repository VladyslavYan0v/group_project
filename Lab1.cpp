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
// Раціональний клас
class Rational {
private:
	long long num;  // номератор
	long long den;  // деномінатор

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


// Матриці
class Matrix {
private:
    size_t n;
    Rational** data;

public:
    // Конструктори
    Matrix(size_t n_ = 0) : n(n_) {
        data = new Rational * [n];
        for (size_t i = 0; i < n; ++i)
            data[i] = new Rational[n];
    }

    Matrix(const Matrix& other) : n(other.n) {
        data = new Rational * [n];
        for (size_t i = 0; i < n; ++i) {
            data[i] = new Rational[n];
            for (size_t j = 0; j < n; ++j)
                data[i][j] = other.data[i][j];
        }
    }

    ~Matrix() {
        for (size_t i = 0; i < n; ++i) delete[] data[i];
        delete[] data;
    }

    Rational* operator[](size_t i) { return data[i]; }
    const Rational* operator[](size_t i) const { return data[i]; }

    size_t size() const { return n; }

    void input() {
        cout << "Enter matrix elements (" << n << "x" << n << "):\n";
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) {
                long long a, b;
                cout << "Element [" << i << "][" << j << "] numerator: ";
                cin >> a;
                cout << "Element [" << i << "][" << j << "] denominator (1 if integer): ";
                cin >> b;
                data[i][j] = Rational(a, b);
            }
    }

    // Виведення матриць
    void print() const {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j)
                cout << data[i][j] << "\t";
            cout << "\n";
        }
    }

	// Базові операції над матрицями
    static Matrix add(const Matrix& A, const Matrix& B) {
        size_t n = A.size();
        Matrix C(n);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    static Matrix subtract(const Matrix& A, const Matrix& B) {
        size_t n = A.size();
        Matrix C(n);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    static Matrix classicalMultiply(const Matrix& A, const Matrix& B) {
        size_t n = A.size();
        Matrix C(n);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) {
                C[i][j] = Rational(0);
                for (size_t k = 0; k < n; ++k)
                    C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        return C;
    }
};


//Алгоритм Штрассена
Matrix strassenMultiply(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    if (n <= 2) return Matrix::classicalMultiply(A, B);

    size_t k = n / 2;
    Matrix A11(k), A12(k), A21(k), A22(k);
    Matrix B11(k), B12(k), B21(k), B22(k);

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

    Matrix C(n);
    for (size_t i = 0; i < k; ++i)
        for (size_t j = 0; j < k; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + k] = C12[i][j];
            C[i + k][j] = C21[i][j];
            C[i + k][j + k] = C22[i][j];
        }

    return C;
}


// Мейн функция
int main() {
    cout << "Strassen Algorithm for Rational Matrices (T**)\n\n";

    size_t n;
    cout << "Enter the size of square matrices (power of 2: 2, 4, 8, ...): ";
    cin >> n;

    Matrix A(n);
    Matrix B(n);

    cout << "\n Input Matrix A\n";
    A.input();

    cout << "\nInput Matrix B\n";
    B.input();

    cout << "\nMatrix A:\n";
    A.print();
    cout << "\nMatrix B:\n";
    B.print();

    Matrix C = strassenMultiply(A, B);

    cout << "\nResult :\n";
    C.print();

    
    // Побудова оберненої матриці методом Гауса-Жордана;
    // Побудова оберненої матриці методом LU-розкладання;
    // Побудова лінійної регресії за відомими даними та мітками;

    return 0;
}
