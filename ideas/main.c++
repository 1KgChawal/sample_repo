#include <complex>
#include <iostream>
#include <vector>

std::vector<std::vector<std::complex<long double>>> operator*(const std::vector<std::vector<std::complex<long double>>>& a, const std::vector<std::vector<std::complex<long double>>>& b) {
    if (a[0].size() != b.size()) {
        std::cerr << "ERROR: Matrix not multiplicable" << std::endl;
    }
    std::vector<std::vector<std::complex<long double>>> v(a.size());
    for (int i = 0; i < v.size(); i++) {
        v[i].resize(b[0].size(), std::complex<long double>(0, 0));
    }
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            for (int k = 0; k < a[0].size(); k++) {
                v[i][j] = v[i][j] + ((a[i][k]) * (b[k][j]));
            }
        }
    }
    return v;
}

std::complex<long double> _2x2_det_(const std::vector<std::vector<std::complex<long double>>>& v) {
    return v[0][0] * v[1][1] - v[0][1] * v[1][0];
}

std::vector<std::vector<std::complex<long double>>> operator*(const std::complex<long double>& a, const std::vector<std::vector<std::complex<long double>>>& b) {
    std::vector<std::vector<std::complex<long double>>> v(b.size());
    for (int i = 0; i < b.size(); i++) {
        v[i].resize(b[i].size());
        for (int j = 0; j < b[i].size(); j++) {
            v[i][j] = b[i][j] * a;
        }
    }
    return v;
}

std::vector<std::vector<std::complex<long double>>> _2x2_inverse_(const std::vector<std::vector<std::complex<long double>>>& a) {
    std::vector<std::vector<std::complex<long double>>> v = a;
    swap(v[0][0], v[1][1]);
    v[0][1] = -v[0][1];
    v[1][0] = -v[1][0];
    if (_2x2_det_(v) == std::complex<long double>(0, 0)) {
        std::cerr << "ERROR: Zero determinant" << std::endl;
        exit(1);
    }
    v = (std::complex<long double>(1, 0) / _2x2_det_(v)) * v;
    return v;
}

class _2nd_order_linear_homogeneous_ {
   private:
    int a, b, c;
    std::pair<int, int> f1, f2;
    std::vector<std::vector<std::complex<long double>>> v, M, A;
    std::complex<long double> a_, b_;

   public:
    _2nd_order_linear_homogeneous_(int a, int b, int c, int n1, int f1, int n2, int f2) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->f1.first = n1;
        this->f1.second = f1;
        this->f2.first = n2;
        this->f2.second = f2;
        M.resize(2);
        for (int i = 0; i < M.size(); i++) {
            M[i].resize(2);
        }
        v.resize(2);
        for (int i = 0; i < v.size(); i++) {
            v[i].resize(1);
        }
        A.resize(2);
        for (int i = 0; i < v.size(); i++) {
            A[i].resize(1);
        }
        A[0][0] = std::complex<long double>(this->f1.second, 0);
        A[1][0] = std::complex<long double>(this->f2.second, 0);
    }
    void solve() {
        a_ = std::complex<long double>(b * b - 4 * a * c, 0);
        b_ = std::complex<long double>(b * b - 4 * a * c, 0);
        a_ = sqrt(a_);
        b_ = sqrt(a_);
        a_ = a_ + std::complex<long double>(-b, 0);
        b_ = std::complex<long double>(-b, 0) - b_;
        a_ /= std::complex<long double>(2 * a, 0);
        b_ /= std::complex<long double>(2 * a, 0);
        if (b * b != 4 * a * c) {
            M[0][0] = pow(a_, f1.first);
            M[0][1] = pow(b_, f1.first);
            M[1][0] = pow(a_, f2.first);
            M[1][1] = pow(b_, f2.first);
        } else {
            M[0][0] = pow(a_, f1.first);
            M[0][1] = M[0][0] * std::complex<long double>(f1.first, 0);
            M[1][0] = pow(a_, f2.first);
            M[1][1] = M[1][0] * std::complex<long double>(f2.first, 0);
        }
        M = _2x2_inverse_(M);
        v = M * A;
    }
    void display() {
        if (b * b != 4 * a * c) {
            std::cout << "a_n = ( " << std::real(a_) << " + " << std::imag(a_) << "i )^n*( " << std::real(v[0][0]) << " + " << std::imag(v[0][0]) << "i ) + ( " << std::real(b_) << " + " << std::imag(b_) << "i )^n*( " << std::real(v[1][0]) << " + " << std::imag(v[1][0]) << "i )" << std::endl;
        } else {
            std::cout << "( ( " << std::real(v[0][0]) << " + " << std::imag(v[0][0]) << "i ) + n*( " << std::real(v[1][0]) << " + " << std::imag(v[1][0]) << "i ) )*( " << std::real(a_) << " + " << std::imag(a_) << "i )^n" << std::endl;
        }
    }
};

int main() {
    _2nd_order_linear_homogeneous_ a(1, -5, 6, 0, 2, 1, 5);
    a.solve();
    a.display();
}