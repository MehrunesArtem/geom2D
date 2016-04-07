#include "GMatrix.h"
#include <cmath>
#include <cassert>

using namespace std;

GMatrix::GMatrix(const dbl& k) {
    x = GVector(k, 0.0L);
    y = GVector(0.0L, k);
}

GMatrix::GMatrix(const GVector& x, const GVector& y) : x(x), y(y) { }

GMatrix::GMatrix(const dbl& d1, const dbl& d2, const dbl& d3, const dbl& d4) : x(d1, d3), y(d2, d4) { }

GMatrix GMatrix::zero() {
    return GMatrix(0.0L);
}

bool GMatrix::isZero() const {
    return x.isZero() && y.isZero();
}

GMatrix GMatrix::identity(){
    return GMatrix(1.0L);
}

bool GMatrix::isIdentity() const {
    return abs(x.x - 1.0L) < EPS && abs(x.y) < EPS && abs(y.x) < EPS && abs(y.y - 1.0L) < EPS;
}

dbl const& GMatrix::operator()(const size_t& r, const size_t& c) const {
    assert(!((r | c) >> 1));
    return r ? c ? y.y : y.x : c ? x.y : x.x;
}

dbl& GMatrix::operator()(const size_t& r, const size_t& c) {
    assert(!((r | c) >> 1));
    return r ? c ? y.y : y.x : c ? x.y : x.x;
}

GMatrix& GMatrix::operator+=(const GMatrix& m) {
    x += m.x;
    y += m.y;
    return *this;
}

GMatrix& GMatrix::operator-=(const GMatrix& m) {
    x -= m.x;
    y -= m.y;
    return *this;;
}

GMatrix& GMatrix::operator*=(const GMatrix& m) {
    return *this = GMatrix(x.x * m.x.x + y.x * m.x.y, x.x * m.y.x + y.x * m.y.y,
                           x.y * m.x.x + y.y * m.x.y, x.y * m.y.x + y.y * m.y.y);
}

GMatrix& GMatrix::operator/=(const GMatrix& m) {
    return operator*=(m.inverse());
}

GMatrix& GMatrix::operator*=(const dbl& k) {
    x *= k;
    y *= k;
    return *this;
}

GMatrix& GMatrix::operator/=(const dbl& k) {
    x /= k;
    y /= k;
    return *this;
}

GMatrix GMatrix::operator-() const {
    return GMatrix(-x, -y);
}

GVector const& GMatrix::first() const {
    return x;
}

GVector& GMatrix::first() {
    return x;
}

GVector const& GMatrix::second() const {
    return y;
}

GVector& GMatrix::second() {
    return y;
}

GMatrix GMatrix::adjoint() const {
    return GMatrix(y.y, -x.y, -y.x, y.y);
}

GMatrix GMatrix::transpose() const {
    return GMatrix(x.x, x.y, y.x, y.y);
}

dbl GMatrix::determinant() const {
    return x.x * y.y - y.x * x.y;
}

int GMatrix::rank() const {
    return 2 - isZero() - !invertible();
}

dbl GMatrix::trace() const {
    return x.x * y.y;
}

GMatrix GMatrix::inverse() const {
    assert(invertible());
    return adjoint() / determinant();
}

GMatrix GMatrix::pow(const int64_t& n) const {
    return n >= 0 ? bin_pow(*this, n) : bin_pow(inverse(), -n);
}

bool GMatrix::invertible() const {
    return abs(determinant()) > EPS;
}

bool GMatrix::normal() const {
    return *this * transpose() == transpose() * *this;
}

bool GMatrix::orthogonal() const {
    return normal() && *this * transpose() == identity();
}

bool GMatrix::diagonal() const {
    return abs(x.y) < EPS && abs(y.x) < EPS;
}

bool GMatrix::antidiagonal() const {
    return abs(x.x) < EPS && abs(y.y) < EPS;
}

bool GMatrix::symmetric() const {
    return abs(x.y - y.x) < EPS;
}

bool GMatrix::persymmetric() const {
    return abs(x.x - y.y) < EPS;
}

bool GMatrix::bisymmetric() const {
    return symmetric() && persymmetric();
}

bool GMatrix::skew_symmetric() const {
    return transpose() == operator-();
}

bool GMatrix::involutory() const {
    return *this * *this == identity();
}

bool GMatrix::idempotent() const {
    return *this * *this == *this;
}

GMatrix GMatrix::bin_pow(GMatrix a, int64_t n, GMatrix r) {
    for (; n; a *= a, n >>= 1)
        if (n & 1) r *= a;
    return r;
}

inline GMatrix operator+(const GMatrix& m1, const GMatrix& m2) {
    return GMatrix(m1.x + m2.x, m1.y + m2.y);
}

inline GMatrix operator-(const GMatrix& m1, const GMatrix& m2) {
    return GMatrix(m1.x - m2.x, m1.y - m2.y);
}

inline GMatrix operator*(const GMatrix& m1, const GMatrix& m2) {
    return GMatrix(m1.x.x * m2.x.x + m1.y.x * m2.x.y, m1.x.x * m2.y.x + m1.y.x * m2.y.y,
                   m1.x.y * m2.x.x + m1.y.y * m2.x.y, m1.x.y * m2.y.x + m1.y.y * m2.y.y);
}

inline GMatrix operator/(const GMatrix& m1, const GMatrix& m2) {
    return m1 * m2.inverse();
}

inline GMatrix operator*(const GMatrix& m, const dbl& k) {
    return GMatrix(m.x * k, m.y * k);
}

inline GMatrix operator*(const dbl& k, const GMatrix& m) {
    return m * k;
}

inline GMatrix operator/(const GMatrix& m, const dbl& k) {
    return GMatrix(m.x / k, m.y / k);
}

inline bool operator==(const GMatrix& m1, const GMatrix& m2) {
    return (m1 - m2).isZero();
}

inline bool operator!=(const GMatrix& m1, const GMatrix& m2) {
    return !(m1 == m2);
}

inline bool operator<(const GMatrix& m1, const GMatrix& m2) {
    return m1.x == m2.x ? m1.y < m2.y : m1.x < m2.x;
}

inline bool operator>=(const GMatrix& m1, const GMatrix& m2) {
    return !(m1 < m2);
}

inline bool operator>(const GMatrix& m1, const GMatrix& m2) {
    return m2 < m1;
}

inline bool operator<=(const GMatrix& m1, const GMatrix& m2) {
    return !(m2 < m1);
}

inline GMatrix pow(const GMatrix& m, const int64_t& n) {
    return m.pow(n);
}

inline GVector operator*(const GMatrix& m, const GVector& v) {
    return GVector(m.x.x * v.x + m.y.x * v.y, m.x.y * v.x + m.y.y * v.y);
}