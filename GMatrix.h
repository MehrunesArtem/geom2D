#ifndef GEOM2D_GMATRIX_H
#define GEOM2D_GMATRIX_H

#include "GVector.h"
#include <cstddef>
#include <cinttypes>

typedef long double dbl;

extern const dbl EPS;

class GMatrix {
    friend GMatrix operator+(const GMatrix&, const GMatrix&);
    friend GMatrix operator-(const GMatrix&, const GMatrix&);
    friend GMatrix operator*(const GMatrix&, const GMatrix&);
    friend GMatrix operator/(const GMatrix&, const GMatrix&);
    friend GMatrix operator*(const GMatrix&, const dbl&);
    friend GMatrix operator*(const dbl&, const GMatrix&);
    friend GMatrix operator/(const GMatrix&, const dbl&);
    friend GVector operator*(const GMatrix&, const GVector&);
    friend bool operator==(const GMatrix&, const GMatrix&);
    friend bool operator!=(const GMatrix&, const GMatrix&);
    friend bool operator<(const GMatrix&, const GMatrix&);
    friend bool operator>=(const GMatrix&, const GMatrix&);
    friend bool operator>(const GMatrix&, const GMatrix&);
    friend bool operator<=(const GMatrix&, const GMatrix&);
public:
    GMatrix(const dbl& k = 0.0L);
    GMatrix(const GVector&, const GVector&);
    GMatrix(const dbl&, const dbl&, const dbl&, const dbl&);

    static GMatrix zero();
    bool isZero() const;
    static GMatrix identity();
    bool isIdentity() const;
    dbl const& operator()(const size_t& r, const size_t& c) const;
    dbl& operator()(const size_t& r, const size_t& c);
    GMatrix& operator+=(const GMatrix&);
    GMatrix& operator-=(const GMatrix&);
    GMatrix& operator*=(const GMatrix&);
    GMatrix& operator/=(const GMatrix&);
    GMatrix& operator*=(const dbl& k);
    GMatrix& operator/=(const dbl& k);
    GMatrix operator-() const;
    GVector const& first() const;
    GVector& first();
    GVector const& second() const;
    GVector& second();
    GMatrix adjoint() const;
    GMatrix transpose() const;
    dbl determinant() const;
    int rank() const;
    dbl trace() const;
    GMatrix inverse() const;
    GMatrix pow(const int64_t&) const;
    bool invertible() const;
    bool normal() const;
    bool orthogonal() const;
    bool diagonal() const;
    bool antidiagonal() const;
    bool symmetric() const;
    bool persymmetric() const;
    bool bisymmetric() const;
    bool skew_symmetric() const;
    bool involutory() const;
    bool idempotent() const;

private:
    static GMatrix bin_pow(GMatrix, int64_t, GMatrix = 1);
    GVector x;
    GVector y;
};

GMatrix operator+(const GMatrix&, const GMatrix&);
GMatrix operator-(const GMatrix&, const GMatrix&);
GMatrix operator*(const GMatrix&, const GMatrix&);
GMatrix operator/(const GMatrix&, const GMatrix&);
GMatrix operator*(const GMatrix&, const dbl&);
GMatrix operator*(const dbl&, const GMatrix&);
GMatrix operator/(const GMatrix&, const dbl&);
bool operator==(const GMatrix&, const GMatrix&);
bool operator!=(const GMatrix&, const GMatrix&);
bool operator<(const GMatrix&, const GMatrix&);
bool operator>=(const GMatrix&, const GMatrix&);
bool operator>(const GMatrix&, const GMatrix&);
bool operator<=(const GMatrix&, const GMatrix&);
GMatrix pow(const GMatrix&, const int64_t&);

#endif //GEOM2D_GMATRIX_H
