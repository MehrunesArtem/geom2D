#ifndef GEOM2D_GVECTOR_H
#define GEOM2D_GVECTOR_H

#include <utility>

typedef long double dbl;

extern const dbl PI;
extern const dbl EPS;

class GMatrix;

class GVector {
public:
    dbl x;
    dbl y;

    GVector(const dbl& x = 0.0L, const dbl& y = 0.0L);

    static GVector zero();
    bool isZero() const;
    GVector normal() const;
    dbl norm() const;
    dbl angle() const;
    GVector unit() const;
    bool isUnit() const;
    static GVector fromPolar(const dbl& r = 0.0L, const dbl& a = 0.0L);
    static std::pair<dbl, dbl> toPolar(const dbl& x = 0.0L, const dbl& y = 0.0L);
    std::pair<dbl, dbl> toPolar() const;
    dbl const& first() const;
    dbl& first();
    dbl const& second() const;
    dbl& second();
    GVector& operator+=(const GVector&);
    GVector& operator-=(const GVector&);
    GVector& operator*=(const dbl&);
    GVector& operator/=(const dbl&);
    GVector operator-() const;
    dbl dot(const GVector&) const;
    dbl cross(const GVector&) const;
    dbl angle(const GVector&) const;
    bool opposite(const GVector&) const;
    bool parallel(const GVector&) const;
    bool antiparallel(const GVector&) const;
    bool collinear(const GVector&) const;
    bool orthogonal(const GVector&) const;
    GVector transform(const GMatrix&) const;
    GVector rotateCCW(const dbl&) const;
    GVector rotateCW(const dbl&) const;
};

GVector operator+(const GVector&, const GVector&);
GVector operator-(const GVector&, const GVector&);
GVector operator*(const GVector&, const dbl&);
GVector operator*(const dbl&, const GVector&);
GVector operator/(const GVector&, const dbl&);
bool operator==(const GVector&, const GVector&);
bool operator!=(const GVector&, const GVector&);
bool operator<(const GVector&, const GVector&);
bool operator>=(const GVector&, const GVector&);
bool operator>(const GVector&, const GVector&);
bool operator<=(const GVector&, const GVector&);
dbl dot(const GVector&, const GVector&);
dbl cross(const GVector&, const GVector&);
dbl angle(const GVector&, const GVector&);
bool opposite(const GVector&, const GVector&);
bool parallel(const GVector&, const GVector&);
bool antiparallel(const GVector&, const GVector&);
bool collinear(const GVector&, const GVector&);
bool orthogonal(const GVector&, const GVector&);
GVector operator*(const GMatrix&, const GVector&);

#endif //GEOM2D_GVECTOR_H
