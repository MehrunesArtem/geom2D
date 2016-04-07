#include "GVector.h"
#include <cmath>

using namespace std;

GVector::GVector(const dbl& x, const dbl& y) : x(x), y(y) { }

GVector GVector::zero() {
    return GVector(0.0L, 0.0L);
}

bool GVector::isZero() const {
    return abs(x) < EPS && abs(y) < EPS;
}

GVector GVector::normal() const {
    return GVector(-y, x);
}

dbl GVector::norm() const {
    return sqrt(x * x + y * y);
}

dbl GVector::angle() const {
    return atan2(y, x);
}

GVector GVector::unit() const {
    return GVector(x / norm(), y / norm());
}

bool GVector::isUnit() const {
    return abs(norm() - 1.0L) < EPS;
}

GVector GVector::fromPolar(const dbl& r, const dbl& a) {
    return GVector(r * cos(a), r * sin(a));
}

std::pair<dbl, dbl> GVector::toPolar(const dbl& x, const dbl& y) {
    return make_pair(sqrt(x * x + y * y), atan2(y, x));
}

std::pair<dbl, dbl> GVector::toPolar() const {
    return make_pair(norm(), angle());
}

dbl const& GVector::first() const {
    return x;
}

dbl& GVector::first() {
    return x;
}

dbl const& GVector::second() const {
    return y;
}

dbl& GVector::second() {
    return y;
}

GVector& GVector::operator+=(const GVector& v) {
    x += v.x;
    y += v.y;
    return *this;
}

GVector& GVector::operator-=(const GVector& v) {
    x -= v.x;
    y -= v.y;
    return *this;
}

GVector& GVector::operator*=(const dbl& k) {
    x *= k;
    y *= k;
    return *this;
}

GVector& GVector::operator/=(const dbl& k) {
    x /= k;
    y /= k;
    return *this;
}

GVector GVector::operator-() const {
    return GVector(-x, -y);
}

dbl GVector::dot(const GVector& v) const {
    return x * v.x + y * v.y;
}

dbl GVector::cross(const GVector& v) const {
    return x * v.y - y * v.x;
}

dbl GVector::angle(const GVector& v) const {
    return atan2(cross(v), dot(v));
}

bool GVector::opposite(const GVector& v) const {
    return abs(x + v.x) < EPS && abs(y + v.y) < EPS;
}

bool GVector::parallel(const GVector& v) const {
    return isZero() || v.isZero() || abs(angle(v)) < EPS;
}

bool GVector::antiparallel(const GVector& v) const {
    return abs(PI - abs(angle(v))) < EPS;
}

bool GVector::collinear(const GVector& v) const {
    return parallel(v) || antiparallel(v);
}

bool GVector::orthogonal(const GVector& v) const {
    return collinear(v.normal());
}

GVector GVector::transform(const GMatrix& m) const {
    return m * *this;
}

GVector GVector::rotateCCW(const dbl& a) const {
    return GVector(cos(a) * x - sin(a) * y, sin(a) * x + cos(a) * y);
}

GVector GVector::rotateCW(const dbl& a) const {
    return GVector(cos(a) * x + sin(a) * y, -sin(a) * x + cos(a) * y);
}

inline GVector operator+(const GVector& v1, const GVector& v2) {
    return GVector(v1.x + v2.x, v1.y + v2.y);
}

inline GVector operator-(const GVector& v1, const GVector& v2) {
    return GVector(v1.x - v2.x, v1.y - v2.y);
}

inline GVector operator*(const GVector& v, const dbl& k) {
    return GVector(v.x * k, v.y * k);
}

inline GVector operator*(const dbl& k, const GVector& v) {
    return v * k;
}

inline GVector operator/(const GVector& v, const dbl& k) {
    return GVector(v.x / k, v.y / k);
}

inline bool operator==(const GVector& v1, const GVector& v2) {
    return (v1 - v2).isZero();
}

inline bool operator!=(const GVector& v1, const GVector& v2) {
    return !(v1 == v2);
}

inline bool operator<(const GVector& v1, const GVector& v2) {
    return abs(v1.x - v2.y) < EPS ? v1.y + EPS < v2.y : v1.x < v2.x;
}

inline bool operator>=(const GVector& v1, const GVector& v2) {
    return !(v1 < v2);
}

inline bool operator>(const GVector& v1, const GVector& v2) {
    return v2 < v1;
}

inline bool operator<=(const GVector& v1, const GVector& v2) {
    return !(v2 < v1);
}

inline dbl dot(const GVector& v1, const GVector& v2) {
    return v1.dot(v2);
}

inline dbl cross(const GVector& v1, const GVector& v2) {
    return v1.cross(v2);
}

inline dbl angle(const GVector& v1, const GVector& v2) {
    return v1.angle(v2);
}

inline bool opposite(const GVector& v1, const GVector& v2) {
    return v1.opposite(v2);
}

inline bool parallel(const GVector& v1, const GVector& v2) {
    return v1.parallel(v2);
}

inline bool antiparallel(const GVector& v1, const GVector& v2) {
    return v1.antiparallel(v2);
}

inline bool collinear(const GVector& v1, const GVector& v2) {
    return v1.collinear(v2);
}

inline bool orthogonal(const GVector& v1, const GVector& v2) {
    return v1.orthogonal(v2);
}