#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

template <size_t DIM, typename T>
struct vec {
    vec() {
        for (size_t i = DIM; i--; data_[i] = T())
            ;
    }
    T &operator[](const size_t i) {
        assert(i < DIM);
        return data_[i];
    }
    const T &operator[](const size_t i) const {
        assert(i < DIM);
        return data_[i];
    }

private:
    T data_[DIM];
};

typedef vec<2, double> Vec2f;
typedef vec<3, double> Vec3f;
typedef vec<3, int> Vec3i;
typedef vec<4, double> Vec4f;

template <typename T>
struct vec<2, T> {
    vec() : x(T()), y(T()) {}
    vec(T X, T Y) : x(X), y(Y) {}
    template <class U>
    vec<2, T>(const vec<2, U> &v);
    T &operator[](const size_t i) {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    const T &operator[](const size_t i) const {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    T x, y;
};

template <typename T>
struct vec<3, T> {
    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
    T &operator[](const size_t i) {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }
    const T &operator[](const size_t i) const {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }
    double norm() const { return std::sqrt(x * x + y * y + z * z); }
    vec<3, T> &normalize(T l = 1) {
        *this = (*this) * (l / norm());
        return *this;
    }

    vec<3, T> normalize(T l = 1) const {
        return (*this) * (l / norm());
    }

    vec<3, T> cross(const vec<3, T> &vector) const {
        return vec<3, T>(
            y * vector.z - z * vector.y,
            z * vector.x - x * vector.z,
            x * vector.y - y * vector.x);
    }
    T x, y, z;
};

template <typename T>
struct vec<4, T> {
    vec() : x(T()), y(T()), z(T()), w(T()) {}
    vec(T X, T Y, T Z, T W) : x(X), y(Y), z(Z), w(W) {}
    T &operator[](const size_t i) {
        assert(i < 4);
        return i <= 0 ? x : (1 == i ? y : (2 == i ? z : w));
    }
    const T &operator[](const size_t i) const {
        assert(i < 4);
        return i <= 0 ? x : (1 == i ? y : (2 == i ? z : w));
    }
    T x, y, z, w;
};

template <size_t DIM, typename T>
T operator*(const vec<DIM, T> &lhs, const vec<DIM, T> &rhs) {
    T ret = T();
    for (size_t i = DIM; i--; ret += lhs[i] * rhs[i])
        ;
    return ret;
}

template <size_t DIM, typename T>
vec<DIM, T> operator+(vec<DIM, T> lhs, const vec<DIM, T> &rhs) {
    for (size_t i = DIM; i--; lhs[i] += rhs[i])
        ;
    return lhs;
}

template <size_t DIM, typename T>
vec<DIM, T> operator-(vec<DIM, T> lhs, const vec<DIM, T> &rhs) {
    for (size_t i = DIM; i--; lhs[i] -= rhs[i])
        ;
    return lhs;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator*(const vec<DIM, T> &lhs, const U &rhs) {
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = lhs[i] * rhs)
        ;
    return ret;
}

template <size_t DIM, typename T>
vec<DIM, T> operator-(const vec<DIM, T> &lhs) {
    return lhs * T(-1);
}

template <typename T>
vec<3, T> cross(vec<3, T> v1, vec<3, T> v2) {
    return vec<3, T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

template <size_t DIM, typename T>
std::ostream &operator<<(std::ostream &out, const vec<DIM, T> &v) {
    for (unsigned int i = 0; i < DIM; i++) {
        out << v[i] << " ";
    }
    return out;
}
#endif //__GEOMETRY_H__