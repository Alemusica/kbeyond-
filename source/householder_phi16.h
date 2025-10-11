#pragma once
// householder_phi16.h
// Utility to (implicitly) apply a Householder transform y = x - 2*u*(u^T x)
// where u is built from a φ-weighted vector, then normalized.
//
// We don't materialize the full 16x16 matrix: we keep u and compute in O(N).

#include <array>
#include <cmath>
#include <algorithm>

template <size_t N>
static inline void make_phi_vector(std::array<double, N>& u, double phiWeight) {
    const double phi = 1.6180339887498948482;
    // base: uniform vector; add φ^i trend controlled by phiWeight
    for (size_t i = 0; i < N; ++i) {
        double unif = 1.0;
        double phiPart = std::pow(phi, (double)i / (double)(N - 1)); // gentle ramp
        u[i] = (1.0 - phiWeight) * unif + (phiWeight) * phiPart;
    }
    // normalize to unit length
    double norm2 = 0.0;
    for (auto v : u) norm2 += v * v;
    norm2 = std::sqrt(std::max(1e-12, norm2));
    for (auto &v : u) v /= norm2;
}

template <size_t N>
static inline void apply_householder(const std::array<double, N>& u,
                                     const std::array<double, N>& x,
                                     std::array<double, N>& y) {
    double dot = 0.0;
    for (size_t i = 0; i < N; ++i) dot += u[i] * x[i];
    double s = -2.0 * dot;
    for (size_t i = 0; i < N; ++i) y[i] = x[i] + s * u[i];
}
