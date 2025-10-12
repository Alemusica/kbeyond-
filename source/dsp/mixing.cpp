#include "mixing.h"

#include <cmath>

#include "../householder_phi16.h"

namespace kbeyond::dsp::mixing {

namespace {
inline void hadamard4(double &a0, double &a1, double &a2, double &a3) {
    const double s0 = a0 + a1;
    const double s1 = a0 - a1;
    const double s2 = a2 + a3;
    const double s3 = a2 - a3;
    constexpr double norm = 0.5; // 1 / sqrt(4)
    a0 = (s0 + s2) * norm;
    a1 = (s1 + s3) * norm;
    a2 = (s0 - s2) * norm;
    a3 = (s1 - s3) * norm;
}
} // namespace

void apply_walsh_hadamard16(const std::array<double, kFdnSize> &input,
                            std::array<double, kFdnSize> &output) {
    std::array<double, kFdnSize> tmp = input;
    for (std::size_t len = 1; len < tmp.size(); len <<= 1) {
        const std::size_t step = len << 1;
        for (std::size_t i = 0; i < tmp.size(); i += step) {
            for (std::size_t j = 0; j < len; ++j) {
                const double a = tmp[i + j];
                const double b = tmp[i + j + len];
                tmp[i + j] = a + b;
                tmp[i + j + len] = a - b;
            }
        }
    }
    const double norm = 1.0 / std::sqrt(static_cast<double>(tmp.size()));
    for (std::size_t i = 0; i < tmp.size(); ++i)
        output[i] = tmp[i] * norm;
}

void apply_walsh_hadamard16_2d(const std::array<double, kFdnSize> &input,
                               std::array<double, kFdnSize> &output) {
    std::array<double, kFdnSize> tmp = input;
    constexpr std::size_t kSide = 4;

    for (std::size_t row = 0; row < kSide; ++row) {
        std::size_t base = row * kSide;
        double a0 = tmp[base + 0];
        double a1 = tmp[base + 1];
        double a2 = tmp[base + 2];
        double a3 = tmp[base + 3];
        hadamard4(a0, a1, a2, a3);
        tmp[base + 0] = a0;
        tmp[base + 1] = a1;
        tmp[base + 2] = a2;
        tmp[base + 3] = a3;
    }

    for (std::size_t col = 0; col < kSide; ++col) {
        double a0 = tmp[0 * kSide + col];
        double a1 = tmp[1 * kSide + col];
        double a2 = tmp[2 * kSide + col];
        double a3 = tmp[3 * kSide + col];
        hadamard4(a0, a1, a2, a3);
        tmp[0 * kSide + col] = a0;
        tmp[1 * kSide + col] = a1;
        tmp[2 * kSide + col] = a2;
        tmp[3 * kSide + col] = a3;
    }

    output = tmp;
}

void apply_hybrid_diffusion(const std::array<double, kFdnSize> &u,
                            const std::array<double, kFdnSize> &input,
                            std::array<double, kFdnSize> &output,
                            std::array<double, kFdnSize> &scratch) {
    apply_walsh_hadamard16(input, scratch);
    apply_householder<kFdnSize>(u, scratch, output);
}

} // namespace kbeyond::dsp::mixing

