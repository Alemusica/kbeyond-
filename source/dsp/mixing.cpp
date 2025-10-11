#include "mixing.h"

#include <cmath>

#include "../householder_phi16.h"

namespace kbeyond::dsp::mixing {

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

void apply_hybrid_diffusion(const std::array<double, kFdnSize> &u,
                            const std::array<double, kFdnSize> &input,
                            std::array<double, kFdnSize> &output,
                            std::array<double, kFdnSize> &scratch) {
    apply_walsh_hadamard16(input, scratch);
    apply_householder<kFdnSize>(u, scratch, output);
}

} // namespace kbeyond::dsp::mixing

