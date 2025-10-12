#pragma once

#include <array>

#include "params.h"

namespace kbeyond::dsp::mixing {

void apply_walsh_hadamard16(const std::array<double, kFdnSize> &input,
                            std::array<double, kFdnSize> &output);

void apply_walsh_hadamard16_2d(const std::array<double, kFdnSize> &input,
                               std::array<double, kFdnSize> &output);

void apply_hybrid_diffusion(const std::array<double, kFdnSize> &u,
                            const std::array<double, kFdnSize> &input,
                            std::array<double, kFdnSize> &output,
                            std::array<double, kFdnSize> &scratch);

} // namespace kbeyond::dsp::mixing

