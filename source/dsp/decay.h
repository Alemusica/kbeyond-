#pragma once

#include <array>

#include "params.h"

namespace kbeyond::dsp {

struct DecayState {
    double dampLF = 1.0;
    double dampMF = 1.0;
    double dampHF = 1.0;
    std::array<double, kFdnSize> perLine{};
};

void update_decay(DecayState &state,
                  double sr,
                  double regen,
                  double decaySeconds,
                  double dampLFNorm,
                  double dampMFNorm,
                  double dampHFNorm,
                  const std::array<double, kFdnSize> &fdnLengths);

} // namespace kbeyond::dsp

