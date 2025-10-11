#pragma once

#include <array>
#include <vector>
#include <cstdint>

#include "buffers.h"
#include "filters.h"
#include "params.h"

namespace kbeyond::dsp {

namespace filters = kbeyond::dsp;

struct FdnState {
    std::array<DelayLine, kFdnSize> lines;
    std::array<double, kFdnSize> lengths{};
    std::array<double, kFdnSize> reads{};
    std::array<double, kFdnSize> outputs{};
    std::array<double, kFdnSize> feedback{};
};

void setup_fdn(FdnState &state,
               double sr,
               double size,
               double modDepth,
               const std::vector<double> &latePattern);

void read_lines(FdnState &state,
                double modDepth,
                std::array<double, kFdnSize> &phase,
                const std::array<double, kFdnSize> &phaseInc,
                std::array<filters::Tilt, kFdnSize> &tilts,
                std::array<filters::OnePoleLP, kFdnSize> &lp,
                std::array<double, kFdnSize> &vec);

void write_feedback(FdnState &state,
                    double midIn,
                    double sideIn,
                    uint32_t &rng,
                    const std::array<double, kFdnSize> &inWeights,
                    const std::array<double, kFdnSize> &sideWeights,
                    const std::array<double, kFdnSize> &midWeights,
                    const std::array<double, kFdnSize> &decay,
                    std::array<filters::OnePoleLP, kFdnSize> &low,
                    std::array<filters::OnePoleLP, kFdnSize> &high,
                    double dampLF,
                    double dampMF,
                    double dampHF,
                    double &tailMid,
                    double &tailSide);

} // namespace kbeyond::dsp

