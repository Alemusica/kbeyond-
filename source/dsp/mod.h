#pragma once

#include <array>

#include "params.h"

namespace kbeyond::dsp {

struct ModState {
    std::array<double, kFdnSize> fdnPhase{};
    std::array<double, kFdnSize> fdnPhaseInc{};
    std::array<double, kFdnSize> uwalkPhase{};
    std::array<double, kFdnSize> uwalkPhaseInc{};
    std::array<double, kFdnSize> uwalkState{};
    std::array<double, kFdnSize> qditherPhase{};
    std::array<double, kFdnSize> qditherPhaseInc{};
};

void init_fdn_phases(ModState &state);
void update_modulators(ModState &state, double sr, double modRate);
void update_quantum_walk(ModState &state, double sr, double uwalkRate);
void reset_quantum_walk(ModState &state);
void apply_quantum_dither(ModState &state, std::array<double, kFdnSize> &vector, double coherence);
void apply_quantum_walk(ModState &state, std::array<double, kFdnSize> &feedback, double coherence);

} // namespace kbeyond::dsp

