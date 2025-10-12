#include "mod.h"

#include <algorithm>
#include <cmath>

namespace kbeyond::dsp {

namespace {
constexpr double kTwoPi = 2.0 * M_PI;

double wrap_phase(double phase) {
    if (phase >= kTwoPi)
        phase -= kTwoPi;
    else if (phase < 0.0)
        phase += kTwoPi;
    return phase;
}

} // namespace

void init_fdn_phases(ModState &state) {
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        double phase = std::fmod(static_cast<double>(i + 1) * 1.2345, kTwoPi);
        if (phase < 0.0)
            phase += kTwoPi;
        state.fdnPhase[static_cast<std::size_t>(i)] = phase;
    }
}

void update_modulators(ModState &state, double sr, double modRate) {
    const double rate = clampd(modRate, 0.0, 5.0);
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double idx = static_cast<double>(i + 1) / static_cast<double>(kFdnSize);
        const double warp = lerp(0.4, 1.2, idx);
        state.fdnPhaseInc[static_cast<std::size_t>(i)] = (sr > 0.0) ? (kTwoPi * rate * warp / sr) : 0.0;
        state.fdnPhase[static_cast<std::size_t>(i)] = std::fmod(state.fdnPhase[static_cast<std::size_t>(i)], kTwoPi);
        if (state.fdnPhase[static_cast<std::size_t>(i)] < 0.0)
            state.fdnPhase[static_cast<std::size_t>(i)] += kTwoPi;
    }
}

void update_quantum_walk(ModState &state, double sr, double uwalkRate) {
    const double rate = clampd(uwalkRate, 0.0, 8.0);
    const double incBase = (rate > 0.0 && sr > 0.0) ? (kTwoPi * rate) / sr : 0.0;
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double idx = static_cast<double>(i + 1) / static_cast<double>(kFdnSize);
        const double warp = lerp(0.35, 1.75, idx);
        state.uwalkPhaseInc[static_cast<std::size_t>(i)] = incBase * warp;
        if (rate <= 0.0)
            state.uwalkPhaseInc[static_cast<std::size_t>(i)] = 0.0;
        const double ditherWarp = lerp(0.55, 1.45, idx);
        state.qditherPhaseInc[static_cast<std::size_t>(i)] = incBase * ditherWarp;
        if (rate <= 0.0)
            state.qditherPhaseInc[static_cast<std::size_t>(i)] = 0.0;
    }
}

void reset_quantum_walk(ModState &state) {
    constexpr double g = 0.6180339887498948482;
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        double seed = ((static_cast<double>(i) + 1.0) * g + 0.37) * kTwoPi;
        double phase = std::fmod(seed, kTwoPi);
        if (phase < 0.0)
            phase += kTwoPi;
        state.uwalkPhase[static_cast<std::size_t>(i)] = phase;
        state.uwalkState[static_cast<std::size_t>(i)] = 0.0;
        double ditherSeed = ((static_cast<double>(i) + 1.0) * g + 0.11) * kTwoPi;
        double ditherPhase = std::fmod(ditherSeed, kTwoPi);
        if (ditherPhase < 0.0)
            ditherPhase += kTwoPi;
        state.qditherPhase[static_cast<std::size_t>(i)] = ditherPhase;
    }
}

void apply_quantum_dither(ModState &state, std::array<double, kFdnSize> &vector, double coherence) {
    const double coherenceAmt = clampd(coherence, 0.0, 1.0);
    const double maxAngle = lerp(0.0, 0.42, coherenceAmt * coherenceAmt);
    const bool advancePhaseOnly = (maxAngle <= 0.0);

    const auto rotate_pair = [&](int idx, int jdx, double angle) {
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        const double a = vector[static_cast<std::size_t>(idx)];
        const double b = vector[static_cast<std::size_t>(jdx)];
        vector[static_cast<std::size_t>(idx)] = c * a - s * b;
        vector[static_cast<std::size_t>(jdx)] = s * a + c * b;
    };

    const auto advance_phase = [&](int idx) {
        state.qditherPhase[static_cast<std::size_t>(idx)] =
            wrap_phase(state.qditherPhase[static_cast<std::size_t>(idx)] +
                       state.qditherPhaseInc[static_cast<std::size_t>(idx)]);
    };

    if (!advancePhaseOnly) {
        for (int start = 0; start < 2; ++start) {
            for (int i = start; i < static_cast<int>(kFdnSize); i += 2) {
                const int j = (i + 1) % static_cast<int>(kFdnSize);
                const double phase = state.qditherPhase[static_cast<std::size_t>(i)];
                const double angle = maxAngle * std::sin(phase);
                rotate_pair(i, j, angle);
                advance_phase(i);
            }
        }
        for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
            const int j = (i + 5) % static_cast<int>(kFdnSize);
            const double phase = state.qditherPhase[static_cast<std::size_t>((i + 3) % static_cast<int>(kFdnSize))];
            const double angle = 0.5 * maxAngle * std::sin(phase * 0.5);
            rotate_pair(i, j, angle);
        }
    } else {
        for (int i = 0; i < static_cast<int>(kFdnSize); ++i)
            advance_phase(i);
    }
}

void apply_quantum_walk(ModState &state, std::array<double, kFdnSize> &feedback, double coherence) {
    const double coherenceAmt = clampd(coherence, 0.0, 1.0);
    if (coherenceAmt <= 0.0)
        return;

    auto rotate_pair = [](std::array<double, kFdnSize> &vec, int idx, int jdx, double angle) {
        if (idx == jdx)
            return;
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        const std::size_t i = static_cast<std::size_t>(idx);
        const std::size_t j = static_cast<std::size_t>(jdx);
        const double a = vec[i];
        const double b = vec[j];
        vec[i] = c * a - s * b;
        vec[j] = s * a + c * b;
    };

    std::array<double, kFdnSize> rotated = feedback;

    const double baseAngleScale = lerp(0.0, 0.35, coherenceAmt * coherenceAmt);
    const double crossAngleScale = baseAngleScale * 0.6;

    std::array<double, kFdnSize> primaryAngles{};
    std::array<double, kFdnSize> secondaryAngles{};
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double phase = state.uwalkPhase[static_cast<std::size_t>(i)];
        primaryAngles[static_cast<std::size_t>(i)] = baseAngleScale * std::sin(phase);
        secondaryAngles[static_cast<std::size_t>(i)] =
            crossAngleScale * std::sin(phase * 1.7320508075688772 + static_cast<double>(i) * 0.4115);
        state.uwalkPhase[static_cast<std::size_t>(i)] =
            wrap_phase(phase + state.uwalkPhaseInc[static_cast<std::size_t>(i)]);
    }

    for (int start = 0; start < 2; ++start) {
        for (int i = start; i < static_cast<int>(kFdnSize); i += 2) {
            const int neighbor = (i + 1) % static_cast<int>(kFdnSize);
            rotate_pair(rotated, i, neighbor, primaryAngles[static_cast<std::size_t>(i)]);
        }
    }

    for (int offset = 0; offset < 4; ++offset) {
        for (int i = offset; i < static_cast<int>(kFdnSize); i += 4) {
            const int neighbor = (i + 5) % static_cast<int>(kFdnSize);
            rotate_pair(rotated, i, neighbor, secondaryAngles[static_cast<std::size_t>(i)]);
        }
    }

    const double stateMix = lerp(0.04, 0.16, coherenceAmt);
    const double modulationScale = lerp(0.0, 0.18, coherenceAmt);

    std::array<double, kFdnSize> crossCache{};
    double crossMean = 0.0;
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const int partner = (i + 9) % static_cast<int>(kFdnSize);
        const double crossEnergy = rotated[static_cast<std::size_t>(i)] * rotated[static_cast<std::size_t>(partner)];
        crossCache[static_cast<std::size_t>(i)] = crossEnergy;
        crossMean += crossEnergy;
    }
    crossMean /= static_cast<double>(kFdnSize);

    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const int partner = (i + 9) % static_cast<int>(kFdnSize);
        const double centered = crossCache[static_cast<std::size_t>(i)] - crossMean;
        const double target = clampd(centered * 1.35, -1.0, 1.0);
        double &stateValue = state.uwalkState[static_cast<std::size_t>(i)];
        stateValue = lerp(stateValue, target, stateMix);
        const double modAngle = modulationScale * stateValue;
        rotate_pair(rotated, i, partner, modAngle);
    }

    feedback = rotated;
}

} // namespace kbeyond::dsp

