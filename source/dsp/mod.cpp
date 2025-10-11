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

    std::array<double, kFdnSize> neighborMix{};
    const double stateMix = lerp(0.08, 0.25, coherenceAmt);
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double phase = state.uwalkPhase[static_cast<std::size_t>(i)];
        const double sinA = std::sin(phase);
        const double sinB = std::sin(phase * 1.7320508075688772 + static_cast<double>(i) * 0.4115);
        const double blend = 0.5 * (sinA + sinB);
        const int a = (i + 1) % static_cast<int>(kFdnSize);
        const int b = (i + 5) % static_cast<int>(kFdnSize);
        const double neighbor = lerp(feedback[static_cast<std::size_t>(a)],
                                     feedback[static_cast<std::size_t>(b)],
                                     0.5 * (blend + 1.0));
        neighborMix[static_cast<std::size_t>(i)] =
            lerp(feedback[static_cast<std::size_t>(i)], neighbor, coherenceAmt);

        const double drift = neighborMix[static_cast<std::size_t>(i)] * (0.6 * sinA + 0.4 * sinB);
        state.uwalkState[static_cast<std::size_t>(i)] =
            lerp(state.uwalkState[static_cast<std::size_t>(i)], drift, stateMix);

        state.uwalkPhase[static_cast<std::size_t>(i)] += state.uwalkPhaseInc[static_cast<std::size_t>(i)];
        if (state.uwalkPhase[static_cast<std::size_t>(i)] > kTwoPi)
            state.uwalkPhase[static_cast<std::size_t>(i)] -= kTwoPi;
    }

    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double modulation = state.uwalkState[static_cast<std::size_t>(i)] * 0.2;
        feedback[static_cast<std::size_t>(i)] =
            lerp(feedback[static_cast<std::size_t>(i)],
                 neighborMix[static_cast<std::size_t>(i)] + modulation,
                 coherenceAmt);
    }
}

} // namespace kbeyond::dsp

