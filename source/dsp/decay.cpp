#include "decay.h"

#include <algorithm>
#include <cmath>

namespace kbeyond::dsp {

void update_decay(DecayState &state,
                  double sr,
                  double regen,
                  double decaySeconds,
                  double dampLFNorm,
                  double dampMFNorm,
                  double dampHFNorm,
                  const std::array<double, kFdnSize> &fdnLengths) {
    const double regenNorm = clampd(regen, 0.0, 0.999);
    const double dampLFScale = clampd(dampLFNorm, 0.0, 1.0);
    const double dampMFScale = clampd(dampMFNorm, 0.0, 1.0);
    const double dampHFScale = clampd(dampHFNorm, 0.0, 1.0);

    state.dampLF = lerp(1.0, 0.35, dampLFScale);
    state.dampMF = lerp(1.0, 0.45, dampMFScale);
    state.dampHF = lerp(1.0, 0.12, dampHFScale);

    const double minDecay = 0.05;
    const double rt60 = std::max(decaySeconds, 0.0);
    const double baseTau = std::max(rt60, minDecay);

    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double delaySamples = fdnLengths[static_cast<std::size_t>(i)];
        const double delaySeconds = sr > 0.0 ? delaySamples / sr : 0.0;
        double gain = 0.0;
        if (delaySeconds > 0.0) {
            if (regenNorm <= 0.0) {
                gain = 0.0;
            } else if (rt60 > 0.0) {
                const double regenWarp = std::pow(regenNorm, 0.65);
                const double timeScale = lerp(0.45, 1.0, regenWarp);
                const double tau = baseTau * timeScale;
                const double exponent = (-3.0 * delaySeconds) / tau;
                const double timeGain = std::pow(10.0, exponent);
                const double blend = lerp(0.3, 0.85, regenWarp);
                const double mixedGain = lerp(regenNorm, timeGain, blend);
                const double floorGain = lerp(0.02, 0.72, regenWarp);
                gain = std::max(mixedGain, floorGain);
            } else {
                const double regenWarp = std::pow(regenNorm, 0.65);
                const double floorGain = lerp(0.02, 0.72, regenWarp);
                gain = std::max(regenNorm, floorGain);
            }
        }
        state.perLine[static_cast<std::size_t>(i)] = clampd(gain, 0.0, 0.99995);
    }
}

} // namespace kbeyond::dsp

