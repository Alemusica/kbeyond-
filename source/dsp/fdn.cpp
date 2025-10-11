#include "fdn.h"

#include <algorithm>
#include <cmath>

namespace kbeyond::dsp {

namespace {
constexpr double kGoldenRatio = 1.6180339887498948482;
}

void setup_fdn(FdnState &state,
               double sr,
               double size,
               double modDepth,
               const std::vector<double> &latePattern) {
    const double baseMs = lerp(15.0, 55.0, size);
    const double spread = lerp(0.65, 1.35, size);
    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double idx = (i < static_cast<int>(latePattern.size()))
                               ? latePattern[static_cast<std::size_t>(i)]
                               : static_cast<double>(i) / static_cast<double>(std::max(1, static_cast<int>(kFdnSize) - 1));
        double ms = baseMs * std::pow(kGoldenRatio, (idx - 0.5) * spread);
        ms = clampd(ms, 8.0, 400.0);
        const double len = ms2samp(ms, sr);
        std::size_t bufLen = static_cast<std::size_t>(std::ceil(len + modDepth + 16.0));
        bufLen = std::max<std::size_t>(bufLen, 64);
        state.lines[static_cast<std::size_t>(i)].setup(bufLen);
        state.lengths[static_cast<std::size_t>(i)] = clampd(len, 8.0, static_cast<double>(bufLen) - 4.0);
        state.reads[static_cast<std::size_t>(i)] = state.lengths[static_cast<std::size_t>(i)];
        state.outputs[static_cast<std::size_t>(i)] = 0.0;
        state.feedback[static_cast<std::size_t>(i)] = 0.0;
    }
}

void read_lines(FdnState &state,
                double modDepth,
                std::array<double, kFdnSize> &phase,
                const std::array<double, kFdnSize> &phaseInc,
                std::array<filters::Tilt, kFdnSize> &tilts,
                std::array<filters::OnePoleLP, kFdnSize> &lp,
                std::array<double, kFdnSize> &vec) {
    for (int l = 0; l < static_cast<int>(kFdnSize); ++l) {
        double mod = 0.0;
        if (modDepth > 0.0) {
            mod = std::sin(phase[static_cast<std::size_t>(l)]) * modDepth;
            phase[static_cast<std::size_t>(l)] += phaseInc[static_cast<std::size_t>(l)];
            if (phase[static_cast<std::size_t>(l)] > 2.0 * M_PI)
                phase[static_cast<std::size_t>(l)] -= 2.0 * M_PI;
        }
        const double read = clampd(state.lengths[static_cast<std::size_t>(l)] + mod,
                                   2.0,
                                   static_cast<double>(state.lines[static_cast<std::size_t>(l)].size()) - 3.0);
        state.reads[static_cast<std::size_t>(l)] = read;
        double sig = state.lines[static_cast<std::size_t>(l)].readFrac(read);
        sig = tilts[static_cast<std::size_t>(l)].process(sig);
        sig = lp[static_cast<std::size_t>(l)].process(sig);
        vec[static_cast<std::size_t>(l)] = sig;
        state.outputs[static_cast<std::size_t>(l)] = sig;
    }
}

void write_feedback(FdnState &state,
                    double midIn,
                    double sideIn,
                    uint32_t &rng,
                    const std::array<double, kFdnSize> &inWeights,
                    const std::array<double, kFdnSize> &outWeightsL,
                    const std::array<double, kFdnSize> &outWeightsR,
                    const std::array<double, kFdnSize> &decay,
                    std::array<filters::OnePoleLP, kFdnSize> &low,
                    std::array<filters::OnePoleLP, kFdnSize> &high,
                    double dampLF,
                    double dampMF,
                    double dampHF,
                    double &tailL,
                    double &tailR) {
    tailL = 0.0;
    tailR = 0.0;
    for (int l = 0; l < static_cast<int>(kFdnSize); ++l) {
        const double fb = state.feedback[static_cast<std::size_t>(l)];
        const double lowBand = low[static_cast<std::size_t>(l)].process(fb);
        const double band = high[static_cast<std::size_t>(l)].process(fb);
        const double highBand = fb - band;
        const double midBand = band - lowBand;
        const double damped = lowBand * dampLF + midBand * dampMF + highBand * dampHF;
        const double feedback = damped * decay[static_cast<std::size_t>(l)];
        const double injection = midIn * inWeights[static_cast<std::size_t>(l)] +
                                 sideIn * outWeightsL[static_cast<std::size_t>(l)] * 0.15;
        state.lines[static_cast<std::size_t>(l)].write(feedback + injection + tiny_noise(rng));
        tailL += state.outputs[static_cast<std::size_t>(l)] * outWeightsL[static_cast<std::size_t>(l)];
        tailR += state.outputs[static_cast<std::size_t>(l)] * outWeightsR[static_cast<std::size_t>(l)];
    }
}

} // namespace kbeyond::dsp

