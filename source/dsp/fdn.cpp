#include "fdn.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <optional>
#include <vector>

#include "dsp_assert.h"
#include "prime_modes.h"

namespace kbeyond::dsp {

namespace filters = kbeyond::dsp;

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
    std::vector<long> usedLengths;
    usedLengths.reserve(static_cast<std::size_t>(kFdnSize));
    const long minSpacing = 2;
    const long minLenSamples = 8;
    const long maxLenSamples = static_cast<long>(std::floor(ms2samp(400.0, sr)));
    auto altPrime = prime_modes::generate_pattern(prime_modes::Pattern::Prime, kFdnSize * 2, 0x91A5u);
    auto altPlastica = prime_modes::generate_pattern(prime_modes::Pattern::Plastica, kFdnSize * 2, 0x91A6u);
    std::size_t altPrimeIndex = 0;
    std::size_t altPlasticaIndex = 0;

    const auto compute_ms = [&](double idxNorm) {
        double ms = baseMs * std::pow(kGoldenRatio, (idxNorm - 0.5) * spread);
        return clampd(ms, 8.0, 400.0);
    };

    for (int i = 0; i < static_cast<int>(kFdnSize); ++i) {
        const double idx = (i < static_cast<int>(latePattern.size()))
                               ? latePattern[static_cast<std::size_t>(i)]
                               : static_cast<double>(i) / static_cast<double>(std::max(1, static_cast<int>(kFdnSize) - 1));
        double ms = compute_ms(idx);
        double lenCandidate = ms2samp(ms, sr);

        const auto conflict_free = [&](long candidate) {
            for (long prev : usedLengths) {
                if (std::llabs(candidate - prev) < minSpacing)
                    return false;
                long g = std::gcd(std::llabs(candidate), std::llabs(prev));
                if (g > 1 && g <= 5)
                    return false;
            }
            return true;
        };

        const auto try_length = [&](double lengthSamples) -> std::optional<long> {
            long base = clampl(static_cast<long>(std::llround(lengthSamples)), minLenSamples, maxLenSamples);
            constexpr int offsets[] = {0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5};
            for (int offset : offsets) {
                long candidate = clampl(base + offset, minLenSamples, maxLenSamples);
                if (conflict_free(candidate))
                    return candidate;
            }
            for (long delta = minSpacing; delta < minSpacing + 96; ++delta) {
                long up = clampl(base + delta, minLenSamples, maxLenSamples);
                if (conflict_free(up))
                    return up;
                long down = clampl(base - delta, minLenSamples, maxLenSamples);
                if (conflict_free(down))
                    return down;
            }
            return std::nullopt;
        };

        auto resolved = try_length(lenCandidate);
        double chosenIdx = idx;
        double chosenLen = lenCandidate;

        if (!resolved) {
            while (altPrimeIndex < altPrime.size() && !resolved) {
                chosenIdx = altPrime[altPrimeIndex++];
                ms = compute_ms(chosenIdx);
                chosenLen = ms2samp(ms, sr);
                resolved = try_length(chosenLen);
            }
            while (altPlasticaIndex < altPlastica.size() && !resolved) {
                chosenIdx = altPlastica[altPlasticaIndex++];
                ms = compute_ms(chosenIdx);
                chosenLen = ms2samp(ms, sr);
                resolved = try_length(chosenLen);
            }
            if (!resolved) {
                long base = clampl(static_cast<long>(std::llround(chosenLen)), minLenSamples, maxLenSamples);
                for (long delta = minSpacing; delta < minSpacing + 128 && !resolved; ++delta) {
                    long up = clampl(base + delta, minLenSamples, maxLenSamples);
                    if (conflict_free(up)) {
                        resolved = up;
                        break;
                    }
                    long down = clampl(base - delta, minLenSamples, maxLenSamples);
                    if (conflict_free(down)) {
                        resolved = down;
                        break;
                    }
                }
            }
        }

        long finalSamples = resolved.value_or(clampl(static_cast<long>(std::llround(chosenLen)), minLenSamples, maxLenSamples));
        dsp_assert_msg(conflict_free(finalSamples), "FDN length spacing/gcd violation");
        usedLengths.push_back(finalSamples);

        std::size_t bufLen = static_cast<std::size_t>(std::ceil(static_cast<double>(finalSamples) + modDepth + 16.0));
        bufLen = std::max<std::size_t>(bufLen, 64);
        const std::size_t idxLine = static_cast<std::size_t>(i);
        state.lines[idxLine].setup(bufLen);
        const double finalLength = clampd(static_cast<double>(finalSamples), 8.0, static_cast<double>(bufLen) - 4.0);
        state.lengths[idxLine] = finalLength;
        state.reads[idxLine] = finalLength;
        state.outputs[idxLine] = 0.0;
        state.feedback[idxLine] = 0.0;
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

