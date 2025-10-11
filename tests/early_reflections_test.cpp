#include "kbeyond_tilde.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

namespace {

bool test_early_reflections() {
    t_kbeyond x{};
    x.setup_sr(48000.0);
    x.width = 1.0;
    x.early = 1.0;
    x.focus = 1.0;
    x.predelay = 0.0;
    x.setup_predelay();

    constexpr int frames = 4096;
    std::vector<double> left(frames, 0.0);
    std::vector<double> right(frames, 0.0);
    auto renderImpulse = [&](t_kbeyond& inst, double leftIn, double rightIn, std::vector<double>& outL, std::vector<double>& outR) {
        for (int n = 0; n < frames; ++n) {
            double inL = (n == 0) ? leftIn : 0.0;
            double inR = (n == 0) ? rightIn : 0.0;
            double l = 0.0;
            double r = 0.0;
            inst.render_early(inL, inR, inst.width, inst.early, inst.focus, l, r);
            outL[n] = l;
            outR[n] = r;
        }
    };

    renderImpulse(x, 1.0, 1.0, left, right);

    auto findPeaks = [&](const std::vector<double>& lVec, const std::vector<double>& rVec) {
        double lp = 0.0;
        double rp = 0.0;
        for (int n = 0; n < frames; ++n) {
            lp = std::max(lp, std::abs(lVec[n]));
            rp = std::max(rp, std::abs(rVec[n]));
        }
        return std::make_pair(lp, rp);
    };

    // Ensure a left-only impulse produces a wider response on the left channel.
    std::vector<double> leftOnly(frames, 0.0);
    std::vector<double> rightOnly(frames, 0.0);
    renderImpulse(x, 1.0, 0.0, leftOnly, rightOnly);

    auto stereoPeaks = findPeaks(left, right);
    double peakL = stereoPeaks.first;
    double peakR = stereoPeaks.second;

    const double eps = 1.0e-12;
    double diffDb = 20.0 * std::log10((peakL + eps) / (peakR + eps));

    if (std::abs(diffDb) > 0.1) {
        std::cerr << "Mono impulse early response imbalance: " << diffDb << " dB" << std::endl;
        return false;
    }

    auto monoPeaks = findPeaks(leftOnly, rightOnly);
    double leftPeak = monoPeaks.first;
    double rightPeak = monoPeaks.second;

    if (!(leftPeak > rightPeak * 1.1)) {
        std::cerr << "Left-only impulse did not produce wider left response" << std::endl;
        return false;
    }

    double baselineRatio = leftPeak / (rightPeak + eps);

    t_kbeyond focusWide{};
    focusWide.setup_sr(48000.0);
    focusWide.width = 1.6;
    focusWide.early = 1.0;
    focusWide.focus = 1.0;
    focusWide.predelay = 0.0;
    focusWide.setup_predelay();

    std::vector<double> wideLeft(frames, 0.0);
    std::vector<double> wideRight(frames, 0.0);
    renderImpulse(focusWide, 1.0, 0.0, wideLeft, wideRight);
    auto widePeaks = findPeaks(wideLeft, wideRight);
    double wideRatio = widePeaks.first / (widePeaks.second + eps);

    t_kbeyond focusTight{};
    focusTight.setup_sr(48000.0);
    focusTight.width = focusWide.width;
    focusTight.early = focusWide.early;
    focusTight.focus = 0.25;
    focusTight.predelay = 0.0;
    focusTight.setup_predelay();

    std::vector<double> tightLeft(frames, 0.0);
    std::vector<double> tightRight(frames, 0.0);
    renderImpulse(focusTight, 1.0, 0.0, tightLeft, tightRight);
    auto tightPeaks = findPeaks(tightLeft, tightRight);
    double tightRatio = tightPeaks.first / (tightPeaks.second + eps);

    if (!(wideRatio > baselineRatio)) {
        std::cerr << "Wide focus ratio did not exceed baseline ratio" << std::endl;
        return false;
    }
    if (!(tightRatio < wideRatio)) {
        std::cerr << "Tight focus ratio not reduced" << std::endl;
        return false;
    }
    if (tightRatio > 1.05) {
        std::cerr << "Tight focus remained too wide: ratio=" << tightRatio << std::endl;
        return false;
    }

    return true;
}

double energy(const std::array<double, t_kbeyond::N> &v) {
    double sum = 0.0;
    for (double s : v)
        sum += s * s;
    return sum;
}

bool test_mix_modes() {
    t_kbeyond x{};
    x.setup_sr(48000.0);

    std::array<double, t_kbeyond::N> input{};
    input[0] = 1.0;

    std::array<double, t_kbeyond::N> outHouse{};
    std::array<double, t_kbeyond::N> outWht{};
    std::array<double, t_kbeyond::N> outHybrid{};

    x.mixMode = t_kbeyond::MixMode::Householder;
    x.update_diffusion();
    x.apply_diffusion(input, outHouse);

    x.mixMode = t_kbeyond::MixMode::WHT;
    x.update_diffusion();
    x.apply_diffusion(input, outWht);

    x.mixMode = t_kbeyond::MixMode::Hybrid;
    x.update_diffusion();
    x.apply_diffusion(input, outHybrid);

    double eHouse = energy(outHouse);
    double eWht = energy(outWht);
    double eHybrid = energy(outHybrid);

    const double tol = 1.0e-9;
    if (std::abs(eHouse - 1.0) > tol) {
        std::cerr << "Householder energy mismatch: " << eHouse << std::endl;
        return false;
    }
    if (std::abs(eWht - 1.0) > tol) {
        std::cerr << "WHT energy mismatch: " << eWht << std::endl;
        return false;
    }
    if (std::abs(eHybrid - 1.0) > tol) {
        std::cerr << "Hybrid energy mismatch: " << eHybrid << std::endl;
        return false;
    }

    double diffHW = 0.0;
    for (std::size_t i = 0; i < outHouse.size(); ++i)
        diffHW += std::abs(outHouse[i] - outWht[i]);
    if (diffHW < 1.0e-3) {
        std::cerr << "WHT and Householder responses are unexpectedly similar" << std::endl;
        return false;
    }

    double diffHouseHybrid = 0.0;
    double diffWhtHybrid = 0.0;
    for (std::size_t i = 0; i < outHybrid.size(); ++i) {
        diffHouseHybrid += std::abs(outHouse[i] - outHybrid[i]);
        diffWhtHybrid += std::abs(outWht[i] - outHybrid[i]);
    }
    if (diffHouseHybrid < 1.0e-3 || diffWhtHybrid < 1.0e-3) {
        std::cerr << "Hybrid response too close to pure modes" << std::endl;
        return false;
    }

    return true;
}

bool test_quantum_dither_energy() {
    t_kbeyond x{};
    x.setup_sr(48000.0);

    std::array<double, t_kbeyond::N> state{};
    for (int i = 0; i < t_kbeyond::N; ++i)
        state[i] = (i % 3 == 0) ? 0.75 : ((i % 3 == 1) ? -0.5 : 0.25);

    auto stateCopy = state;

    x.coherence = 0.0;
    x.update_quantum_walk();
    x.apply_quantum_dither(stateCopy);
    if (stateCopy != state) {
        std::cerr << "Quantum dither modified state when coherence was zero" << std::endl;
        return false;
    }

    x.coherence = 0.9;
    x.uwalkRate = 0.0;
    x.update_quantum_walk();
    stateCopy = state;
    x.apply_quantum_dither(stateCopy);
    double eBase = energy(state);
    double eAfter = energy(stateCopy);
    if (std::abs(eAfter - eBase) > 1.0e-9) {
        std::cerr << "Quantum dither energy mismatch at zero rate" << std::endl;
        return false;
    }

    x.uwalkRate = 0.75;
    x.update_quantum_walk();
    double targetEnergy = energy(state);
    for (int n = 0; n < 256; ++n) {
        x.apply_quantum_dither(state);
        double e = energy(state);
        if (std::abs(e - targetEnergy) > 1.0e-9) {
            std::cerr << "Quantum dither energy drift at iteration " << n << std::endl;
            return false;
        }
    }

    return true;
}

struct DecayResponse {
    double minCoeff = 0.0;
    double maxCoeff = 0.0;
    double energy = 0.0;
};

DecayResponse measure_decay_response(double regenValue, double decayValue) {
    t_kbeyond x{};
    x.setup_sr(48000.0);
    x.mix = 1.0;
    x.early = 0.0;
    x.predelay = 0.0;
    x.setup_predelay();
    x.regen = regenValue;
    x.decay = decayValue;
    x.update_decay();

    DecayResponse response{};
    response.minCoeff = *std::min_element(x.decayState.perLine.begin(), x.decayState.perLine.end());
    response.maxCoeff = *std::max_element(x.decayState.perLine.begin(), x.decayState.perLine.end());

    constexpr int frames = 48000;
    double widthNorm = clampd(x.width, 0.0, 2.0);
    double earlyAmt = clampd(x.early, 0.0, 1.0);
    double focusAmt = clampd(x.focus, 0.0, 1.0);
    double moddepth = clampd(x.moddepth, 0.0, 32.0);

    for (int n = 0; n < frames; ++n) {
        double inL = (n == 0) ? 1.0 : 0.0;
        double inR = inL;

        double predLen = x.predSamps;
        double predOutL = x.predL.readFrac(predLen);
        double predOutR = x.predR.readFrac(predLen);
        x.predL.write(inL + kbeyond::dsp::tiny_noise(x.rng));
        x.predR.write(inR + kbeyond::dsp::tiny_noise(x.rng));

        double midIn = 0.5 * (predOutL + predOutR);
        double sideIn = 0.5 * (predOutL - predOutR);

        double earlyL = 0.0;
        double earlyR = 0.0;
        double clusterAmt = clampd(x.laser, 0.0, 1.0);
        x.earlySection.render(predOutL,
                              predOutR,
                              widthNorm,
                              earlyAmt,
                              focusAmt,
                              clusterAmt,
                              earlyL,
                              earlyR,
                              x.rng);

        std::array<double, t_kbeyond::N> vec{};
        kbeyond::dsp::read_lines(x.fdnState,
                                 moddepth,
                                 x.modState.fdnPhase,
                                 x.modState.fdnPhaseInc,
                                 x.fdn_tilt,
                                 x.fdn_lp,
                                 vec);

        x.apply_diffusion(vec, x.fdnState.feedback);
        x.apply_quantum_walk(x.fdnState.feedback);

        double tailL = 0.0;
        double tailR = 0.0;
        kbeyond::dsp::write_feedback(x.fdnState,
                                     midIn,
                                     sideIn,
                                     x.rng,
                                     x.inWeights,
                                     x.outWeightsL,
                                     x.outWeightsR,
                                     x.decayState.perLine,
                                     x.fdn_low,
                                     x.fdn_high,
                                     x.decayState.dampLF,
                                     x.decayState.dampMF,
                                     x.decayState.dampHF,
                                     tailL,
                                     tailR);

        double wetL = earlyL + tailL;
        double wetR = earlyR + tailR;
        response.energy += wetL * wetL + wetR * wetR;
    }

    return response;
}

bool test_decay_regen_response() {
    std::array<double, 3> values{0.1, 0.5, 0.9};
    std::array<DecayResponse, 3> responses{};

    for (std::size_t i = 0; i < values.size(); ++i)
        responses[i] = measure_decay_response(values[i], values[i]);

    for (const auto &resp : responses) {
        if (resp.minCoeff < 0.015) {
            std::cerr << "FDN decay floor too low: " << resp.minCoeff << std::endl;
            return false;
        }
    }

    if (!(responses[1].energy > responses[0].energy * 1.03)) {
        std::cerr << "Energy at regen/decay 0.5 did not exceed 0.1 setting" << std::endl;
        return false;
    }
    if (!(responses[2].energy > responses[1].energy * 1.06)) {
        std::cerr << "Energy at regen/decay 0.9 did not exceed 0.5 setting" << std::endl;
        return false;
    }

    if (!(responses[2].maxCoeff - responses[0].maxCoeff > 0.05)) {
        std::cerr << "FDN decay coefficients lack sufficient spread" << std::endl;
        return false;
    }

    return true;
}

} // namespace

int main() {
    if (!test_early_reflections())
        return 1;
    if (!test_mix_modes())
        return 2;
    if (!test_quantum_dither_energy())
        return 3;
    if (!test_decay_regen_response())
        return 4;
    return 0;
}
