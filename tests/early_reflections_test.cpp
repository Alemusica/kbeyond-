#include "kbeyond_tilde.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

using kbeyond::dsp::clampd;

namespace motion_tests {
bool run_motion_width_response();
bool run_motion_moddepth_response();
bool run_mid_side_mix_normalization();
}

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

    t_kbeyond spreadWide{};
    spreadWide.setup_sr(48000.0);
    spreadWide.width = 1.6;
    spreadWide.early = 1.0;
    spreadWide.focus = 1.0;
    spreadWide.predelay = 0.0;
    spreadWide.setup_predelay();

    std::vector<double> wideLeft(frames, 0.0);
    std::vector<double> wideRight(frames, 0.0);
    renderImpulse(spreadWide, 1.0, 0.0, wideLeft, wideRight);
    auto widePeaks = findPeaks(wideLeft, wideRight);
    double wideRatio = widePeaks.first / (widePeaks.second + eps);

    t_kbeyond spreadTight{};
    spreadTight.setup_sr(48000.0);
    spreadTight.width = spreadWide.width;
    spreadTight.early = 0.2;
    spreadTight.focus = spreadWide.focus;
    spreadTight.predelay = 0.0;
    spreadTight.setup_predelay();

    std::vector<double> tightLeft(frames, 0.0);
    std::vector<double> tightRight(frames, 0.0);
    renderImpulse(spreadTight, 1.0, 0.0, tightLeft, tightRight);
    auto tightPeaks = findPeaks(tightLeft, tightRight);
    double tightRatio = tightPeaks.first / (tightPeaks.second + eps);

    if (!(wideRatio > baselineRatio)) {
        std::cerr << "Wide spread ratio did not exceed baseline ratio" << std::endl;
        return false;
    }
    if (!(tightRatio < wideRatio)) {
        std::cerr << "Tight spread ratio not reduced" << std::endl;
        return false;
    }
    if (tightRatio > 1.05) {
        std::cerr << "Tight spread remained too wide: ratio=" << tightRatio << std::endl;
        return false;
    }

    auto energySum = [&](const std::vector<double> &lVec, const std::vector<double> &rVec) {
        double acc = 0.0;
        for (int n = 0; n < frames; ++n)
            acc += lVec[n] * lVec[n] + rVec[n] * rVec[n];
        return acc;
    };

    t_kbeyond focusFull{};
    focusFull.setup_sr(48000.0);
    focusFull.width = spreadWide.width;
    focusFull.early = spreadWide.early;
    focusFull.focus = 1.0;
    focusFull.predelay = 0.0;
    focusFull.setup_predelay();

    t_kbeyond focusNarrow{};
    focusNarrow.setup_sr(48000.0);
    focusNarrow.width = focusFull.width;
    focusNarrow.early = focusFull.early;
    focusNarrow.focus = 0.25;
    focusNarrow.predelay = 0.0;
    focusNarrow.setup_predelay();

    std::vector<double> focusFullLeft(frames, 0.0);
    std::vector<double> focusFullRight(frames, 0.0);
    std::vector<double> focusSoftLeft(frames, 0.0);
    std::vector<double> focusSoftRight(frames, 0.0);
    renderImpulse(focusFull, 1.0, 0.0, focusFullLeft, focusFullRight);
    renderImpulse(focusNarrow, 1.0, 0.0, focusSoftLeft, focusSoftRight);

    double fullEnergy = energySum(focusFullLeft, focusFullRight);
    double reducedEnergy = energySum(focusSoftLeft, focusSoftRight);

    if (!(reducedEnergy < fullEnergy * 0.35)) {
        std::cerr << "Focus control did not sufficiently reduce early energy" << std::endl;
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

bool test_prime_mode_spacing() {
    struct PresetCase {
        const char* name;
        double size;
        double focus;
        double width;
        double early;
        double moddepth;
        double laser;
        double laserFocus;
        double laserGate;
        double laserWindow;
        double laserDiffusion;
    };

    const std::vector<PresetCase> presets = {
        {"Factory Default", kbeyond::dsp::AttributeDefaults::size, kbeyond::dsp::AttributeDefaults::focus,
         kbeyond::dsp::AttributeDefaults::width, kbeyond::dsp::AttributeDefaults::early,
         kbeyond::dsp::AttributeDefaults::moddepth, kbeyond::dsp::AttributeDefaults::laser,
         kbeyond::dsp::AttributeDefaults::laserFocus, kbeyond::dsp::AttributeDefaults::laserGate,
         kbeyond::dsp::AttributeDefaults::laserWindow, kbeyond::dsp::AttributeDefaults::laserDiffusion},
        {"Hall Larga", 0.92, 1.0, 1.6, 0.35, 6.0, 0.0, kbeyond::dsp::AttributeDefaults::laserFocus,
         kbeyond::dsp::AttributeDefaults::laserGate, kbeyond::dsp::AttributeDefaults::laserWindow,
         kbeyond::dsp::AttributeDefaults::laserDiffusion},
        {"Plate Densa", 0.58, 0.95, 1.3, 0.28, 4.0, 0.0, kbeyond::dsp::AttributeDefaults::laserFocus,
         kbeyond::dsp::AttributeDefaults::laserGate, kbeyond::dsp::AttributeDefaults::laserWindow,
         kbeyond::dsp::AttributeDefaults::laserDiffusion},
        {"Room Intima", 0.35, 0.6, 0.9, 0.45, 2.5, 0.0, kbeyond::dsp::AttributeDefaults::laserFocus,
         kbeyond::dsp::AttributeDefaults::laserGate, kbeyond::dsp::AttributeDefaults::laserWindow,
         kbeyond::dsp::AttributeDefaults::laserDiffusion},
        {"Laser Sweep", 0.7, 0.85, 1.45, 0.38, 5.5, 0.85, 0.65, 0.32, 0.42, 0.72},
    };

    for (const auto& preset : presets) {
        t_kbeyond inst{};
        inst.setup_sr(48000.0);
        inst.size = preset.size;
        inst.focus = preset.focus;
        inst.width = preset.width;
        inst.early = preset.early;
        inst.moddepth = preset.moddepth;
        inst.laser = preset.laser;
        inst.laserFocus = preset.laserFocus;
        inst.laserGate = preset.laserGate;
        inst.laserWindow = preset.laserWindow;
        inst.laserDiffusion = preset.laserDiffusion;
        inst.setup_early();
        inst.setup_fdn();

        const auto& delays = inst.earlySection.debugTapDelays();
        std::vector<long> uniqueDelays;
        uniqueDelays.reserve(delays.size());
        for (long d : delays) {
            if (d > 0)
                uniqueDelays.push_back(d);
        }
        std::sort(uniqueDelays.begin(), uniqueDelays.end());
        uniqueDelays.erase(std::unique(uniqueDelays.begin(), uniqueDelays.end()), uniqueDelays.end());
        for (std::size_t i = 1; i < uniqueDelays.size(); ++i) {
            if (uniqueDelays[i] - uniqueDelays[i - 1] < 2) {
                std::cerr << "Early tap spacing violation in preset " << preset.name << std::endl;
                return false;
            }
        }

        const auto& gains = inst.earlySection.debugTapGains();
        const auto& cosines = inst.earlySection.debugTapCos();
        const auto& sines = inst.earlySection.debugTapSin();
        double energy = 0.0;
        for (int t = 0; t < t_kbeyond::kEarlyTaps; ++t) {
            double wL = gains[static_cast<std::size_t>(t)] * cosines[static_cast<std::size_t>(t)];
            double wR = gains[static_cast<std::size_t>(t)] * sines[static_cast<std::size_t>(t)];
            energy += wL * wL + wR * wR;
        }
        if (std::abs(energy - 1.0) > 1.0e-6) {
            std::cerr << "Early gain normalization failed for preset " << preset.name << std::endl;
            return false;
        }

        const auto& lengths = inst.fdnState.lengths;
        for (std::size_t i = 0; i < lengths.size(); ++i) {
            long a = static_cast<long>(std::llround(lengths[i]));
            for (std::size_t j = i + 1; j < lengths.size(); ++j) {
                long b = static_cast<long>(std::llround(lengths[j]));
                if (std::llabs(a - b) < 2) {
                    std::cerr << "FDN spacing violation between lines " << i << " and " << j
                              << " in preset " << preset.name << std::endl;
                    return false;
                }
                long g = std::gcd(std::llabs(a), std::llabs(b));
                if (g > 1 && g <= 5) {
                    std::cerr << "FDN gcd violation (" << g << ") between lines " << i << " and " << j
                              << " in preset " << preset.name << std::endl;
                    return false;
                }
            }
        }
    }

    return true;
}

struct DecayResponse {
    double minCoeff = 0.0;
    double maxCoeff = 0.0;
    double energy = 0.0;
    double tailEnergy = 0.0;
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

        double tailMid = 0.0;
        double tailSide = 0.0;
        kbeyond::dsp::write_feedback(x.fdnState,
                                     midIn,
                                     sideIn,
                                     x.rng,
                                     x.inWeights,
                                     x.outWeightsSide,
                                     x.outMidBasis,
                                     x.decayState.perLine,
                                     x.fdn_low,
                                     x.fdn_high,
                                     x.decayState.dampLF,
                                     x.decayState.dampMF,
                                     x.decayState.dampHF,
                                     tailMid,
                                     tailSide);

        double tailL = 0.0;
        double tailR = 0.0;
        x.mix_mid_side_to_lr(tailMid, tailSide, widthNorm, tailL, tailR);

        double wetL = earlyL + tailL;
        double wetR = earlyR + tailR;
        response.energy += wetL * wetL + wetR * wetR;
        response.tailEnergy += tailL * tailL + tailR * tailR;
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

    if (!(responses[1].tailEnergy > responses[0].tailEnergy * 1.03)) {
        std::cerr << "Tail energy at regen/decay 0.5 did not exceed 0.1 setting: low="
                  << responses[0].tailEnergy << " mid=" << responses[1].tailEnergy
                  << " high=" << responses[2].tailEnergy << std::endl;
        return false;
    }
    if (!(responses[2].tailEnergy > responses[1].tailEnergy * 1.06)) {
        std::cerr << "Tail energy at regen/decay 0.9 did not exceed 0.5 setting" << std::endl;
        return false;
    }

    if (!(responses[2].maxCoeff - responses[0].maxCoeff > 0.05)) {
        std::cerr << "FDN decay coefficients lack sufficient spread" << std::endl;
        return false;
    }

    return true;
}

bool test_side_impulse_width_balance() {
    t_kbeyond inst{};
    inst.setup_sr(48000.0);
    inst.mix = 1.0;
    inst.early = 0.0;
    inst.focus = 0.5;
    inst.predelay = 0.0;
    inst.setup_predelay();
    inst.regen = 0.7;
    inst.decay = 4.0;
    inst.update_decay();
    inst.width = 1.0;
    inst.apply_width(inst.width);

    constexpr int impulseFrames = 4096;
    constexpr int tailFrames = 4096;
    constexpr int totalFrames = impulseFrames + tailFrames;
    const double widthRate = 2.0 * M_PI / static_cast<double>(impulseFrames);
    const double baseWidth = 1.0;
    const double sweep = 0.85;

    double widthPhase = 0.0;
    double energyL = 0.0;
    double energyR = 0.0;
    double energyMid = 0.0;
    double energySide = 0.0;

    double dampLFValue = inst.decayState.dampLF;
    double dampMFValue = inst.decayState.dampMF;
    double dampHFValue = inst.decayState.dampHF;
    double moddepth = clampd(inst.moddepth, 0.0, 32.0);

    for (int n = 0; n < totalFrames; ++n) {
        double widthVal = clampd(baseWidth + sweep * std::sin(widthPhase), 0.0, 2.0);
        widthPhase += widthRate;
        inst.width = widthVal;
        inst.apply_width(widthVal);

        double inL = (n == 0) ? 1.0 : 0.0;
        double inR = (n == 0) ? -1.0 : 0.0;

        double predLen = inst.predSamps;
        double predOutL = inst.predL.readFrac(predLen);
        double predOutR = inst.predR.readFrac(predLen);
        inst.predL.write(inL + kbeyond::dsp::tiny_noise(inst.rng));
        inst.predR.write(inR + kbeyond::dsp::tiny_noise(inst.rng));

        double midIn = 0.5 * (predOutL + predOutR);
        double sideIn = 0.5 * (predOutL - predOutR);

        double earlyL = 0.0;
        double earlyR = 0.0;
        double clusterAmt = clampd(inst.laser, 0.0, 1.0);
        inst.earlySection.render(predOutL,
                                 predOutR,
                                 widthVal,
                                 0.0,
                                 clampd(inst.focus, 0.0, 1.0),
                                 clusterAmt,
                                 earlyL,
                                 earlyR,
                                 inst.rng);

        std::array<double, t_kbeyond::N> vec{};
        kbeyond::dsp::read_lines(inst.fdnState,
                                 moddepth,
                                 inst.modState.fdnPhase,
                                 inst.modState.fdnPhaseInc,
                                 inst.fdn_tilt,
                                 inst.fdn_lp,
                                 vec);

        inst.apply_diffusion(vec, inst.fdnState.feedback);
        inst.apply_quantum_walk(inst.fdnState.feedback);

        double tailMid = 0.0;
        double tailSide = 0.0;
        kbeyond::dsp::write_feedback(inst.fdnState,
                                     midIn,
                                     sideIn,
                                     inst.rng,
                                     inst.inWeights,
                                     inst.outWeightsSide,
                                     inst.outMidBasis,
                                     inst.decayState.perLine,
                                     inst.fdn_low,
                                     inst.fdn_high,
                                     dampLFValue,
                                     dampMFValue,
                                     dampHFValue,
                                     tailMid,
                                     tailSide);
        double tailL = 0.0;
        double tailR = 0.0;
        inst.mix_mid_side_to_lr(tailMid, tailSide, widthVal, tailL, tailR);
        energyMid += tailMid * tailMid;
        energySide += tailSide * tailSide;

        double wetL = earlyL + tailL;
        double wetR = earlyR + tailR;
        energyL += wetL * wetL;
        energyR += wetR * wetR;
    }

    const double eps = 1.0e-12;
    double totalEnergy = std::max(energyL + energyR, eps);
    double balance = std::abs(energyL - energyR) / totalEnergy;
    if (balance > 0.03) {
        std::cerr << "Side impulse energy imbalance under width modulation: balance="
                  << balance << std::endl;
        double midSideDot = 0.0;
        double midNorm = 0.0;
        double sideNorm = 0.0;
        for (int i = 0; i < t_kbeyond::N; ++i) {
            double mid = inst.outMidBasis[i];
            double side = inst.outWeightsSide[i];
            midSideDot += mid * side;
            midNorm += mid * mid;
            sideNorm += side * side;
        }
        double lrDot = 0.0;
        double diffNorm = 0.0;
        for (int i = 0; i < t_kbeyond::N; ++i) {
            double diff = inst.outWeightsL[i] - inst.outWeightsR[i];
            diffNorm += diff * diff;
            lrDot += inst.outWeightsL[i] * inst.outWeightsR[i];
        }
        std::cerr << "midSideDot=" << midSideDot << " midNorm=" << midNorm
                  << " sideNorm=" << sideNorm << " lrDot=" << lrDot
                  << " diffNorm=" << diffNorm << std::endl;
        std::cerr << "energyMid=" << energyMid << " energySide=" << energySide << std::endl;
        return false;
    }

    return true;
}

bool run_side_width_energy_test() {
    t_kbeyond inst{};
    inst.setup_sr(48000.0);
    inst.mix = 1.0;
    inst.motion = 0.0;
    inst.moddepth = 0.0;
    inst.laser = 0.0;
    inst.early = 0.0;
    inst.focus = 0.5;
    inst.predelay = 0.0;
    inst.setup_predelay();
    inst.regen = 0.7;
    inst.decay = 4.0;
    inst.update_decay();
    inst.width = 0.0;
    inst.apply_width(inst.width);
    inst.motionDetector.reset();

    constexpr long blockSize = 32;
    constexpr int settleBlocks = 128;
    constexpr int widthSteps = 41;
    constexpr double widthStep = 2.0 / static_cast<double>(widthSteps - 1);

    std::array<double, blockSize> inL{};
    std::array<double, blockSize> inR{};
    std::array<double, blockSize> outL{};
    std::array<double, blockSize> outR{};
    double* ins[2] = {inL.data(), inR.data()};
    double* outs[2] = {outL.data(), outR.data()};

    double energyL = 0.0;
    double energyR = 0.0;

    auto accumulateBlock = [&]() {
        for (long n = 0; n < blockSize; ++n) {
            double l = outL[static_cast<std::size_t>(n)];
            double r = outR[static_cast<std::size_t>(n)];
            energyL += l * l;
            energyR += r * r;
        }
    };

    for (int step = 0; step < widthSteps; ++step) {
        double widthTarget = widthStep * static_cast<double>(step);
        inst.width = widthTarget;
        inst.apply_width(widthTarget);

        for (int b = 0; b < settleBlocks; ++b) {
            std::fill(inL.begin(), inL.end(), 0.0);
            std::fill(inR.begin(), inR.end(), 0.0);
            kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);
            accumulateBlock();
        }

        std::fill(inL.begin(), inL.end(), 0.0);
        std::fill(inR.begin(), inR.end(), 0.0);
        inL[0] = 1.0;
        inR[0] = -1.0;
        kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);
        accumulateBlock();
    }

    for (int b = 0; b < settleBlocks; ++b) {
        std::fill(inL.begin(), inL.end(), 0.0);
        std::fill(inR.begin(), inR.end(), 0.0);
        kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);
        accumulateBlock();
    }

    const double eps = 1.0e-12;
    double totalEnergy = std::max(energyL + energyR, eps);
    double imbalance = std::abs(energyL - energyR) / totalEnergy;
    if (imbalance > 0.03) {
        std::cerr << "Side width sweep energy imbalance: imbalance=" << imbalance
                  << " energyL=" << energyL << " energyR=" << energyR << std::endl;
        return false;
    }

    return true;
}

bool test_wet_tail_makeup_balance() {
    constexpr long blockSize = 256;
    constexpr int tailBlocks = 6;
    constexpr double impulseAmp = 1.0e-3;
    const std::array<double, 4> mixes{0.25, 0.5, 0.75, 1.0};
    std::array<double, mixes.size()> tailRms{};

    for (std::size_t idx = 0; idx < mixes.size(); ++idx) {
        t_kbeyond inst{};
        inst.setup_sr(48000.0);
        inst.mix = mixes[idx];
        inst.early = 0.0;
        inst.predelay = 0.0;
        inst.setup_predelay();
        inst.decay = 4.0;
        inst.update_decay();

        std::array<double, blockSize> inL{};
        std::array<double, blockSize> inR{};
        std::array<double, blockSize> outL{};
        std::array<double, blockSize> outR{};
        double *ins[2] = {inL.data(), inR.data()};
        double *outs[2] = {outL.data(), outR.data()};

        for (long n = 0; n < blockSize; ++n) {
            inL[n] = (n == 0) ? impulseAmp : 0.0;
            inR[n] = inL[n];
        }
        kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);

        std::array<double, blockSize> tailL{};
        std::array<double, blockSize> tailR{};
        for (int b = 0; b < tailBlocks; ++b) {
            std::fill(inL.begin(), inL.end(), 0.0);
            std::fill(inR.begin(), inR.end(), 0.0);
            kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);
            if (b == tailBlocks - 1) {
                tailL = outL;
                tailR = outR;
            }
        }

        double energy = 0.0;
        for (long n = 0; n < blockSize; ++n) {
            double l = tailL[n];
            double r = tailR[n];
            energy += l * l + r * r;
        }
        double mean = energy / (2.0 * static_cast<double>(blockSize));
        tailRms[idx] = std::sqrt(std::max(mean, 0.0));
    }

    constexpr double eps = 1.0e-12;
    for (std::size_t idx = 1; idx < mixes.size(); ++idx) {
        double prev = tailRms[idx - 1];
        double current = tailRms[idx];
        if (current + eps < prev * 0.6) {
            std::cerr << "Wet tail RMS collapsed between mix=" << mixes[idx - 1]
                      << " (" << prev << ") and mix=" << mixes[idx]
                      << " (" << current << ")" << std::endl;
            return false;
        }
    }

    double fullWet = tailRms.back();
    double refWet = tailRms[mixes.size() - 2];
    double diffDb = 20.0 * std::log10((fullWet + eps) / (refWet + eps));
    if (diffDb < -1.0) {
        std::cerr << "Full wet tail dropped by " << diffDb
                  << " dB versus mix=" << mixes[mixes.size() - 2] << std::endl;
        return false;
    }

    // Regression: ensure that rebuilding the FDN (@size change) keeps the wet
    // compensation bounded and the tail energy consistent with a fresh instance.
    {
        constexpr double initialSize = 0.25;
        constexpr double targetSize = 0.8;
        constexpr int preChangeBlocks = 2;
        constexpr int postSettleBlocks = tailBlocks;

        auto computeRms = [&](const std::array<double, blockSize> &lVec, const std::array<double, blockSize> &rVec) {
            double energy = 0.0;
            for (long n = 0; n < blockSize; ++n)
                energy += lVec[n] * lVec[n] + rVec[n] * rVec[n];
            double mean = energy / (2.0 * static_cast<double>(blockSize));
            return std::sqrt(std::max(mean, 0.0));
        };

        t_kbeyond inst{};
        inst.setup_sr(48000.0);
        inst.mix = 1.0;
        inst.early = 0.0;
        inst.predelay = 0.0;
        inst.setup_predelay();
        inst.decay = 4.0;
        inst.update_decay();
        inst.size = initialSize;
        inst.setup_fdn();
        inst.setup_early();
        inst.update_modulators();

        std::array<double, blockSize> inL{};
        std::array<double, blockSize> inR{};
        std::array<double, blockSize> outL{};
        std::array<double, blockSize> outR{};
        double *ins[2] = {inL.data(), inR.data()};
        double *outs[2] = {outL.data(), outR.data()};

        for (long n = 0; n < blockSize; ++n) {
            inL[n] = (n == 0) ? impulseAmp : 0.0;
            inR[n] = inL[n];
        }
        kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);

        for (int b = 0; b < preChangeBlocks; ++b) {
            std::fill(inL.begin(), inL.end(), 0.0);
            std::fill(inR.begin(), inR.end(), 0.0);
            kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);
        }

        inst.size = targetSize;
        inst.setup_fdn();
        inst.setup_early();
        inst.update_modulators();
        if (inst.wetMakeup > 1.001) {
            std::cerr << "Wet makeup did not reset after @size rebuild: " << inst.wetMakeup << std::endl;
            return false;
        }

        std::fill(inL.begin(), inL.end(), 0.0);
        std::fill(inR.begin(), inR.end(), 0.0);
        kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);

        for (long n = 0; n < blockSize; ++n) {
            inL[n] = (n == 0) ? impulseAmp : 0.0;
            inR[n] = inL[n];
        }
        kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);

        std::array<double, blockSize> postTailL{};
        std::array<double, blockSize> postTailR{};
        for (int b = 0; b < postSettleBlocks; ++b) {
            std::fill(inL.begin(), inL.end(), 0.0);
            std::fill(inR.begin(), inR.end(), 0.0);
            kbeyond_perform64(&inst, nullptr, ins, 2, outs, 2, blockSize, 0, nullptr);
            if (b == postSettleBlocks - 1) {
                postTailL = outL;
                postTailR = outR;
            }
        }

        double postRms = computeRms(postTailL, postTailR);
        if (inst.wetMakeup > 2.0) {
            std::cerr << "Wet makeup exploded after @size change: " << inst.wetMakeup << std::endl;
            return false;
        }

        t_kbeyond reference{};
        reference.setup_sr(48000.0);
        reference.mix = 1.0;
        reference.early = 0.0;
        reference.predelay = 0.0;
        reference.setup_predelay();
        reference.decay = 4.0;
        reference.update_decay();
        reference.size = targetSize;
        reference.setup_fdn();
        reference.setup_early();
        reference.update_modulators();

        std::array<double, blockSize> refInL{};
        std::array<double, blockSize> refInR{};
        std::array<double, blockSize> refOutL{};
        std::array<double, blockSize> refOutR{};
        double *refIns[2] = {refInL.data(), refInR.data()};
        double *refOuts[2] = {refOutL.data(), refOutR.data()};

        for (long n = 0; n < blockSize; ++n) {
            refInL[n] = (n == 0) ? impulseAmp : 0.0;
            refInR[n] = refInL[n];
        }
        kbeyond_perform64(&reference, nullptr, refIns, 2, refOuts, 2, blockSize, 0, nullptr);

        std::array<double, blockSize> refTailL{};
        std::array<double, blockSize> refTailR{};
        for (int b = 0; b < postSettleBlocks; ++b) {
            std::fill(refInL.begin(), refInL.end(), 0.0);
            std::fill(refInR.begin(), refInR.end(), 0.0);
            kbeyond_perform64(&reference, nullptr, refIns, 2, refOuts, 2, blockSize, 0, nullptr);
            if (b == postSettleBlocks - 1) {
                refTailL = refOutL;
                refTailR = refOutR;
            }
        }

        double refRms = computeRms(refTailL, refTailR);
        double rebuildDiffDb = 20.0 * std::log10((postRms + eps) / (refRms + eps));
        if (std::abs(rebuildDiffDb) > 1.0) {
            std::cerr << "Post-rebuild wet RMS diverged by " << rebuildDiffDb << " dB" << std::endl;
            return false;
        }
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
    if (!test_prime_mode_spacing())
        return 4;
    if (!test_decay_regen_response())
        return 5;
    if (!motion_tests::run_motion_width_response())
        return 6;
    if (!test_side_impulse_width_balance())
        return 7;
    if (!run_side_width_energy_test())
        return 8;
    if (!test_wet_tail_makeup_balance())
        return 9;
    if (!motion_tests::run_motion_moddepth_response())
        return 10;
    if (!motion_tests::run_mid_side_mix_normalization())
        return 11;
    return 0;
}
