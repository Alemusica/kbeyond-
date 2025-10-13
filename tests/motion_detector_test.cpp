#include "kbeyond_tilde.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

namespace motion_tests {
namespace {

double weights_energy(const std::array<double, t_kbeyond::N> &weights) {
    double sum = 0.0;
    for (double w : weights)
        sum += w * w;
    return sum;
}

} // namespace

bool run_motion_width_response() {
    t_kbeyond x{};
    x.setup_sr(48000.0);
    x.motion = 1.0;
    x.width = 1.0;
    x.mix = 0.0;
    x.early = 0.0;
    x.predelay = 0.0;
    x.setup_predelay();
    x.regen = 0.0;
    x.decay = 0.0;
    x.update_decay();
    x.motionDetector.reset();
    x.rangeEnv = x.motionDetector.range();
    x.dopplerEnv = x.motionDetector.doppler();
    x.spreadEnv = x.motionDetector.spread();

    constexpr long frames = 64;
    std::array<double, frames> inL{};
    std::array<double, frames> inR{};
    std::array<double, frames> outL{};
    std::array<double, frames> outR{};
    double *ins[2] = {inL.data(), inR.data()};
    double *outs[2] = {outL.data(), outR.data()};

    auto runBlock = [&](double amp) {
        for (long n = 0; n < frames; ++n) {
            inL[n] = amp;
            inR[n] = amp;
        }
        kbeyond_perform64(&x, nullptr, ins, 2, outs, 2, frames, 0, nullptr);
        return x.debugWidthTarget;
    };

    runBlock(0.0);

    const int settleBlocks = 128;
    double lowWidth = 0.0;
    for (int n = 0; n < settleBlocks; ++n)
        lowWidth = runBlock(0.05);

    double highWidth = 0.0;
    for (int n = 0; n < settleBlocks; ++n)
        highWidth = runBlock(1.0);

    double energyL = weights_energy(x.outWeightsL);
    double energyR = weights_energy(x.outWeightsR);

    double widthBack = 0.0;
    for (int n = 0; n < settleBlocks; ++n)
        widthBack = runBlock(0.05);

    if (!(highWidth > lowWidth + 0.05)) {
        std::cerr << "Dynamic width did not grow with higher amplitude: low=" << lowWidth
                  << " high=" << highWidth << std::endl;
        return false;
    }
    if (!(widthBack < highWidth - 0.05)) {
        std::cerr << "Dynamic width did not relax after amplitude drop: back=" << widthBack
                  << " high=" << highWidth << std::endl;
        return false;
    }
    const double tol = 1.0e-9;
    if (std::abs(energyL - 1.0) > tol || std::abs(energyR - 1.0) > tol) {
        std::cerr << "Output weights lost unit energy: L=" << energyL << " R=" << energyR << std::endl;
        return false;
    }

    return true;
}

bool run_motion_moddepth_response() {
    t_kbeyond x{};
    x.setup_sr(48000.0);
    x.motion = 1.0;
    x.mix = 0.0;
    x.early = 0.0;
    x.predelay = 0.0;
    x.setup_predelay();
    x.regen = 0.0;
    x.decay = 0.0;
    x.update_decay();
    x.moddepth = 2.0;

    constexpr long frames = 64;
    std::array<double, frames> inL{};
    std::array<double, frames> inR{};
    std::array<double, frames> outL{};
    std::array<double, frames> outR{};
    double *ins[2] = {inL.data(), inR.data()};
    double *outs[2] = {outL.data(), outR.data()};

    auto runSineBlock = [&](double freq, double &phase) {
        const double inc = 2.0 * M_PI * freq / x.sr;
        for (long n = 0; n < frames; ++n) {
            double sample = std::sin(phase);
            phase += inc;
            if (phase > 2.0 * M_PI)
                phase -= 2.0 * M_PI;
            inL[n] = sample;
            inR[n] = sample;
        }
        kbeyond_perform64(&x, nullptr, ins, 2, outs, 2, frames, 0, nullptr);
        return x.debugModDepthTarget;
    };

    double phase = 0.0;
    double lowDepth = 0.0;
    for (int n = 0; n < 6; ++n)
        lowDepth = runSineBlock(200.0, phase);

    double highDepth = 0.0;
    for (int n = 0; n < 6; ++n)
        highDepth = runSineBlock(2000.0, phase);

    if (!(highDepth > lowDepth + 0.1)) {
        std::cerr << "Mod depth did not react to HF content: low=" << lowDepth
                  << " high=" << highDepth << std::endl;
        return false;
    }
    if (highDepth + 1.0e-9 < x.moddepth) {
        std::cerr << "HF modulation reduced below base moddepth" << std::endl;
        return false;
    }
    if (x.debugDampHFValue + 1.0e-9 < x.decayState.dampHF || x.debugDampMFValue + 1.0e-9 < x.decayState.dampMF) {
        std::cerr << "Dynamic damping dropped below baseline" << std::endl;
        return false;
    }

    return true;
}

bool run_mid_side_mix_normalization() {
    t_kbeyond x{};
    x.setup_sr(48000.0);

    constexpr double width = 1.0;
    constexpr double eps = 1.0e-12;
    auto apply_and_expect = [&](double widthNorm, double tailMid, double tailSide, double &expectedL, double &expectedR) -> bool {
        double clamped = std::max(0.0, std::min(widthNorm, 2.0));
        x.apply_width(clamped);
        double midCoeffL = 0.0;
        double sideCoeffL = 0.0;
        double midCoeffR = 0.0;
        double sideCoeffR = 0.0;
        for (int i = 0; i < t_kbeyond::N; ++i) {
            double mid = x.outMidBasis[i];
            double side = x.outWeightsSide[i];
            midCoeffL += x.outWeightsL[i] * mid;
            sideCoeffL += x.outWeightsL[i] * side;
            midCoeffR += x.outWeightsR[i] * mid;
            sideCoeffR += x.outWeightsR[i] * side;
        }
        const double tol = 1.0e-9;
        if (std::abs(x.outMidToL - midCoeffL) > tol || std::abs(x.outSideToL - sideCoeffL) > tol ||
            std::abs(x.outMidToR - midCoeffR) > tol || std::abs(x.outSideToR - sideCoeffR) > tol) {
            std::cerr << "Stored mid/side projection coefficients drifted from dot-product baseline" << std::endl;
            std::cerr << "midToL diff=" << std::abs(x.outMidToL - midCoeffL)
                      << " sideToL diff=" << std::abs(x.outSideToL - sideCoeffL)
                      << " midToR diff=" << std::abs(x.outMidToR - midCoeffR)
                      << " sideToR diff=" << std::abs(x.outSideToR - sideCoeffR) << std::endl;
            return false;
        }
        expectedL = tailMid * midCoeffL + tailSide * sideCoeffL;
        expectedR = tailMid * midCoeffR + tailSide * sideCoeffR;
        return true;
    };

    {
        const double tailMid = 1.0;
        const double tailSide = 0.05;

        double outL = 0.0;
        double outR = 0.0;
        double expectedL = 0.0;
        double expectedR = 0.0;
        if (!apply_and_expect(width, tailMid, tailSide, expectedL, expectedR))
            return false;
        x.mix_mid_side_to_lr(tailMid, tailSide, outL, outR);

        const double errL = std::abs(outL - expectedL);
        const double errR = std::abs(outR - expectedR);
        if (errL > 1.0e-9 || errR > 1.0e-9) {
            std::cerr << "Mid/side mix deviated from orthonormal expectation: L err=" << errL
                      << " R err=" << errR << std::endl;
            return false;
        }
    }

    {
        const double tailMid = 1.0e-3;
        const double tailSide = 1.0;

        double outL = 0.0;
        double outR = 0.0;
        double expectedL = 0.0;
        double expectedR = 0.0;
        if (!apply_and_expect(width, tailMid, tailSide, expectedL, expectedR))
            return false;
        x.mix_mid_side_to_lr(tailMid, tailSide, outL, outR);
        const double outputEnergy = outL * outL + outR * outR;
        const double expectedEnergy = expectedL * expectedL + expectedR * expectedR;
        const double energyDiffDb = 10.0 * std::log10((outputEnergy + eps) / (expectedEnergy + eps));
        if (std::abs(energyDiffDb) > 1.0e-6) {
            std::cerr << "Mid/side mix deviated from expected energy by " << energyDiffDb
                      << " dB" << std::endl;
            return false;
        }

        const double rmsExpected = std::sqrt(0.5 * expectedEnergy);
        const double rmsActual = std::sqrt(0.5 * outputEnergy);
        const double rmsDiffDb = 20.0 * std::log10((rmsActual + eps) / (rmsExpected + eps));
        if (std::abs(rmsDiffDb) > 1.0e-6) {
            std::cerr << "Mid/side mix RMS deviated from orthonormal target by " << rmsDiffDb
                      << " dB" << std::endl;
            return false;
        }
    }

    {
        const double tailMid = 0.1;
        const double tailSide = 0.8;
        const std::array<double, 3> widths{0.5, 1.0, 2.0};
        constexpr double leakEps = 1.0e-12;

        for (double widthNorm : widths) {
            double outL = 0.0;
            double outR = 0.0;
            double expectedL = 0.0;
            double expectedR = 0.0;
            if (!apply_and_expect(widthNorm, tailMid, tailSide, expectedL, expectedR))
                return false;
            x.mix_mid_side_to_lr(tailMid, tailSide, outL, outR);

            const double errL = std::abs(outL - expectedL);
            const double errR = std::abs(outR - expectedR);
            if (errL > 1.0e-9 || errR > 1.0e-9) {
                std::cerr << "Mid/side mix deviated from expected unitary projection at width=" << widthNorm
                          << ": L err=" << errL << " R err=" << errR << std::endl;
                return false;
            }

            const double expectedEnergyL = expectedL * expectedL;
            const double expectedEnergyR = expectedR * expectedR;
            const double energyL = outL * outL;
            const double energyR = outR * outR;
            const double rmsDiffLdB = 10.0 * std::log10((energyL + leakEps) / (expectedEnergyL + leakEps));
            const double rmsDiffRdB = 10.0 * std::log10((energyR + leakEps) / (expectedEnergyR + leakEps));
            if (std::abs(rmsDiffLdB) > 1.0e-9 || std::abs(rmsDiffRdB) > 1.0e-9) {
                std::cerr << "Mid/side mix per-channel RMS deviated beyond tolerance at width=" << widthNorm
                          << ": L diff=" << rmsDiffLdB << " dB R diff=" << rmsDiffRdB << " dB" << std::endl;
                return false;
            }
        }
    }

    return true;
}

bool run_width_projection_energy_sweep() {
    t_kbeyond x{};
    x.setup_sr(48000.0);

    const double tailMid = 0.6;
    const double tailSide = 0.4;
    const int steps = 256;
    double maxEnergyError = 0.0;
    double maxStereoImbalance = 0.0;

    for (int i = 0; i < steps; ++i) {
        double phase = static_cast<double>(i) / static_cast<double>(steps);
        double widthValue = 1.0 + std::sin(phase * 8.0 * M_PI);
        widthValue = std::max(0.0, std::min(widthValue, 2.0));
        x.apply_width(widthValue);

        double outL = 0.0;
        double outR = 0.0;
        x.mix_mid_side_to_lr(tailMid, tailSide, outL, outR);

        double energy = outL * outL + outR * outR;
        double expectedL = tailMid * x.outMidToL + tailSide * x.outSideToL;
        double expectedR = tailMid * x.outMidToR + tailSide * x.outSideToR;
        double expectedEnergy = expectedL * expectedL + expectedR * expectedR;
        maxEnergyError = std::max(maxEnergyError, std::abs(energy - expectedEnergy));

        std::array<double, t_kbeyond::N> fdnVec{};
        for (int n = 0; n < t_kbeyond::N; ++n)
            fdnVec[n] = tailMid * x.outMidBasis[n] + tailSide * x.outWeightsSide[n];

        double directL = 0.0;
        double directR = 0.0;
        for (int n = 0; n < t_kbeyond::N; ++n) {
            directL += fdnVec[n] * x.outWeightsL[n];
            directR += fdnVec[n] * x.outWeightsR[n];
        }

        if (std::abs(directL - outL) > 1.0e-12 || std::abs(directR - outR) > 1.0e-12) {
            std::cerr << "Direct projection mismatch under width sweep" << std::endl;
            return false;
        }

        double sideL = 0.0;
        double sideR = 0.0;
        x.mix_mid_side_to_lr(0.0, 1.0, sideL, sideR);
        double sideEnergyDiff = std::abs(sideL * sideL - sideR * sideR);
        maxStereoImbalance = std::max(maxStereoImbalance, sideEnergyDiff);

        double energyL = weights_energy(x.outWeightsL);
        double energyR = weights_energy(x.outWeightsR);
        if (std::abs(energyL - 1.0) > 1.0e-9 || std::abs(energyR - 1.0) > 1.0e-9) {
            std::cerr << "Width sweep broke unit-energy constraint: L=" << energyL << " R=" << energyR
                      << std::endl;
            return false;
        }
    }

    if (maxEnergyError > 1.0e-9) {
        std::cerr << "Width sweep energy deviated from analytic expectation: error=" << maxEnergyError << std::endl;
        return false;
    }
    if (maxStereoImbalance > 1.0e-9) {
        std::cerr << "Width sweep broke stereo balance: imbalance=" << maxStereoImbalance << std::endl;
        return false;
    }

    return true;
}

} // namespace motion_tests

