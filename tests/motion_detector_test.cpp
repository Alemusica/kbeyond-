#include "kbeyond_tilde.h"

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

} // namespace motion_tests

