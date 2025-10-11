#pragma once

#include <array>
#include <vector>
#include <cstdint>

#include "buffers.h"
#include "params.h"
#ifdef KBEYOND_UNIT_TEST
#include <cstddef>
#endif

namespace kbeyond::dsp {

struct EarlySection {
    void setup(double sr,
               double size,
               double focus,
               const std::vector<double> &tapPattern,
               const std::vector<double> &laserMidPattern,
               const std::vector<double> &laserSidePattern);

    void resetState();
    void updateGate(double gate);
    void updateWindow(double windowSeconds, double sr);
    void updatePhaseIncrement(double focus, double sr);
    void updateEnvelopeCoefficients(double sr);

    void render(double inL,
                double inR,
                double widthNorm,
                double earlyAmt,
                double focusAmt,
                double clusterAmt,
                double &earlyL,
                double &earlyR,
                uint32_t &rng);

    double excite() const { return laserExcite_; }

    double computeQSwitchMix(double clusterAmt, double diffusion);

#ifdef KBEYOND_UNIT_TEST
    const std::array<long, kEarlyTaps> &debugTapDelays() const { return earlyDel_; }
    const std::array<double, kEarlyTaps> &debugTapGains() const { return earlyGain_; }
    const std::array<double, kEarlyTaps> &debugTapCos() const { return earlyCos_; }
    const std::array<double, kEarlyTaps> &debugTapSin() const { return earlySin_; }
#endif

private:
    DelayLine earlyBufMid_{};
    DelayLine earlyBufSide_{};

    std::array<long, kEarlyTaps> earlyDel_{};
    std::array<double, kEarlyTaps> earlyGain_{};
    std::array<double, kEarlyTaps> earlyCos_{};
    std::array<double, kEarlyTaps> earlySin_{};

#if KBEYOND_ENABLE_LASER
    std::array<double, kLaserTaps> laserDelay_{};
    std::array<double, kLaserTaps> laserMidGain_{};
    std::array<double, kLaserTaps> laserSideGain_{};
    std::array<double, kLaserTaps> laserCos_{};
    std::array<double, kLaserTaps> laserSin_{};
    std::array<double, kLaserTaps> laserShape_{};

    double laserEnv_ = 0.0;
    double laserExcite_ = 0.0;
    double laserGateScaled_ = 0.02;
    double laserPhase_ = 0.0;
    double laserPhaseInc_ = 0.0;
    double laserEnvAttack_ = 0.0;
    double laserEnvRelease_ = 0.0;
    double qswitchEnv_ = 0.0;
    double qswitchAttack_ = 0.0;
    double qswitchRelease_ = 0.0;
    long qswitchWindowSamples_ = 0;
    long qswitchCounter_ = 0;
#else
    double laserExcite_ = 0.0;
#endif
    double sampleRate_ = 48000.0;
};

} // namespace kbeyond::dsp

