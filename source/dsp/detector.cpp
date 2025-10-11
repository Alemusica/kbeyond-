#include "detector.h"

#include <algorithm>
#include <cmath>

namespace kbeyond::dsp {

namespace {
constexpr double kSlowWindowMs = 100.0;
constexpr double kFastWindowMs = 15.0;
constexpr double kRangeSlewMs = 45.0;
constexpr double kDopplerSlewMs = 25.0;
constexpr double kSpreadSlewMs = 50.0;
constexpr double kMinMs = 1.0;
constexpr double kHpCutHz = 700.0;
constexpr double kLpCutHz = 3000.0;
constexpr double kEpsilon = 1.0e-12;
}

double RangeDopplerDetector::BandFilter::process(double x) {
    lowState += hpCoef * (x - lowState);
    const double high = x - lowState;
    bandState += lpCoef * (high - bandState);
    return bandState;
}

void RangeDopplerDetector::BandFilter::reset() {
    lowState = 0.0;
    bandState = 0.0;
}

double RangeDopplerDetector::smoothing_from_ms(double sr, double ms) {
    const double srSafe = sr > 1.0 ? sr : 48000.0;
    const double tau = std::max(ms, kMinMs) * 0.001;
    const double denom = std::max(tau * srSafe, 1.0e-6);
    return std::exp(-1.0 / denom);
}

double RangeDopplerDetector::pole_from_hz(double sr, double hz) {
    if (sr <= 0.0 || hz <= 0.0)
        return 0.0;
    const double exponent = -2.0 * M_PI * hz / sr;
    return 1.0 - std::exp(exponent);
}

void RangeDopplerDetector::setSampleRate(double sr) {
    sr_ = sr > 1.0 ? sr : 48000.0;
    slowCoef_ = smoothing_from_ms(sr_, kSlowWindowMs);
    fastCoef_ = smoothing_from_ms(sr_, kFastWindowMs);
    rangeSlew_ = smoothing_from_ms(sr_, kRangeSlewMs);
    dopplerSlew_ = smoothing_from_ms(sr_, kDopplerSlewMs);
    spreadSlew_ = smoothing_from_ms(sr_, kSpreadSlewMs);
    const double hpCoef = pole_from_hz(sr_, kHpCutHz);
    const double lpCoef = pole_from_hz(sr_, kLpCutHz);
    midFilter_.hpCoef = hpCoef;
    midFilter_.lpCoef = lpCoef;
    sideFilter_.hpCoef = hpCoef;
    sideFilter_.lpCoef = lpCoef;
}

void RangeDopplerDetector::reset() {
    slowMid_ = 0.0;
    slowSide_ = 0.0;
    bandMid_ = 0.0;
    bandSide_ = 0.0;
    rangeRms_ = 0.0;
    rangeOut_ = 0.0;
    dopplerOut_ = 0.0;
    spreadOut_ = 0.5;
    midFilter_.reset();
    sideFilter_.reset();
}

void RangeDopplerDetector::process(double mid, double side) {
    const double midSq = mid * mid;
    const double sideSq = side * side;

    slowMid_ = slowMid_ * slowCoef_ + midSq * (1.0 - slowCoef_);
    slowSide_ = slowSide_ * slowCoef_ + sideSq * (1.0 - slowCoef_);

    const double rangeEnergy = slowMid_ + slowSide_;
    rangeRms_ = std::sqrt(std::max(rangeEnergy, 0.0));
    const double rangeNorm = std::clamp(rangeRms_, 0.0, 1.0);
    rangeOut_ = rangeOut_ * rangeSlew_ + rangeNorm * (1.0 - rangeSlew_);

    const double midBand = midFilter_.process(mid);
    const double sideBand = sideFilter_.process(side);

    bandMid_ = bandMid_ * fastCoef_ + midBand * midBand * (1.0 - fastCoef_);
    bandSide_ = bandSide_ * fastCoef_ + sideBand * sideBand * (1.0 - fastCoef_);

    const double bandEnergy = bandMid_ + bandSide_;
    const double denom = rangeEnergy + kEpsilon;
    const double dopplerRaw = std::clamp(bandEnergy / denom, 0.0, 1.0);
    dopplerOut_ = dopplerOut_ * dopplerSlew_ + dopplerRaw * (1.0 - dopplerSlew_);

    const double spreadRaw = std::clamp(slowSide_ / denom, 0.0, 1.0);
    spreadOut_ = spreadOut_ * spreadSlew_ + spreadRaw * (1.0 - spreadSlew_);
}

} // namespace kbeyond::dsp

