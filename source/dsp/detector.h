#pragma once

#include <cstddef>

namespace kbeyond::dsp {

class RangeDopplerDetector {
public:
    void setSampleRate(double sr);
    void reset();
    void process(double mid, double side);

    double range() const { return rangeOut_; }
    double doppler() const { return dopplerOut_; }
    double spread() const { return spreadOut_; }
    double slowEnergy() const { return slowMid_ + slowSide_; }
    double rangeRms() const { return rangeRms_; }

private:
    struct BandFilter {
        double lowState = 0.0;
        double bandState = 0.0;
        double hpCoef = 0.0;
        double lpCoef = 0.0;

        double process(double x);
        void reset();
    };

    static double smoothing_from_ms(double sr, double ms);
    static double pole_from_hz(double sr, double hz);

    double sr_ = 48000.0;
    double slowCoef_ = 0.0;
    double fastCoef_ = 0.0;
    double rangeSlew_ = 0.0;
    double dopplerSlew_ = 0.0;
    double spreadSlew_ = 0.0;

    double slowMid_ = 0.0;
    double slowSide_ = 0.0;
    double bandMid_ = 0.0;
    double bandSide_ = 0.0;

    double rangeRms_ = 0.0;
    double rangeOut_ = 0.0;
    double dopplerOut_ = 0.0;
    double spreadOut_ = 0.5;

    BandFilter midFilter_{};
    BandFilter sideFilter_{};
};

} // namespace kbeyond::dsp

