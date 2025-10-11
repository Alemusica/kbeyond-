#pragma once

#include "params.h"

namespace kbeyond::dsp {

class Tilt {
public:
    void set(double coef);
    double process(double x);
    void reset() { z_ = 0.0; }

private:
    double z_ = 0.0;
    double a_ = 0.0;
};

class OnePoleLP {
public:
    void setCutoff(double sr, double cutoffHz);
    double process(double x);
    void reset() { z_ = 0.0; }

private:
    double a_ = 0.0;
    double b_ = 0.0;
    double z_ = 0.0;
};

} // namespace kbeyond::dsp

