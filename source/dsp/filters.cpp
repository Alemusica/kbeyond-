#include "filters.h"

#include <cmath>

namespace kbeyond::dsp {

void Tilt::set(double coef) {
    a_ = clampd(coef, -0.999, 0.999);
}

double Tilt::process(double x) {
    const double y = x + a_ * (x - z_);
    z_ = x;
    return y;
}

void OnePoleLP::setCutoff(double sr, double cutoffHz) {
    cutoffHz = clampd(cutoffHz, 1.0, sr * 0.45);
    const double x = std::exp(-2.0 * M_PI * cutoffHz / sr);
    a_ = 1.0 - x;
    b_ = x;
}

double OnePoleLP::process(double x) {
    z_ = a_ * x + b_ * z_;
    return z_;
}

} // namespace kbeyond::dsp

