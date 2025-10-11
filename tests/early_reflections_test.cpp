#include "kbeyond_tilde.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool test_early_reflections() {
    t_kbeyond x{};
    x.setup_sr(48000.0);
    x.width = 1.0;
    x.early = 1.0;
    x.predelay = 0.0;
    x.setup_predelay();

    constexpr int frames = 4096;
    std::vector<double> left(frames, 0.0);
    std::vector<double> right(frames, 0.0);
    for (int n = 0; n < frames; ++n) {
        double midIn = (n == 0) ? 1.0 : 0.0;
        double l = 0.0;
        double r = 0.0;
        x.render_early(midIn, x.width, x.early, l, r);
        left[n] = l;
        right[n] = r;
    }

    double peakL = 0.0;
    double peakR = 0.0;
    for (int n = 0; n < frames; ++n) {
        peakL = std::max(peakL, std::abs(left[n]));
        peakR = std::max(peakR, std::abs(right[n]));
    }

    const double eps = 1.0e-12;
    double diffDb = 20.0 * std::log10((peakL + eps) / (peakR + eps));

    if (std::abs(diffDb) > 0.1) {
        std::cerr << "Mono impulse early response imbalance: " << diffDb << " dB" << std::endl;
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

} // namespace

int main() {
    if (!test_early_reflections())
        return 1;
    if (!test_mix_modes())
        return 2;
    return 0;
}
