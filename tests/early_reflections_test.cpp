#include "kbeyond_tilde.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
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
        return 1;
    }

    return 0;
}
