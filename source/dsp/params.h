#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>

namespace kbeyond::dsp {

constexpr double kPredelayMaxSeconds = 0.5;
constexpr double kPredelaySafetySamples = 16.0;
constexpr std::size_t kFdnSize = 16;
constexpr std::size_t kEarlyTaps = 12;
constexpr std::size_t kLaserGroups = 3;
constexpr std::size_t kLaserTaps = 24;

enum class MixMode { Householder = 0, WHT, Hybrid };

struct AttributeDefaults {
    static constexpr double regen = 0.7;
    static constexpr double derez = 1.0;
    static constexpr double filter = 0.6;
    static constexpr double early = 0.3;
    static constexpr double focus = 1.0;
    static constexpr double predelay = 0.05; // seconds
    static constexpr double mix = 0.5;
    static constexpr double laser = 0.0;
    static constexpr double laserFocus = 0.55;
    static constexpr double laserGate = 0.35;
    static constexpr double laserWindow = 0.35;
    static constexpr double laserDiffusion = 0.65;
    static constexpr double width = 1.2;
    static constexpr double size = 0.6;
    static constexpr double color = 0.0;
    static constexpr double modrate = 0.15; // Hz
    static constexpr double moddepth = 3.0; // samples
    static constexpr double phiweight = 0.7;
    static constexpr double coherence = 0.8;
    static constexpr double uwalkRate = 0.25;
    static constexpr double decay = 0.0; // seconds
    static constexpr double dampLF = 0.2;
    static constexpr double dampMF = 0.0;
    static constexpr double dampHF = 0.6;
};

inline double clampd(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

inline long clampl(long v, long lo, long hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

inline double lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

inline double ms2samp(double ms, double sr) {
    return ms * 0.001 * sr;
}

inline double tiny_noise(uint32_t &rng) {
    rng ^= rng << 13;
    rng ^= rng >> 17;
    rng ^= rng << 5;
    return static_cast<double>(rng & 0xFFFFFF) * 1.0e-12 * (1.0 / 16777216.0);
}

} // namespace kbeyond::dsp

