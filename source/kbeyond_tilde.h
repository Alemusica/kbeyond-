#pragma once
// kbeyond_tilde.h
// Max/MSP external: "kbeyond~" — a spacious FDN reverb inspired by kBeyond
// MIT-like license for this example. Use at your own risk.
//
// Build: Max 8 SDK (msp), C++17 (or C++11), x64
// Object name: kbeyond~
// Inlets:  2 (signal L/R)
// Outlets: 2 (signal L/R)
// Attributes (0..1 unless specified):
//   @regen, @derez, @filter, @early, @predelay, @mix,
//   @width(0.5..2.0), @size(0..1), @color(-1..+1), @modrate(Hz), @moddepth(samples),
//   @phiweight(0..1)
// Notes:
//  - Delay network is a 16x16 FDN with Householder(φ) mixing (computed on the fly).
//  - Early reflections are φ-spaced; width via M/S shuffler.
//  - "derez" here controls pre-FDN bandwidth (brickwall-lite). A more advanced
//    decimator+Bézier+sinc can be dropped in later if desired.

extern "C" {
#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
}

#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <random>

// Forward decls
typedef struct _kbeyond t_kbeyond;
void *kbeyond_new(t_symbol *s, long argc, t_atom *argv);
void kbeyond_free(t_kbeyond *x);
void kbeyond_assist(t_kbeyond *x, void *b, long m, long a, char *s);
void kbeyond_dsp64(t_kbeyond *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void kbeyond_perform64(t_kbeyond *x, t_object *dsp64, double **ins, long nin, double **outs, long nout, long sampleframes, long flags, void *userparam);
void ext_main(void *r);

// ------------------------------- Helpers
static inline double clampd(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}
static inline long clampl(long v, long lo, long hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}
static inline double lerp(double a, double b, double t) {
    return a + (b - a) * t;
}
static inline double ms2samp(double ms, double sr) { return ms * 0.001 * sr; }

// ------------------------------- DSP core state
struct DelayLine {
    std::vector<double> buf;
    long w = 0;
    DelayLine() {}
    void setup(size_t len) { buf.assign(len, 0.0); w = 0; }
    inline double readFrac(double delaySamp) const {
        // Read delay using linear interpolation, delaySamp in [1, len-2]
        size_t len = buf.size();
        double readPos = (double)w - delaySamp;
        // wrap to [0, len)
        readPos -= std::floor(readPos / (double)len) * (double)len;
        long i0 = (long)readPos;
        long i1 = (i0 + 1) % len;
        double frac = readPos - (double)i0;
        return buf[i0] + (buf[i1] - buf[i0]) * frac;
    }
    inline double readInt(long delaySamp) const {
        size_t len = buf.size();
        long rp = w - delaySamp;
        while (rp < 0) rp += (long)len;
        return buf[rp % len];
    }
    inline void write(double x) {
        buf[w] = x;
        w = (w + 1) % (long)buf.size();
    }
    inline size_t size() const { return buf.size(); }
};

// Small utility: one-pole tone tilt (acts as gentle high/low tilt depending on sign)
struct Tilt {
    double z = 0.0;
    double a = 0.0;
    void set(double coef) { a = clampd(coef, -0.999, 0.999); }
    inline double process(double x) {
        // y = x + a*(x - z); z = x; => high-shelf-ish tilt around Nyquist/4-ish
        double y = x + a * (x - z);
        z = x;
        return y;
    }
};

struct OnePoleLP {
    double a = 0.0, b = 0.0, z = 0.0;
    void setCutoff(double sr, double cutoffHz) {
        cutoffHz = clampd(cutoffHz, 1.0, sr * 0.45);
        double x = std::exp(-2.0 * M_PI * cutoffHz / sr);
        a = 1.0 - x;
        b = x;
    }
    inline double process(double x) {
        z = a * x + b * z;
        return z;
    }
};

// ------------------------------- Main object
struct _kbeyond {
    t_pxobject      ob;

    // Parameters
    double regen     = 0.7;
    double derez     = 1.0;
    double filter    = 0.6;
    double early     = 0.3;
    double predelay  = 0.05;
    double mix       = 0.5;

    double width     = 1.2;
    double size      = 0.6;
    double color     = 0.0;
    double modrate   = 0.15; // Hz
    double moddepth  = 3.0;  // samples
    double phiweight = 0.7;  // 0..1

    // SR / vector
    double sr = 48000.0;
    long   vs = 64;

    // Predelay
    DelayLine predL, predR;
    long      maxPred = 16384;
    double    predSamps = 0.0;

    // Early reflections
    static const int kEarlyTaps = 12;
    DelayLine earlyBufL, earlyBufR;
    std::array<long, kEarlyTaps> earlyDelL {};
    std::array<long, kEarlyTaps> earlyDelR {};
    std::array<double, kEarlyTaps> earlyGain {};

    // FDN
    static const int N = 16;
    std::array<DelayLine, N> fdn;
    std::array<double, N>    fdn_len {};          // nominal lengths (samples)
    std::array<double, N>    fdn_read {};         // current read length with modulation
    std::array<double, N>    fdn_phase {};        // LFO phase
    std::array<double, N>    fdn_phaseInc {};     // LFO increment
    std::array<double, N>    fdn_out {};          // outputs of lines (before mix)
    std::array<double, N>    fdn_fb {};           // feedback vector (after Householder)
    std::array<Tilt, N>      fdn_tilt {};         // per-line tilt in feedback
    std::array<OnePoleLP, N> fdn_lp {};           // per-line LP for "filter/derez"

    // Output mapping
    std::array<double, N> outWeightsL {};
    std::array<double, N> outWeightsR {};

    // Householder vector u (normalized)
    std::array<double, N> u {};

    // RNG for tiny noise to avoid denormals
    uint32_t rng = 0x1234567u;

    // Methods
    void setup_sr(double newsr);
    void setup_predelay();
    void setup_early();
    void setup_fdn();
    void update_householder();
    void update_output_weights();
    void update_modulators();
    inline double tiny() {
        rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
        return (double)(rng & 0xFFFFFF) * 1.0e-12 * (1.0/16777216.0);
    }
};
