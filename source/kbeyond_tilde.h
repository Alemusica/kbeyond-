#pragma once
// kbeyond_tilde.h
// Max/MSP external: "kbeyond~" — a spacious FDN reverb inspired by kBeyond
// MIT License.
//
// Build: Max 8 SDK (msp), C++17, x64
// Object name: kbeyond~
// Inlets:  2 (signal L/R)
// Outlets: 2 (signal L/R)
// Attributes (0..1 unless specified):
//   @regen, @decay(seconds), @derez, @filter, @early, @predelay (0..0.5 seconds), @mix,
//   @width(0..2), @size(0..1), @color(-1..+1), @modrate(Hz), @moddepth(samples),
//   @damplf, @dampmf, @damphf, @phiweight(0..1)
// Notes:
//  - Delay network is a 16x16 FDN with Householder(φ) mixing (computed on the fly).
//  - Early reflections are φ-spaced; width via M/S shuffler.
//  - "derez" here controls pre-FDN bandwidth.

#ifndef KBEYOND_UNIT_TEST
extern "C" {
#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
}

#ifndef C74_EXPORT
#define C74_EXPORT
#endif
#else
struct t_pxobject {};
struct t_symbol {};
struct t_atom {};
struct t_object {};
struct t_class {};
using t_max_err = int;

#ifndef C74_EXPORT
#define C74_EXPORT
#endif
#endif

#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <algorithm>

struct t_kbeyond;

void *kbeyond_new(t_symbol *s, long argc, t_atom *argv);
void kbeyond_free(t_kbeyond *x);
void kbeyond_assist(t_kbeyond *x, void *b, long m, long a, char *s);
void kbeyond_dsp64(t_kbeyond *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void kbeyond_perform64(t_kbeyond *x, t_object *dsp64, double **ins, long nin, double **outs, long nout, long sampleframes, long flags, void *userparam);
void C74_EXPORT ext_main(void *r);

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
        size_t len = buf.size();
        if (len == 0) return 0.0;
        double readPos = (double)w - delaySamp;
        readPos -= std::floor(readPos / (double)len) * (double)len;
        long i0 = (long)readPos;
        long i1 = (i0 + 1) % (long)len;
        double frac = readPos - (double)i0;
        return buf[(size_t)i0] + (buf[(size_t)i1] - buf[(size_t)i0]) * frac;
    }
    inline double readInt(long delaySamp) const {
        size_t len = buf.size();
        if (len == 0) return 0.0;
        long rp = w - delaySamp;
        while (rp < 0) rp += (long)len;
        return buf[(size_t)(rp % (long)len)];
    }
    inline void write(double x) {
        if (buf.empty()) return;
        buf[(size_t)w] = x;
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
        double y = x + a * (x - z); // Householder tail tone tilt
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
    void reset() { z = 0.0; }
};

// ------------------------------- Main object
struct t_kbeyond {
    static constexpr double kPredelayMaxSeconds = 0.5;
    static constexpr double kPredelaySafetySamples = 16.0;

    t_pxobject      ob;

    // Parameters
    double regen     = 0.7;
    double derez     = 1.0;
    double filter    = 0.6;
    double early     = 0.3;
    double predelay  = 0.05; // seconds
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
    DelayLine earlyBufMid;
    DelayLine earlyBufSide;
    std::array<long, kEarlyTaps> earlyDel {};
    std::array<double, kEarlyTaps> earlyGain {};
    std::array<double, kEarlyTaps> earlyCos {};
    std::array<double, kEarlyTaps> earlySin {};

    enum class MixMode { Householder = 0, WHT, Hybrid };

    // FDN
    static const int N = 16;
    std::array<DelayLine, N> fdn;
    std::array<double, N>    fdn_len {};
    std::array<double, N>    fdn_read {};
    std::array<double, N>    fdn_phase {};
    std::array<double, N>    fdn_phaseInc {};
    std::array<double, N>    fdn_out {};
    std::array<double, N>    fdn_fb {};
    std::array<Tilt, N>      fdn_tilt {};
    std::array<OnePoleLP, N> fdn_lp {};
    std::array<OnePoleLP, N> fdn_low {};
    std::array<OnePoleLP, N> fdn_high {};
    std::array<double, N>    inWeights {};

    // Output mapping
    std::array<double, N> outBaseMid {};
    std::array<double, N> outBaseSide {};
    std::array<double, N> outWeightsL {};
    std::array<double, N> outWeightsR {};

    // Householder vector u (normalized)
    std::array<double, N> u {};
    std::array<double, N> diffusionScratch {};

    MixMode   mixMode = MixMode::Householder;
    t_symbol *modeMixSym = nullptr;

    // Decay / damping
    double decay    = 0.0;  // seconds, 0 disables RT60 mapping
    double dampLF   = 0.2;
    double dampMF   = 0.0;
    double dampHF   = 0.6;
    double dampLF_mul = 1.0;
    double dampMF_mul = 1.0;
    double dampHF_mul = 1.0;
    std::array<double, N> fdn_decay {};

    // RNG for tiny noise to avoid denormals
    uint32_t rng = 0x1234567u;

    // Methods
    void setup_sr(double newsr);
    void setup_predelay();
    void setup_early();
    void setup_fdn();
    void refresh_filters();
    void update_diffusion();
    void update_output_weights();
    void update_modulators();
    void update_decay();
    void apply_diffusion(const std::array<double, N> &input, std::array<double, N> &output);
    void render_early(double inL, double inR, double widthNorm, double earlyAmt, double &earlyL, double &earlyR);
    inline double tiny() {
        rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
        return (double)(rng & 0xFFFFFF) * 1.0e-12 * (1.0/16777216.0);
    }
};

