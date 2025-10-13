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
//   @regen, @decay(seconds), @derez, @filter, @early, @predelay (0..0.5 seconds, runtime ≥1 sample), @mix,
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
struct t_pxobject {
    int z_disabled = 0;
};
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

#ifndef KBEYOND_ENABLE_LASER
#define KBEYOND_ENABLE_LASER 0
#endif

#include "dsp/prime_modes.h"
#include "dsp/params.h"
#include "dsp/buffers.h"
#include "dsp/filters.h"
#include "dsp/early.h"
#include "dsp/fdn.h"
#include "dsp/decay.h"
#include "dsp/mod.h"
#include "dsp/detector.h"

struct t_kbeyond;

void *kbeyond_new(t_symbol *s, long argc, t_atom *argv);
void kbeyond_free(t_kbeyond *x);
void kbeyond_assist(t_kbeyond *x, void *b, long m, long a, char *s);
void kbeyond_dsp64(t_kbeyond *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void kbeyond_perform64(t_kbeyond *x, t_object *dsp64, double **ins, long nin, double **outs, long nout, long sampleframes, long flags, void *userparam);
void C74_EXPORT ext_main(void *r);

// ------------------------------- Helpers
// ------------------------------- Main object
struct t_kbeyond {
    static constexpr double kPredelayMaxSeconds = kbeyond::dsp::kPredelayMaxSeconds;
    static constexpr double kPredelaySafetySamples = kbeyond::dsp::kPredelaySafetySamples;

    t_pxobject      ob;

    // Parameters
    double regen     = kbeyond::dsp::AttributeDefaults::regen;
    double derez     = kbeyond::dsp::AttributeDefaults::derez;
    double filter    = kbeyond::dsp::AttributeDefaults::filter;
    double early     = kbeyond::dsp::AttributeDefaults::early;
    double focus     = kbeyond::dsp::AttributeDefaults::focus;
    double predelay  = kbeyond::dsp::AttributeDefaults::predelay;
    double mix       = kbeyond::dsp::AttributeDefaults::mix;

    double laser          = kbeyond::dsp::AttributeDefaults::laser;
    double laserFocus     = kbeyond::dsp::AttributeDefaults::laserFocus;
    double laserGate      = kbeyond::dsp::AttributeDefaults::laserGate;
    double laserWindow    = kbeyond::dsp::AttributeDefaults::laserWindow;
    double laserDiffusion = kbeyond::dsp::AttributeDefaults::laserDiffusion;

    double width     = kbeyond::dsp::AttributeDefaults::width;
    double size      = kbeyond::dsp::AttributeDefaults::size;
    double color     = kbeyond::dsp::AttributeDefaults::color;
    double modrate   = kbeyond::dsp::AttributeDefaults::modrate;
    double moddepth  = kbeyond::dsp::AttributeDefaults::moddepth;
    double phiweight = kbeyond::dsp::AttributeDefaults::phiweight;
    double coherence = kbeyond::dsp::AttributeDefaults::coherence;
    double uwalkRate = kbeyond::dsp::AttributeDefaults::uwalkRate;
    double motion    = kbeyond::dsp::AttributeDefaults::motion;

    // SR / vector
    double sr = 48000.0;
    long   vs = 64;

    // Predelay
    kbeyond::dsp::DelayLine predL, predR;
    long      maxPred = 16384;
    double    predSamps = 0.0;

    kbeyond::dsp::EarlySection earlySection;

    using MixMode = kbeyond::dsp::MixMode;

    // FDN
    static const int N = static_cast<int>(kbeyond::dsp::kFdnSize);
    static const int kEarlyTaps = static_cast<int>(kbeyond::dsp::kEarlyTaps);
    static const int kLaserTaps = static_cast<int>(kbeyond::dsp::kLaserTaps);
    kbeyond::dsp::FdnState fdnState;
    std::array<kbeyond::dsp::Tilt, N>      fdn_tilt {};
    std::array<kbeyond::dsp::OnePoleLP, N> fdn_lp {};
    std::array<kbeyond::dsp::OnePoleLP, N> fdn_low {};
    std::array<kbeyond::dsp::OnePoleLP, N> fdn_high {};
    std::array<double, N>    inWeights {};

    // Output mapping
    std::array<double, N> outBaseMid {};
    std::array<double, N> outBaseSide {};
    std::array<double, N> outMidBasis {};
    std::array<double, N> outWeightsSide {};
    std::array<double, N> outWeightsL {};
    std::array<double, N> outWeightsR {};
    double widthApplied = kbeyond::dsp::AttributeDefaults::width;
    double outMidToL = 0.0;
    double outSideToL = 0.0;
    double outMidToR = 0.0;
    double outSideToR = 0.0;

    // Householder vector u (normalized)
    std::array<double, N> u {};
    std::array<double, N> diffusionScratch {};

    MixMode   mixMode = MixMode::Householder;
    t_symbol *modeMixSym = nullptr;

    prime_modes::Pattern modeER = prime_modes::Pattern::Aureo;
    prime_modes::Pattern modeLate = prime_modes::Pattern::PrimeAureo;
    prime_modes::Pattern modeMid = prime_modes::Pattern::Prime;
    prime_modes::Pattern modeSide = prime_modes::Pattern::Aureo;
    long patternSeed = 1337;
    t_symbol *modeERSym = nullptr;
    t_symbol *modeLateSym = nullptr;
    t_symbol *modeMidSym = nullptr;
    t_symbol *modeSideSym = nullptr;

    // Decay / damping
    double decay    = 0.0;  // seconds, 0 disables RT60 mapping
    double dampLF   = 0.2;
    double dampMF   = 0.0;
    double dampHF   = 0.6;
    kbeyond::dsp::DecayState decayState {};

    // RNG for tiny noise to avoid denormals
    uint32_t rng = 0x1234567u;

    kbeyond::dsp::ModState modState {};
    kbeyond::dsp::RangeDopplerDetector motionDetector {};
    double rangeEnv = 0.0;
    double dopplerEnv = 0.0;
    double spreadEnv = 0.5;
    double dryEnergyEnv = 0.0;
    double wetEnergyEnv = 0.0;
    double dryLevelEnv = 0.0;
    double wetLevelEnv = 0.0;
    double wetMakeup = 1.0;

    // Methods
    void setup_sr(double newsr);
    void setup_predelay();
    void setup_early();
    void setup_fdn();
    void update_laser_gate();
    void update_laser_window();
    void update_laser_phase_inc();
    void update_laser_envelope();
    void refresh_filters();
    void update_diffusion();
    void update_output_weights();
    void update_modulators();
    void update_decay();
    void update_quantum_walk();
    void reset_quantum_walk();
    void reset_tail_gain_state();
    void apply_quantum_dither(std::array<double, N> &vector);
    void apply_diffusion(const std::array<double, N> &input, std::array<double, N> &output);
    void apply_quantum_walk(std::array<double, N> &feedback);
    void update_injection_weights();
    void apply_width(double widthNorm);
    void update_output_basis();
    void mix_mid_side_to_lr(double tailMid, double tailSide, double &outL, double &outR) const;
    std::vector<double> make_pattern(prime_modes::Pattern mode, std::size_t count, std::uint32_t salt) const;

#ifdef KBEYOND_UNIT_TEST
    void render_early(double inL,
                      double inR,
                      double widthValue,
                      double earlyAmt,
                      double focusAmt,
                      double &outL,
                      double &outR);
    double debugWidthTarget = 0.0;
    double debugModDepthTarget = 0.0;
    double debugDampMFValue = 0.0;
    double debugDampHFValue = 0.0;
#endif
};

