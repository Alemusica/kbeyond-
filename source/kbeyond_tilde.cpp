#include "kbeyond_tilde.h"
#include "householder_phi16.h"
#include "dsp/mixing.h"

#include <algorithm>
#include <numeric>
#include <cstdio>
#include <cstddef>
#include <cmath>

static constexpr std::size_t kAssistStringMax = 256;

static t_class* s_kbeyond_class = nullptr;

#ifndef KBEYOND_UNIT_TEST
namespace {

double atom_to_double(long argc, t_atom* argv, double fallback) {
    if (argc <= 0 || argv == nullptr)
        return fallback;
    if (atom_gettype(argv) == A_LONG)
        return static_cast<double>(atom_getlong(argv));
    if (atom_gettype(argv) == A_FLOAT)
        return atom_getfloat(argv);
    return fallback;
}

} // namespace
#endif

using kbeyond::dsp::clampd;
using kbeyond::dsp::clampl;
using kbeyond::dsp::lerp;
using kbeyond::dsp::ms2samp;

void t_kbeyond::setup_sr(double newsr) {
    sr = newsr > 1.0 ? newsr : 48000.0;
    vs = std::max<long>(64, vs);
    setup_predelay();
    setup_early();
    update_laser_gate();
    update_laser_window();
    update_laser_phase_inc();
    update_laser_envelope();
    setup_fdn();
    refresh_filters();
    update_diffusion();
    update_output_weights();
    update_modulators();
    update_quantum_walk();
    motionDetector.setSampleRate(sr);
    motionDetector.reset();
    rangeEnv = motionDetector.range();
    dopplerEnv = motionDetector.doppler();
    spreadEnv = motionDetector.spread();
    dryEnergyEnv = 0.0;
    wetEnergyEnv = 0.0;
    dryLevelEnv = 0.0;
    wetLevelEnv = 0.0;
    wetMakeup = 1.0;
    apply_width(clampd(width, 0.0, 2.0));
}

std::vector<double> t_kbeyond::make_pattern(prime_modes::Pattern mode, std::size_t count, std::uint32_t salt) const {
    std::uint32_t base = patternSeed >= 0 ? static_cast<std::uint32_t>(patternSeed)
                                          : static_cast<std::uint32_t>(-patternSeed);
    std::uint32_t mix = salt * 0x9E3779B9u + 0xA511E9B5u;
    std::uint32_t seed = base ^ mix;
    if (seed == 0)
        seed = 0x9E3779B9u;
    return prime_modes::generate_pattern(mode, count, seed);
}

void t_kbeyond::setup_predelay() {
    double maxSamples = kPredelayMaxSeconds * sr;
    double desired = std::ceil(maxSamples + kPredelaySafetySamples);
    maxPred = (long)std::max<double>(desired, 256.0);
    size_t bufferLen = (size_t)maxPred + 8;
    predL.setup(bufferLen);
    predR.setup(bufferLen);
    double maxDelay = std::max(0.0, (double)bufferLen - 4.0);
    predSamps = clampd(predelay * sr, 0.0, maxDelay);
    if (predSamps < 1.0)
        predSamps = 1.0;
}

void t_kbeyond::setup_early() {
    int pairs = kEarlyTaps / 2;
    int patternCount = pairs;
    if (kEarlyTaps % 2 != 0)
        ++patternCount;
    if (patternCount <= 0)
        patternCount = 1;
    auto tapPattern = make_pattern(modeER, static_cast<std::size_t>(patternCount), 0x1001u);
    std::vector<double> laserMidPattern;
    std::vector<double> laserSidePattern;
#if KBEYOND_ENABLE_LASER
    laserMidPattern = make_pattern(modeMid, kbeyond::dsp::kLaserTaps, 0x5001u);
    laserSidePattern = make_pattern(modeSide, kbeyond::dsp::kLaserTaps, 0x5002u);
#endif
    earlySection.setup(sr, size, laserFocus, tapPattern, laserMidPattern, laserSidePattern);
    earlySection.resetState();
}

#ifdef KBEYOND_UNIT_TEST
void t_kbeyond::render_early(double inL,
                             double inR,
                             double widthValue,
                             double earlyAmt,
                             double focusAmt,
                             double &outL,
                             double &outR) {
    double clusterAmt = clampd(laser, 0.0, 1.0);
    earlySection.render(inL,
                        inR,
                        clampd(widthValue, 0.0, 2.0),
                        clampd(earlyAmt, 0.0, 1.0),
                        clampd(focusAmt, 0.0, 1.0),
                        clusterAmt,
                        outL,
                        outR,
                        rng);
}
#endif

void t_kbeyond::update_laser_gate() {
#if KBEYOND_ENABLE_LASER
    earlySection.updateGate(laserGate);
#endif
}

void t_kbeyond::update_laser_window() {
#if KBEYOND_ENABLE_LASER
    earlySection.updateWindow(laserWindow, sr);
#endif
}

void t_kbeyond::update_laser_phase_inc() {
#if KBEYOND_ENABLE_LASER
    earlySection.updatePhaseIncrement(laserFocus, sr);
#endif
}

void t_kbeyond::update_laser_envelope() {
#if KBEYOND_ENABLE_LASER
    earlySection.updateEnvelopeCoefficients(sr);
#endif
}

void t_kbeyond::setup_fdn() {
    auto latePattern = make_pattern(modeLate, N, 0x2001u);
    kbeyond::dsp::setup_fdn(fdnState, sr, size, clampd(moddepth, 0.0, 32.0), latePattern);
    kbeyond::dsp::init_fdn_phases(modState);
    update_decay();
    reset_quantum_walk();
    update_injection_weights();
}

void t_kbeyond::refresh_filters() {
    double hf = lerp(1200.0, sr * 0.45, clampd(filter, 0.0, 1.0));
    double derezCut = lerp(1600.0, sr * 0.38, clampd(derez, 0.0, 1.0));
    double cutoff = std::min(hf, derezCut);
    double tiltBase = clampd(color, -1.0, 1.0) * 0.35;
    for (int i = 0; i < N; ++i) {
        double idx = (double)i / (double)std::max(1, N - 1);
        double tilt = tiltBase * (idx * 2.0 - 1.0);
        fdn_tilt[i].set(tilt);
        double bias = lerp(0.75, 1.25, idx);
        fdn_lp[i].setCutoff(sr, cutoff * bias);
    }
}

void t_kbeyond::update_injection_weights() {
    auto midPattern = make_pattern(modeMid, N, 0x3001u);
    for (int i = 0; i < N; ++i) {
        double angle = (2.0 * M_PI * (double)i) / (double)N;
        double base = 0.9 / (double)N * (1.0 + 0.35 * std::sin(angle * 0.73 + 0.2));
        double scale = (i < (int)midPattern.size()) ? lerp(0.65, 1.45, midPattern[(std::size_t)i]) : 1.0;
        inWeights[i] = base * scale;
    }
}

void t_kbeyond::update_diffusion() {
    make_phi_vector<N>(u, clampd(phiweight, 0.0, 1.0));
}

void t_kbeyond::update_output_basis() {
    double midNorm = 0.0;
    for (int i = 0; i < N; ++i) {
        double mid = outBaseMid[i];
        midNorm += mid * mid;
    }

    if (midNorm <= 0.0) {
        std::fill(outMidBasis.begin(), outMidBasis.end(), 0.0);
        outMidBasis[0] = 1.0;
    } else {
        double invMid = 1.0 / std::sqrt(midNorm);
        for (int i = 0; i < N; ++i)
            outMidBasis[i] = outBaseMid[i] * invMid;
    }

    double dot = 0.0;
    for (int i = 0; i < N; ++i)
        dot += outMidBasis[i] * outBaseSide[i];

    double sideNorm = 0.0;
    for (int i = 0; i < N; ++i) {
        double rawSide = outBaseSide[i] - dot * outMidBasis[i];
        outWeightsSide[i] = rawSide;
        sideNorm += rawSide * rawSide;
    }

    if (sideNorm <= 0.0) {
        std::fill(outWeightsSide.begin(), outWeightsSide.end(), 0.0);
        int idx0 = 0;
        int idx1 = (N > 1) ? 1 : 0;
        outWeightsSide[idx0] = outMidBasis[idx1];
        if (idx1 != idx0)
            outWeightsSide[idx1] = -outMidBasis[idx0];
        else
            outWeightsSide[idx0] = 1.0;
        sideNorm = 0.0;
        for (double v : outWeightsSide)
            sideNorm += v * v;
    }

    double invSide = sideNorm > 0.0 ? 1.0 / std::sqrt(sideNorm) : 1.0;
    for (double &v : outWeightsSide)
        v *= invSide;
}

void t_kbeyond::apply_width(double widthNorm) {
    update_output_basis();
    widthNorm = clampd(widthNorm, 0.0, 2.0);
    double normL = 0.0;
    double normR = 0.0;
    for (int i = 0; i < N; ++i) {
        double mid = outMidBasis[i];
        double side = outWeightsSide[i];
        double wL = mid + widthNorm * side;
        double wR = mid - widthNorm * side;
        outWeightsL[i] = wL;
        outWeightsR[i] = wR;
        normL += wL * wL;
        normR += wR * wR;
    }
    double scaleL = normL > 0.0 ? 1.0 / std::sqrt(normL) : 1.0;
    double scaleR = normR > 0.0 ? 1.0 / std::sqrt(normR) : 1.0;
    for (int i = 0; i < N; ++i) {
        outWeightsL[i] *= scaleL;
        outWeightsR[i] *= scaleR;
    }
}

void t_kbeyond::mix_mid_side_to_lr(double tailMid, double tailSide, double widthNorm, double &outL, double &outR) const {
    widthNorm = clampd(widthNorm, 0.0, 2.0);
    double baseMid = 1.0;
    double baseSide = widthNorm;
    constexpr double eps = 1.0e-12;
    double absMid = std::abs(tailMid);
    double absSide = std::abs(tailSide);
    double ratio = absMid > eps ? absSide / absMid : 0.0;
    constexpr double leakComp = 4.0;
    double comp = ratio > 1.0 ? 1.0 / (1.0 + leakComp * (ratio - 1.0)) : 1.0;
    double midGain = baseMid * comp;
    double sideGain = baseSide;
    double norm = std::sqrt(midGain * midGain + sideGain * sideGain);
    if (norm <= 0.0) {
        outL = tailMid;
        outR = tailMid;
        return;
    }
    midGain /= norm;
    sideGain /= norm;
    outL = midGain * tailMid + sideGain * tailSide;
    outR = midGain * tailMid - sideGain * tailSide;
}

void t_kbeyond::update_output_weights() {
    auto midPattern = make_pattern(modeMid, N, 0x4001u);
    auto sidePattern = make_pattern(modeSide, N, 0x4002u);
    for (int i = 0; i < N; ++i) {
        double angle = (2.0 * M_PI * (double)i) / (double)N;
        double l = std::sin(angle * 0.91 + 0.17);
        double r = std::cos(angle * 1.07 - 0.11);
        double midScale = (i < (int)midPattern.size()) ? lerp(0.6, 1.4, midPattern[(std::size_t)i]) : 1.0;
        double sideScale = (i < (int)sidePattern.size()) ? lerp(0.6, 1.4, sidePattern[(std::size_t)i]) : 1.0;
        double mid = 0.5 * (l + r) * midScale;
        double side = 0.5 * (l - r) * sideScale;
        outBaseMid[i] = mid;
        outBaseSide[i] = side;
    }
    update_output_basis();
    apply_width(clampd(width, 0.0, 2.0));
}

void t_kbeyond::update_modulators() {
    kbeyond::dsp::update_modulators(modState, sr, modrate);
}

void t_kbeyond::update_decay() {
    kbeyond::dsp::update_decay(decayState,
                               sr,
                               regen,
                               decay,
                               dampLF,
                               dampMF,
                               dampHF,
                               fdnState.lengths);
    for (int i = 0; i < N; ++i) {
        double idx = (double)i / (double)std::max(1, N - 1);
        double lfCut = lerp(140.0, 320.0, idx);
        double hfCut = lerp(2800.0, sr * 0.4, idx);
        fdn_low[i].setCutoff(sr, lfCut);
        fdn_high[i].setCutoff(sr, hfCut);
        fdn_low[i].reset();
        fdn_high[i].reset();
    }
}

void t_kbeyond::update_quantum_walk() {
    kbeyond::dsp::update_quantum_walk(modState, sr, uwalkRate);
}

void t_kbeyond::reset_quantum_walk() {
    kbeyond::dsp::reset_quantum_walk(modState);
    update_quantum_walk();
}

void t_kbeyond::apply_quantum_dither(std::array<double, N> &vector) {
    kbeyond::dsp::apply_quantum_dither(modState, vector, coherence);
}

void t_kbeyond::apply_diffusion(const std::array<double, N> &input, std::array<double, N> &output) {
    const std::array<double, N> *source = &input;
    std::array<double, N> dithered {};
    if (clampd(coherence, 0.0, 1.0) > 0.0) {
        dithered = input;
        apply_quantum_dither(dithered);
        source = &dithered;
    }
    switch (mixMode) {
    case MixMode::Householder:
        apply_householder<N>(u, *source, output);
        break;
    case MixMode::WHT:
        kbeyond::dsp::mixing::apply_walsh_hadamard16(*source, output);
        break;
    case MixMode::Hybrid:
        kbeyond::dsp::mixing::apply_hybrid_diffusion(u, *source, output, diffusionScratch);
        break;
    default:
        apply_householder<N>(u, *source, output);
        break;
    }
}

void t_kbeyond::apply_quantum_walk(std::array<double, N> &feedback) {
    kbeyond::dsp::apply_quantum_walk(modState, feedback, coherence);
}

#ifndef KBEYOND_UNIT_TEST
// ------------------------------------------------------------ Attribute helpers

t_max_err kbeyond_attr_set_double(t_kbeyond* x, void*, long argc, t_atom* argv, double* target, double lo, double hi, void (t_kbeyond::*after)()) {
    if (!x || !target)
        return MAX_ERR_GENERIC;
    double v = atom_to_double(argc, argv, *target);
    v = clampd(v, lo, hi);
    *target = v;
    if (after)
        (x->*after)();
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_set_regen(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->regen, 0.0, 0.999, &t_kbeyond::update_decay);
}

t_max_err kbeyond_attr_set_decay(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->decay, 0.0, 20.0, &t_kbeyond::update_decay);
}

t_max_err kbeyond_attr_set_derez(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->derez, 0.0, 1.0, &t_kbeyond::refresh_filters);
}

t_max_err kbeyond_attr_set_filter(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->filter, 0.0, 1.0, &t_kbeyond::refresh_filters);
}

t_max_err kbeyond_attr_set_early(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->early, 0.0, 1.0, nullptr);
}

t_max_err kbeyond_attr_set_focus(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->focus, 0.0, 1.0, nullptr);
}

t_max_err kbeyond_attr_set_predelay(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->predelay, 0.0, t_kbeyond::kPredelayMaxSeconds, &t_kbeyond::setup_predelay);
    if (!err)
        x->predSamps = std::max(1.0, clampd(x->predelay * x->sr, 0.0, std::max(0.0, (double)x->predL.size() - 4.0)));
    return err;
}

t_max_err kbeyond_attr_set_mix(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->mix, 0.0, 1.0, nullptr);
}

t_max_err kbeyond_attr_set_width(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->width, 0.0, 2.0, &t_kbeyond::update_output_weights);
}

t_max_err kbeyond_attr_set_size(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->size, 0.0, 1.0, &t_kbeyond::setup_fdn);
    if (!err) {
        x->setup_early();
        x->update_modulators();
    }
    return err;
}

t_max_err kbeyond_attr_set_color(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->color, -1.0, 1.0, &t_kbeyond::refresh_filters);
}

t_max_err kbeyond_attr_set_damplf(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->dampLF, 0.0, 1.0, &t_kbeyond::update_decay);
}

t_max_err kbeyond_attr_set_dampmf(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->dampMF, 0.0, 1.0, &t_kbeyond::update_decay);
}

t_max_err kbeyond_attr_set_damphf(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->dampHF, 0.0, 1.0, &t_kbeyond::update_decay);
}

t_max_err kbeyond_attr_set_modrate(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->modrate, 0.0, 5.0, &t_kbeyond::update_modulators);
}

t_max_err kbeyond_attr_set_moddepth(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->moddepth, 0.0, 32.0, &t_kbeyond::setup_fdn);
    if (!err)
        x->update_modulators();
    return err;
}

t_max_err kbeyond_attr_set_motion(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->motion, 0.0, 1.0, nullptr);
}

t_max_err kbeyond_attr_set_phiweight(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->phiweight, 0.0, 1.0, &t_kbeyond::update_diffusion);
}

t_max_err kbeyond_attr_set_laser(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->laser, 0.0, 1.0, nullptr);
}

t_max_err kbeyond_attr_set_laserfocus(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->laserFocus, 0.0, 1.0, &t_kbeyond::setup_early);
    if (!err)
        x->update_laser_phase_inc();
    return err;
}

t_max_err kbeyond_attr_set_lasergate(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->laserGate, 0.0, 1.0, nullptr);
    if (!err)
        x->update_laser_gate();
    return err;
}

t_max_err kbeyond_attr_set_laserwindow(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->laserWindow, 0.2, 0.6, nullptr);
    if (!err)
        x->update_laser_window();
    return err;
}

t_max_err kbeyond_attr_set_laserdiffusion(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->laserDiffusion, 0.0, 1.0, nullptr);
}

t_max_err kbeyond_attr_set_coherence(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->coherence, 0.0, 1.0, &t_kbeyond::update_quantum_walk);
}

t_max_err kbeyond_attr_set_uwalkrate(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->uwalkRate, 0.0, 8.0, &t_kbeyond::update_quantum_walk);
}

static t_symbol* symbol_from_pattern(prime_modes::Pattern pattern) {
    switch (pattern) {
    case prime_modes::Pattern::Prime:
        return gensym("prime");
    case prime_modes::Pattern::Aureo:
        return gensym("aureo");
    case prime_modes::Pattern::Plastica:
        return gensym("plastica");
    case prime_modes::Pattern::PrimeAureo:
        return gensym("prime_aureo");
    default:
        return gensym("aureo");
    }
}

static prime_modes::Pattern pattern_from_symbol(t_symbol* sym, prime_modes::Pattern fallback) {
    if (!sym)
        return fallback;
    if (sym == gensym("prime"))
        return prime_modes::Pattern::Prime;
    if (sym == gensym("aureo"))
        return prime_modes::Pattern::Aureo;
    if (sym == gensym("plastica"))
        return prime_modes::Pattern::Plastica;
    if (sym == gensym("prime_aureo"))
        return prime_modes::Pattern::PrimeAureo;
    return fallback;
}

t_max_err kbeyond_attr_set_mode_er(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    if (!x)
        return MAX_ERR_GENERIC;
    (void)attr;
    t_symbol* sym = (argc > 0 && argv) ? atom_getsym(argv) : x->modeERSym;
    if (!sym)
        sym = gensym("aureo");
    x->modeER = pattern_from_symbol(sym, x->modeER);
    x->modeERSym = symbol_from_pattern(x->modeER);
    x->setup_early();
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_get_mode_er(t_kbeyond* x, void* attr, long* argc, t_atom** argv) {
    if (!x || !argc || !argv)
        return MAX_ERR_GENERIC;
    (void)attr;
    if (!*argv)
        *argv = (t_atom*)sysmem_newptr(sizeof(t_atom));
    if (!*argv) {
        *argc = 0;
        return MAX_ERR_GENERIC;
    }
    *argc = 1;
    t_symbol* sym = x->modeERSym ? x->modeERSym : symbol_from_pattern(x->modeER);
    atom_setsym(*argv, sym);
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_set_mode_late(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    if (!x)
        return MAX_ERR_GENERIC;
    (void)attr;
    t_symbol* sym = (argc > 0 && argv) ? atom_getsym(argv) : x->modeLateSym;
    if (!sym)
        sym = gensym("prime_aureo");
    x->modeLate = pattern_from_symbol(sym, x->modeLate);
    x->modeLateSym = symbol_from_pattern(x->modeLate);
    x->setup_fdn();
    x->update_output_weights();
    x->update_modulators();
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_get_mode_late(t_kbeyond* x, void* attr, long* argc, t_atom** argv) {
    if (!x || !argc || !argv)
        return MAX_ERR_GENERIC;
    (void)attr;
    if (!*argv)
        *argv = (t_atom*)sysmem_newptr(sizeof(t_atom));
    if (!*argv) {
        *argc = 0;
        return MAX_ERR_GENERIC;
    }
    *argc = 1;
    t_symbol* sym = x->modeLateSym ? x->modeLateSym : symbol_from_pattern(x->modeLate);
    atom_setsym(*argv, sym);
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_set_mode_mid(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    if (!x)
        return MAX_ERR_GENERIC;
    (void)attr;
    t_symbol* sym = (argc > 0 && argv) ? atom_getsym(argv) : x->modeMidSym;
    if (!sym)
        sym = gensym("prime");
    x->modeMid = pattern_from_symbol(sym, x->modeMid);
    x->modeMidSym = symbol_from_pattern(x->modeMid);
    x->update_injection_weights();
    x->setup_early();
    x->update_output_weights();
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_get_mode_mid(t_kbeyond* x, void* attr, long* argc, t_atom** argv) {
    if (!x || !argc || !argv)
        return MAX_ERR_GENERIC;
    (void)attr;
    if (!*argv)
        *argv = (t_atom*)sysmem_newptr(sizeof(t_atom));
    if (!*argv) {
        *argc = 0;
        return MAX_ERR_GENERIC;
    }
    *argc = 1;
    t_symbol* sym = x->modeMidSym ? x->modeMidSym : symbol_from_pattern(x->modeMid);
    atom_setsym(*argv, sym);
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_set_mode_side(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    if (!x)
        return MAX_ERR_GENERIC;
    (void)attr;
    t_symbol* sym = (argc > 0 && argv) ? atom_getsym(argv) : x->modeSideSym;
    if (!sym)
        sym = gensym("aureo");
    x->modeSide = pattern_from_symbol(sym, x->modeSide);
    x->modeSideSym = symbol_from_pattern(x->modeSide);
    x->setup_early();
    x->update_output_weights();
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_get_mode_side(t_kbeyond* x, void* attr, long* argc, t_atom** argv) {
    if (!x || !argc || !argv)
        return MAX_ERR_GENERIC;
    (void)attr;
    if (!*argv)
        *argv = (t_atom*)sysmem_newptr(sizeof(t_atom));
    if (!*argv) {
        *argc = 0;
        return MAX_ERR_GENERIC;
    }
    *argc = 1;
    t_symbol* sym = x->modeSideSym ? x->modeSideSym : symbol_from_pattern(x->modeSide);
    atom_setsym(*argv, sym);
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_set_seed(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    if (!x)
        return MAX_ERR_GENERIC;
    (void)attr;
    long value = x->patternSeed;
    if (argc > 0 && argv) {
        if (atom_gettype(argv) == A_LONG)
            value = atom_getlong(argv);
        else if (atom_gettype(argv) == A_FLOAT)
            value = (long)std::llround(atom_getfloat(argv));
    }
    if (value < 0)
        value = -value;
    x->patternSeed = value;
    x->setup_early();
    x->setup_fdn();
    x->update_output_weights();
    x->update_modulators();
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_get_seed(t_kbeyond* x, void* attr, long* argc, t_atom** argv) {
    if (!x || !argc || !argv)
        return MAX_ERR_GENERIC;
    (void)attr;
    if (!*argv)
        *argv = (t_atom*)sysmem_newptr(sizeof(t_atom));
    if (!*argv) {
        *argc = 0;
        return MAX_ERR_GENERIC;
    }
    *argc = 1;
    atom_setlong(*argv, x->patternSeed);
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_set_mode_mix(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    if (!x)
        return MAX_ERR_GENERIC;
    (void)attr;
    t_symbol* sym = nullptr;
    if (argc > 0 && argv)
        sym = atom_getsym(argv);
    if (!sym)
        sym = gensym("householder");

    t_symbol* symHouse = gensym("householder");
    t_symbol* symWht = gensym("wht");
    t_symbol* symHybrid = gensym("hybrid");

    t_kbeyond::MixMode newMode = t_kbeyond::MixMode::Householder;
    if (sym == symHouse) {
        newMode = t_kbeyond::MixMode::Householder;
    } else if (sym == symWht) {
        newMode = t_kbeyond::MixMode::WHT;
    } else if (sym == symHybrid) {
        newMode = t_kbeyond::MixMode::Hybrid;
    } else {
        sym = symHouse;
    }

    x->mixMode = newMode;
    x->modeMixSym = sym;
    x->update_diffusion();
    x->update_output_weights();
    return MAX_ERR_NONE;
}

t_max_err kbeyond_attr_get_mode_mix(t_kbeyond* x, void* attr, long* argc, t_atom** argv) {
    if (!x || !argc || !argv)
        return MAX_ERR_GENERIC;
    (void)attr;
    if (!*argv)
        *argv = (t_atom*)sysmem_newptr(sizeof(t_atom));
    if (!*argv) {
        *argc = 0;
        return MAX_ERR_GENERIC;
    }
    *argc = 1;
    t_symbol* sym = x->modeMixSym ? x->modeMixSym : gensym("householder");
    atom_setsym(*argv, sym);
    return MAX_ERR_NONE;
}

// ------------------------------------------------------------ Object lifecycle

void *kbeyond_new(t_symbol *, long argc, t_atom *argv) {
    t_kbeyond *x = (t_kbeyond *)object_alloc(s_kbeyond_class);
    if (!x)
        return nullptr;

    dsp_setup((t_pxobject *)x, 2);
    outlet_new((t_object *)x, "signal");
    outlet_new((t_object *)x, "signal");

    x->mixMode = t_kbeyond::MixMode::Householder;
    x->modeMixSym = gensym("householder");

    x->modeER = prime_modes::Pattern::Aureo;
    x->modeLate = prime_modes::Pattern::PrimeAureo;
    x->modeMid = prime_modes::Pattern::Prime;
    x->modeSide = prime_modes::Pattern::Aureo;
    x->modeERSym = symbol_from_pattern(x->modeER);
    x->modeLateSym = symbol_from_pattern(x->modeLate);
    x->modeMidSym = symbol_from_pattern(x->modeMid);
    x->modeSideSym = symbol_from_pattern(x->modeSide);
    x->patternSeed = 1337;

    x->coherence = 0.8;
    x->uwalkRate = 0.25;
    x->reset_quantum_walk();

    if (argc > 0 && (atom_gettype(argv) == A_LONG || atom_gettype(argv) == A_FLOAT))
        kbeyond_attr_set_phiweight(x, nullptr, 1, argv);

    attr_args_process(x, argc, argv);

    x->setup_sr(sys_getsr() > 0.0 ? sys_getsr() : 48000.0);

    return x;
}

void kbeyond_free(t_kbeyond *x) {
    dsp_free((t_pxobject *)x);
}

void kbeyond_assist(t_kbeyond *x, void *b, long m, long a, char *s) {
    if (m == ASSIST_INLET) {
        if (a == 0)
            std::snprintf(s, kAssistStringMax, "Signal L (with attributes)");
        else
            std::snprintf(s, kAssistStringMax, "Signal R");
    } else {
        if (a == 0)
            std::snprintf(s, kAssistStringMax, "Wet L out");
        else
            std::snprintf(s, kAssistStringMax, "Wet R out");
    }
}

void kbeyond_dsp64(t_kbeyond *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long) {
    x->vs = maxvectorsize;
    x->setup_sr(samplerate);
    object_method(dsp64, gensym("dsp_add64"), x, (t_perfroutine64)kbeyond_perform64, 0, nullptr);
    (void)count;
}

#endif // KBEYOND_UNIT_TEST

void kbeyond_perform64(t_kbeyond *x, t_object *, double **ins, long nin, double **outs, long nout, long sampleframes, long, void *) {
    if (!x || x->ob.z_disabled) {
        for (long ch = 0; ch < nout; ++ch) {
            if (outs[ch])
                std::fill(outs[ch], outs[ch] + sampleframes, 0.0);
        }
        return;
    }

    double *outL = (nout > 0 && outs[0]) ? outs[0] : nullptr;
    double *outR = (nout > 1 && outs[1]) ? outs[1] : nullptr;

    double motionAmt = clampd(x->motion, 0.0, 1.0);
    double rangeEnv = clampd(x->rangeEnv, 0.0, 1.0);
    double dopplerEnv = clampd(x->dopplerEnv, 0.0, 1.0);
    double spreadEnv = clampd(x->spreadEnv, 0.0, 1.0);

    double baseWidth = clampd(x->width, 0.0, 2.0);
    double rangeLift = std::max(rangeEnv - 0.35, 0.0);
    double widthDelta = motionAmt * (0.45 * rangeLift + 0.35 * (spreadEnv - 0.5) * rangeEnv);
    double widthNorm = clampd(baseWidth + widthDelta, 0.0, 2.0);
    x->apply_width(widthNorm);

    double mix = clampd(x->mix, 0.0, 1.0);
    double earlyAmt = clampd(x->early, 0.0, 1.0);
    double focusAmt = clampd(x->focus, 0.0, 1.0);
    double dryGainBase = std::cos(0.5 * M_PI * mix);
    double wetGainBase = std::sin(0.5 * M_PI * mix);

    double baseModDepth = clampd(x->moddepth, 0.0, 32.0);
    double moddepth = clampd(baseModDepth * (1.0 + motionAmt * 0.65 * dopplerEnv), 0.0, 32.0);

    double dampBoost = std::max(rangeEnv - 0.3, 0.0);
    double dampMFNormBase = clampd(x->dampMF, 0.0, 1.0);
    double dampHFNormBase = clampd(x->dampHF, 0.0, 1.0);
    double dampMFNorm = clampd(dampMFNormBase + motionAmt * 0.55 * dampBoost, 0.0, 1.0);
    double dampHFNorm = clampd(dampHFNormBase + motionAmt * (0.45 * dampBoost + 0.55 * dopplerEnv), 0.0, 1.0);
    double dampLFValue = x->decayState.dampLF;
    double dampMFBase = x->decayState.dampMF;
    double dampHFBase = x->decayState.dampHF;
    double dampMFValue = std::max(lerp(1.0, 0.45, dampMFNorm), dampMFBase);
    double dampHFValue = std::max(lerp(1.0, 0.12, dampHFNorm), dampHFBase);

#ifdef KBEYOND_UNIT_TEST
    x->debugWidthTarget = widthNorm;
    x->debugModDepthTarget = moddepth;
    x->debugDampMFValue = dampMFValue;
    x->debugDampHFValue = dampHFValue;
#endif

    for (long i = 0; i < sampleframes; ++i) {
        double inL = (nin > 0 && ins[0]) ? ins[0][i] : 0.0;
        double inR = (nin > 1 && ins[1]) ? ins[1][i] : inL;

        double dryL = inL;
        double dryR = inR;

        double predLen = x->predSamps;
        double predOutL = x->predL.readFrac(predLen);
        double predOutR = x->predR.readFrac(predLen);
        x->predL.write(inL + kbeyond::dsp::tiny_noise(x->rng));
        x->predR.write(inR + kbeyond::dsp::tiny_noise(x->rng));

        double midIn = 0.5 * (predOutL + predOutR);
        double sideIn = 0.5 * (predOutL - predOutR);

        x->motionDetector.process(midIn, sideIn);

        double earlyL = 0.0;
        double earlyR = 0.0;
#if KBEYOND_ENABLE_LASER
        double clusterAmt = clampd(x->laser, 0.0, 1.0);
#else
        double clusterAmt = 0.0;
#endif
        x->earlySection.render(predOutL,
                               predOutR,
                               widthNorm,
                               earlyAmt,
                               focusAmt,
                               clusterAmt,
                               earlyL,
                               earlyR,
                               x->rng);

        std::array<double, t_kbeyond::N> vec {};
        kbeyond::dsp::read_lines(x->fdnState,
                                 moddepth,
                                 x->modState.fdnPhase,
                                 x->modState.fdnPhaseInc,
                                 x->fdn_tilt,
                                 x->fdn_lp,
                                 vec);

        x->apply_diffusion(vec, x->fdnState.feedback);
#if KBEYOND_ENABLE_LASER
        double qMix = x->earlySection.computeQSwitchMix(clusterAmt, x->laserDiffusion);
        if (qMix > 0.0) {
            std::array<double, t_kbeyond::N> alt {};
            kbeyond::dsp::mixing::apply_walsh_hadamard16(vec, alt);
            for (int l = 0; l < t_kbeyond::N; ++l)
                x->fdnState.feedback[l] = lerp(x->fdnState.feedback[l], alt[l], qMix);
        }
#else
        (void)clusterAmt;
#endif
        x->apply_quantum_walk(x->fdnState.feedback);

        double tailMid = 0.0;
        double tailSide = 0.0;
        kbeyond::dsp::write_feedback(x->fdnState,
                                     midIn,
                                     sideIn,
                                     x->rng,
                                     x->inWeights,
                                     x->outWeightsSide,
                                     x->outMidBasis,
                                     x->decayState.perLine,
                                     x->fdn_low,
                                     x->fdn_high,
                                     dampLFValue,
                                     dampMFValue,
                                     dampHFValue,
                                     tailMid,
                                     tailSide);

        double tailL = 0.0;
        double tailR = 0.0;
        x->mix_mid_side_to_lr(tailMid, tailSide, widthNorm, tailL, tailR);

        double wetL = earlyL + tailL;
        double wetR = earlyR + tailR;

        double dryEnergy = 0.5 * (dryL * dryL + dryR * dryR);
        double wetEnergy = 0.5 * (wetL * wetL + wetR * wetR);
        constexpr double envCoef = 0.995;
        x->dryEnergyEnv = lerp(dryEnergy, x->dryEnergyEnv, envCoef);
        x->wetEnergyEnv = lerp(wetEnergy, x->wetEnergyEnv, envCoef);
        x->dryLevelEnv = std::sqrt(std::max(x->dryEnergyEnv, 0.0));
        x->wetLevelEnv = std::sqrt(std::max(x->wetEnergyEnv, 0.0));

        double compTarget = 1.0;
        if (mix > 0.0) {
            constexpr double eps = 1.0e-12;
            constexpr double makeupMin = 1.0;
            constexpr double makeupMax = 8.0;
            double ratio = std::sqrt((x->dryEnergyEnv + eps) / (x->wetEnergyEnv + eps));
            // Never allow the compensation to dip below unity when the wet tail dominates.
            ratio = std::max(ratio, 1.0);
            // Clamp to a safe yet generous range so quiet tails are lifted without runaway gain.
            compTarget = clampd(ratio, makeupMin, makeupMax);
        }
        constexpr double makeupCoef = 0.9975;
        x->wetMakeup = lerp(compTarget, x->wetMakeup, makeupCoef);

        double wetGain = wetGainBase;
        if (mix > 0.0) {
            double makeup = (mix >= 1.0) ? x->wetMakeup : (1.0 + mix * (x->wetMakeup - 1.0));
            wetGain *= makeup;
        }

        double dryOutL = dryGainBase * dryL;
        double dryOutR = dryGainBase * dryR;
        double wetOutL = wetGain * wetL;
        double wetOutR = wetGain * wetR;

        if (outL)
            outL[i] = dryOutL + wetOutL;
        if (outR)
            outR[i] = dryOutR + wetOutR;
    }

    x->rangeEnv = clampd(x->motionDetector.range(), 0.0, 1.0);
    x->dopplerEnv = clampd(x->motionDetector.doppler(), 0.0, 1.0);
    x->spreadEnv = clampd(x->motionDetector.spread(), 0.0, 1.0);
}

#ifndef KBEYOND_UNIT_TEST
extern "C" C74_EXPORT void ext_main(void *r) {
    t_class *c = class_new("kbeyond~", (method)kbeyond_new, (method)kbeyond_free, sizeof(t_kbeyond), 0L, A_GIMME, 0);

    class_addmethod(c, (method)kbeyond_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)kbeyond_assist, "assist", A_CANT, 0);

    CLASS_ATTR_DOUBLE(c, "regen", 0, t_kbeyond, regen);
    CLASS_ATTR_ACCESSORS(c, "regen", NULL, kbeyond_attr_set_regen);
    CLASS_ATTR_FILTER_CLIP(c, "regen", 0.0, 0.999);

    CLASS_ATTR_DOUBLE(c, "decay", 0, t_kbeyond, decay);
    CLASS_ATTR_ACCESSORS(c, "decay", NULL, kbeyond_attr_set_decay);
    CLASS_ATTR_FILTER_CLIP(c, "decay", 0.0, 20.0);

    CLASS_ATTR_DOUBLE(c, "derez", 0, t_kbeyond, derez);
    CLASS_ATTR_ACCESSORS(c, "derez", NULL, kbeyond_attr_set_derez);
    CLASS_ATTR_FILTER_CLIP(c, "derez", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "filter", 0, t_kbeyond, filter);
    CLASS_ATTR_ACCESSORS(c, "filter", NULL, kbeyond_attr_set_filter);
    CLASS_ATTR_FILTER_CLIP(c, "filter", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "early", 0, t_kbeyond, early);
    CLASS_ATTR_ACCESSORS(c, "early", NULL, kbeyond_attr_set_early);
    CLASS_ATTR_FILTER_CLIP(c, "early", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "focus", 0, t_kbeyond, focus);
    CLASS_ATTR_ACCESSORS(c, "focus", NULL, kbeyond_attr_set_focus);
    CLASS_ATTR_FILTER_CLIP(c, "focus", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "predelay", 0, t_kbeyond, predelay);
    CLASS_ATTR_ACCESSORS(c, "predelay", NULL, kbeyond_attr_set_predelay);
    CLASS_ATTR_FILTER_CLIP(c, "predelay", 0.0, t_kbeyond::kPredelayMaxSeconds);

    CLASS_ATTR_DOUBLE(c, "mix", 0, t_kbeyond, mix);
    CLASS_ATTR_ACCESSORS(c, "mix", NULL, kbeyond_attr_set_mix);
    CLASS_ATTR_FILTER_CLIP(c, "mix", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "width", 0, t_kbeyond, width);
    CLASS_ATTR_ACCESSORS(c, "width", NULL, kbeyond_attr_set_width);
    CLASS_ATTR_FILTER_CLIP(c, "width", 0.0, 2.0);

    CLASS_ATTR_DOUBLE(c, "size", 0, t_kbeyond, size);
    CLASS_ATTR_ACCESSORS(c, "size", NULL, kbeyond_attr_set_size);
    CLASS_ATTR_FILTER_CLIP(c, "size", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "color", 0, t_kbeyond, color);
    CLASS_ATTR_ACCESSORS(c, "color", NULL, kbeyond_attr_set_color);
    CLASS_ATTR_FILTER_CLIP(c, "color", -1.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "damplf", 0, t_kbeyond, dampLF);
    CLASS_ATTR_ACCESSORS(c, "damplf", NULL, kbeyond_attr_set_damplf);
    CLASS_ATTR_FILTER_CLIP(c, "damplf", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "dampmf", 0, t_kbeyond, dampMF);
    CLASS_ATTR_ACCESSORS(c, "dampmf", NULL, kbeyond_attr_set_dampmf);
    CLASS_ATTR_FILTER_CLIP(c, "dampmf", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "damphf", 0, t_kbeyond, dampHF);
    CLASS_ATTR_ACCESSORS(c, "damphf", NULL, kbeyond_attr_set_damphf);
    CLASS_ATTR_FILTER_CLIP(c, "damphf", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "modrate", 0, t_kbeyond, modrate);
    CLASS_ATTR_ACCESSORS(c, "modrate", NULL, kbeyond_attr_set_modrate);
    CLASS_ATTR_FILTER_CLIP(c, "modrate", 0.0, 5.0);

    CLASS_ATTR_DOUBLE(c, "moddepth", 0, t_kbeyond, moddepth);
    CLASS_ATTR_ACCESSORS(c, "moddepth", NULL, kbeyond_attr_set_moddepth);
    CLASS_ATTR_FILTER_CLIP(c, "moddepth", 0.0, 32.0);

    CLASS_ATTR_DOUBLE(c, "motion", 0, t_kbeyond, motion);
    CLASS_ATTR_ACCESSORS(c, "motion", NULL, kbeyond_attr_set_motion);
    CLASS_ATTR_FILTER_CLIP(c, "motion", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "phiweight", 0, t_kbeyond, phiweight);
    CLASS_ATTR_ACCESSORS(c, "phiweight", NULL, kbeyond_attr_set_phiweight);
    CLASS_ATTR_FILTER_CLIP(c, "phiweight", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "laser", 0, t_kbeyond, laser);
    CLASS_ATTR_ACCESSORS(c, "laser", NULL, kbeyond_attr_set_laser);
    CLASS_ATTR_FILTER_CLIP(c, "laser", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "laserfocus", 0, t_kbeyond, laserFocus);
    CLASS_ATTR_ACCESSORS(c, "laserfocus", NULL, kbeyond_attr_set_laserfocus);
    CLASS_ATTR_FILTER_CLIP(c, "laserfocus", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "lasergate", 0, t_kbeyond, laserGate);
    CLASS_ATTR_ACCESSORS(c, "lasergate", NULL, kbeyond_attr_set_lasergate);
    CLASS_ATTR_FILTER_CLIP(c, "lasergate", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "laserwindow", 0, t_kbeyond, laserWindow);
    CLASS_ATTR_ACCESSORS(c, "laserwindow", NULL, kbeyond_attr_set_laserwindow);
    CLASS_ATTR_FILTER_CLIP(c, "laserwindow", 0.2, 0.6);

    CLASS_ATTR_DOUBLE(c, "laserdiffusion", 0, t_kbeyond, laserDiffusion);
    CLASS_ATTR_ACCESSORS(c, "laserdiffusion", NULL, kbeyond_attr_set_laserdiffusion);
    CLASS_ATTR_FILTER_CLIP(c, "laserdiffusion", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "coherence", 0, t_kbeyond, coherence);
    CLASS_ATTR_ACCESSORS(c, "coherence", NULL, kbeyond_attr_set_coherence);
    CLASS_ATTR_FILTER_CLIP(c, "coherence", 0.0, 1.0);

    CLASS_ATTR_DOUBLE(c, "uwalkrate", 0, t_kbeyond, uwalkRate);
    CLASS_ATTR_ACCESSORS(c, "uwalkrate", NULL, kbeyond_attr_set_uwalkrate);
    CLASS_ATTR_FILTER_CLIP(c, "uwalkrate", 0.0, 8.0);

    CLASS_ATTR_SYM(c, "mode_er", 0, t_kbeyond, modeERSym);
    CLASS_ATTR_ACCESSORS(c, "mode_er", kbeyond_attr_get_mode_er, kbeyond_attr_set_mode_er);
    CLASS_ATTR_ENUM(c, "mode_er", 0, "prime aureo plastica prime_aureo");
    CLASS_ATTR_LABEL(c, "mode_er", 0, "Early Pattern Mode");

    CLASS_ATTR_SYM(c, "mode_late", 0, t_kbeyond, modeLateSym);
    CLASS_ATTR_ACCESSORS(c, "mode_late", kbeyond_attr_get_mode_late, kbeyond_attr_set_mode_late);
    CLASS_ATTR_ENUM(c, "mode_late", 0, "prime aureo plastica prime_aureo");
    CLASS_ATTR_LABEL(c, "mode_late", 0, "Late Pattern Mode");

    CLASS_ATTR_SYM(c, "mode_mid", 0, t_kbeyond, modeMidSym);
    CLASS_ATTR_ACCESSORS(c, "mode_mid", kbeyond_attr_get_mode_mid, kbeyond_attr_set_mode_mid);
    CLASS_ATTR_ENUM(c, "mode_mid", 0, "prime aureo plastica prime_aureo");
    CLASS_ATTR_LABEL(c, "mode_mid", 0, "Mid Routing Mode");

    CLASS_ATTR_SYM(c, "mode_side", 0, t_kbeyond, modeSideSym);
    CLASS_ATTR_ACCESSORS(c, "mode_side", kbeyond_attr_get_mode_side, kbeyond_attr_set_mode_side);
    CLASS_ATTR_ENUM(c, "mode_side", 0, "prime aureo plastica prime_aureo");
    CLASS_ATTR_LABEL(c, "mode_side", 0, "Side Routing Mode");

    CLASS_ATTR_LONG(c, "seed", 0, t_kbeyond, patternSeed);
    CLASS_ATTR_ACCESSORS(c, "seed", kbeyond_attr_get_seed, kbeyond_attr_set_seed);
    CLASS_ATTR_FILTER_MIN(c, "seed", 0);
    CLASS_ATTR_LABEL(c, "seed", 0, "Pattern Seed");

    CLASS_ATTR_SYM(c, "mode_mix", 0, t_kbeyond, modeMixSym);
    CLASS_ATTR_ACCESSORS(c, "mode_mix", kbeyond_attr_get_mode_mix, kbeyond_attr_set_mode_mix);
    CLASS_ATTR_ENUM(c, "mode_mix", 0, "householder wht hybrid");
    CLASS_ATTR_LABEL(c, "mode_mix", 0, "Diffusion Mode");

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    s_kbeyond_class = c;
}

#endif // KBEYOND_UNIT_TEST

