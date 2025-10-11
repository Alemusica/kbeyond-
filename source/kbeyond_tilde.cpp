#include "kbeyond_tilde.h"
#include "householder_phi16.h"

#include <algorithm>
#include <numeric>
#include <cstdio>
#include <cstddef>
#include <cmath>

static constexpr std::size_t kAssistStringMax = 256;

static t_class* s_kbeyond_class = nullptr;

static void apply_walsh_hadamard16(const std::array<double, t_kbeyond::N> &input,
                                   std::array<double, t_kbeyond::N> &output) {
    std::array<double, t_kbeyond::N> tmp = input;
    for (std::size_t len = 1; len < tmp.size(); len <<= 1) {
        std::size_t step = len << 1;
        for (std::size_t i = 0; i < tmp.size(); i += step) {
            for (std::size_t j = 0; j < len; ++j) {
                double a = tmp[i + j];
                double b = tmp[i + j + len];
                tmp[i + j] = a + b;
                tmp[i + j + len] = a - b;
            }
        }
    }
    double norm = 1.0 / std::sqrt(static_cast<double>(tmp.size()));
    for (std::size_t i = 0; i < tmp.size(); ++i)
        output[i] = tmp[i] * norm;
}

static void apply_hybrid_diffusion(const std::array<double, t_kbeyond::N> &u,
                                   const std::array<double, t_kbeyond::N> &input,
                                   std::array<double, t_kbeyond::N> &output,
                                   std::array<double, t_kbeyond::N> &scratch) {
    apply_walsh_hadamard16(input, scratch);
    apply_householder<t_kbeyond::N>(u, scratch, output);
}

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

void t_kbeyond::setup_sr(double newsr) {
    sr = newsr > 1.0 ? newsr : 48000.0;
    vs = std::max<long>(64, vs);
    setup_predelay();
    setup_early();
    update_laser_gate();
    update_laser_window();
    update_laser_phase_inc();
    update_laser_envelope();
    laserEnv = 0.0;
    laserExcite = 0.0;
    qswitchEnv = 0.0;
    qswitchCounter = 0;
    laserPhase = 0.0;
    setup_fdn();
    refresh_filters();
    update_diffusion();
    update_output_weights();
    update_modulators();
    update_quantum_walk();
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
}

void t_kbeyond::setup_early() {
    double maxMs = 120.0;
    size_t len = (size_t)std::ceil(ms2samp(maxMs * 1.25, sr) + 8.0);
    earlyBufMid.setup(len);
    earlyBufSide.setup(len);
    const double phi = 1.6180339887498948482;
    double baseMs = 5.2;
    double scale = lerp(0.7, 1.45, size);
    int pairs = kEarlyTaps / 2;
    int patternCount = pairs;
    if (kEarlyTaps % 2 != 0)
        ++patternCount;
    if (patternCount <= 0)
        patternCount = 1;
    auto tapPattern = make_pattern(modeER, (std::size_t)patternCount, 0x1001u);
    auto laserMidPattern = make_pattern(modeMid, kLaserTaps, 0x5001u);
    auto laserSidePattern = make_pattern(modeSide, kLaserTaps, 0x5002u);
    auto storeTap = [&](int index, long samp, double gain, double pan) {
        double theta = (pan + 1.0) * (0.25 * M_PI);
        earlyDel[index] = samp;
        earlyGain[index] = gain;
        earlyCos[index] = std::cos(theta);
        earlySin[index] = std::sin(theta);
    };
    for (int p = 0; p < pairs; ++p) {
        double idxNorm = (p < (int)tapPattern.size())
                             ? tapPattern[(std::size_t)p]
                             : ((pairs > 1) ? (double)p / (double)(pairs - 1) : 0.5);
        double ms = baseMs * std::pow(phi, idxNorm * 1.2);
        ms = std::min(ms * scale, maxMs);
        long samp = (long)clampd(std::floor(ms2samp(ms, sr)), 1.0, (double)earlyBufMid.size() - 2.0);
        double gain = std::pow(0.72, (double)p + 1.0);
        double panMag = 1.0 - idxNorm;
        double panLeft = -panMag;
        double panRight = panMag;
        storeTap(p, samp, gain, panLeft);
        storeTap(kEarlyTaps - 1 - p, samp, gain, panRight);
    }
    if (kEarlyTaps % 2 != 0) {
        int center = pairs;
        double idxNorm = (center < (int)tapPattern.size()) ? tapPattern[(std::size_t)center] : 0.65;
        double ms = baseMs * std::pow(phi, idxNorm);
        ms = std::min(ms * scale, maxMs);
        long samp = (long)clampd(std::floor(ms2samp(ms, sr)), 1.0, (double)earlyBufMid.size() - 2.0);
        double gain = std::pow(0.72, (double)pairs + 1.0);
        storeTap(center, samp, gain, 0.0);
    }

    double focus = clampd(laserFocus, 0.0, 1.0);
    double baseLaser = lerp(6.5, 17.5, size);
    double groupSpacing = lerp(3.0, 9.5, size);
    double clusterSpan = lerp(2.2, 7.8, focus);
    double chirpSpan = lerp(0.6, 4.2, focus);
    double jitterAmt = lerp(0.08, 0.35, focus);
    uint32_t lfsr = 0x1Fu;
    auto next_code = [&]() {
        int bit = ((lfsr >> 4) ^ (lfsr >> 2)) & 1;
        lfsr = ((lfsr << 1) & 0x1F) | bit;
        return (lfsr & 1) ? 1.0 : -1.0;
    };
    double maxDelay = (double)earlyBufMid.size() - 2.0;
    for (int idx = 0; idx < kLaserTaps; ++idx) {
        double norm = (kLaserTaps > 1) ? (double)idx / (double)(kLaserTaps - 1) : 0.0;
        double groupFloat = norm * (double)kLaserGroups;
        int group = (int)std::floor(groupFloat);
        if (group >= kLaserGroups)
            group = kLaserGroups - 1;
        double local = groupFloat - (double)group;
        double startMs = baseLaser + groupSpacing * (double)group;
        double sweep = clusterSpan * (1.0 + 0.27 * (double)group);
        double chirp = chirpSpan * local * local;
        double jitter = jitterAmt * next_code() * (1.0 - local);
        double tapMs = clampd(startMs + local * sweep + chirp + jitter, 1.0, maxMs - 1.0);
        double delaySamp = clampd(ms2samp(tapMs, sr), 1.0, maxDelay);
        laserDelay[idx] = delaySamp;

        double shape = next_code();
        laserShape[idx] = shape;
        double baseGain = 0.22 * std::pow(0.78, (double)idx * 0.55);
        double bright = lerp(0.78, 1.32, focus) * (1.0 + 0.22 * shape);
        double midScale = (idx < (int)laserMidPattern.size())
                              ? lerp(0.65, 1.45, laserMidPattern[(std::size_t)idx])
                              : 1.0;
        laserMidGain[idx] = baseGain * bright * midScale;
        double sideSign = next_code();
        double sideScale = (idx < (int)laserSidePattern.size())
                               ? lerp(0.65, 1.45, laserSidePattern[(std::size_t)idx])
                               : 1.0;
        double sideWeight = lerp(0.18, 0.55, focus) * sideSign * sideScale;
        laserSideGain[idx] = baseGain * sideWeight;
        double panBase = lerp(-0.9, 0.9, norm);
        double panJitter = 0.25 * next_code() * (1.0 - local);
        double pan = clampd(panBase + panJitter, -1.0, 1.0);
        double theta = (pan + 1.0) * (0.25 * M_PI);
        laserCos[idx] = std::cos(theta);
        laserSin[idx] = std::sin(theta);
    }
}

void t_kbeyond::update_laser_gate() {
    double g = clampd(laserGate, 0.0, 1.0);
    laserGateScaled = 0.004 + 0.25 * g * g;
}

void t_kbeyond::update_laser_window() {
    double winSec = clampd(laserWindow, 0.2, 0.6);
    qswitchWindowSamples = (long)std::llround(winSec * sr);
    if (qswitchWindowSamples < 1)
        qswitchWindowSamples = 1;
    if (qswitchCounter > qswitchWindowSamples)
        qswitchCounter = qswitchWindowSamples;
}

void t_kbeyond::update_laser_phase_inc() {
    double focus = clampd(laserFocus, 0.0, 1.0);
    double sweepHz = lerp(0.5, 7.0, focus);
    laserPhaseInc = 2.0 * M_PI * sweepHz / std::max(1.0, sr);
}

void t_kbeyond::update_laser_envelope() {
    double attackTime = 0.0035;
    double releaseTime = 0.12;
    laserEnvAttack = std::exp(-1.0 / std::max(1.0, sr * attackTime));
    laserEnvRelease = std::exp(-1.0 / std::max(1.0, sr * releaseTime));
    double qAttackTime = 0.004;
    double qReleaseTime = 0.24;
    qswitchAttack = std::exp(-1.0 / std::max(1.0, sr * qAttackTime));
    qswitchRelease = std::exp(-1.0 / std::max(1.0, sr * qReleaseTime));
}

void t_kbeyond::render_early(double inL, double inR, double widthNorm, double earlyAmt, double focusAmt, double &earlyL, double &earlyR) {
    widthNorm = clampd(widthNorm, 0.0, 2.0);
    double focusNorm = clampd(focusAmt, 0.0, 1.0);
    double wMainBase = 0.5 * (1.0 + widthNorm);
    double wCrossBase = 0.5 * (1.0 - widthNorm);
    double wMain = lerp(1.0, wMainBase, focusNorm);
    double wCross = lerp(0.0, wCrossBase, focusNorm);
    double midIn = 0.5 * (inL + inR);
    double sideIn = 0.5 * (inL - inR);
    double detector = std::max(std::fabs(midIn), std::fabs(sideIn));
    double envCoef = (detector > laserEnv) ? laserEnvAttack : laserEnvRelease;
    laserEnv = lerp(detector, laserEnv, envCoef);
    double gateBase = clampd((laserEnv - laserGateScaled) / 0.3, 0.0, 1.0);
    double clusterGate = std::sqrt(gateBase);
    double focusGate = clusterGate * focusNorm;
    laserExcite = focusGate;
    earlyBufMid.write(midIn + tiny());
    earlyBufSide.write(sideIn + tiny());
    double sideBlendBase = clampd(0.5 * widthNorm, 0.0, 1.0);
    double sideBlend = lerp(0.0, sideBlendBase, focusNorm);
    double left = 0.0;
    double right = 0.0;
    for (int tap = 0; tap < kEarlyTaps; ++tap) {
        double tapMid = earlyBufMid.readInt(earlyDel[tap]);
        double tapSide = earlyBufSide.readInt(earlyDel[tap]);
        double gain = earlyGain[tap];
        double cosTheta = earlyCos[tap];
        double sinTheta = earlySin[tap];
        double baseL = tapMid * gain * cosTheta;
        double baseR = tapMid * gain * sinTheta;
        double sideL = tapSide * gain * cosTheta;
        double sideR = -tapSide * gain * sinTheta;
        double mixL = baseL + sideBlend * sideL;
        double mixR = baseR + sideBlend * sideR;
        left += wMain * mixL + wCross * mixR;
        right += wCross * mixL + wMain * mixR;
    }
    double clusterAmt = clampd(laser, 0.0, 1.0);
    if (clusterAmt > 0.0) {
        double phaseNow = laserPhase;
        laserPhase += laserPhaseInc;
        if (laserPhase >= 2.0 * M_PI)
            laserPhase -= 2.0 * M_PI;
        double swirl = std::sin(phaseNow);
        double swirlB = std::sin(phaseNow * 0.5 + 1.0471975511965976);
        double swirlMix = 0.5 * (swirl + swirlB);
        double modDepth = lerp(0.25, 0.9, clusterAmt * clusterAmt);
        double gate = focusGate * focusGate;
        if (gate < 1.0e-6)
            gate = 0.0;
        if (gate > 0.0) {
            double clusterL = 0.0;
            double clusterR = 0.0;
            for (int tap = 0; tap < kLaserTaps; ++tap) {
                double tapMid = earlyBufMid.readFrac(laserDelay[tap]);
                double tapSide = earlyBufSide.readFrac(laserDelay[tap]);
                double mod = 1.0 + modDepth * swirlMix * laserShape[tap];
                double midGain = laserMidGain[tap] * mod;
                double sideGain = laserSideGain[tap] * mod;
                double cosTheta = laserCos[tap];
                double sinTheta = laserSin[tap];
                double baseL = tapMid * midGain * cosTheta;
                double baseR = tapMid * midGain * sinTheta;
                double sideL = tapSide * sideGain * cosTheta;
                double sideR = -tapSide * sideGain * sinTheta;
                clusterL += baseL + sideL;
                clusterR += baseR + sideR;
            }
            left += clusterAmt * gate * clusterL;
            right += clusterAmt * gate * clusterR;
        }
    } else {
        laserPhase += laserPhaseInc;
        if (laserPhase >= 2.0 * M_PI)
            laserPhase -= 2.0 * M_PI;
    }
    if (focusNorm < 1.0) {
        double center = 0.5 * (left + right);
        double blend = focusNorm * focusNorm;
        left = lerp(center, left, blend);
        right = lerp(center, right, blend);
    }
    earlyL = left * earlyAmt;
    earlyR = right * earlyAmt;
}

void t_kbeyond::setup_fdn() {
    const double phi = 1.6180339887498948482;
    double baseMs = lerp(15.0, 55.0, size);
    double spread = lerp(0.65, 1.35, size);
    auto latePattern = make_pattern(modeLate, N, 0x2001u);
    for (int i = 0; i < N; ++i) {
        double idx = (i < (int)latePattern.size())
                         ? latePattern[(std::size_t)i]
                         : (double)i / (double)std::max(1, N - 1);
        double ms = baseMs * std::pow(phi, (idx - 0.5) * spread);
        ms = clampd(ms, 8.0, 400.0);
        double len = ms2samp(ms, sr);
        size_t bufLen = (size_t)std::ceil(len + moddepth + 16.0);
        bufLen = std::max<size_t>(bufLen, 64);
        fdn[i].setup(bufLen);
        fdn_len[i] = clampd(len, 8.0, (double)bufLen - 4.0);
        fdn_phase[i] = std::fmod((double)(i + 1) * 1.2345, 2.0 * M_PI);
        fdn_out[i] = 0.0;
        fdn_fb[i] = 0.0;
    }
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

void t_kbeyond::update_output_weights() {
    double widthNorm = clampd(width, 0.0, 2.0);
    auto midPattern = make_pattern(modeMid, N, 0x4001u);
    auto sidePattern = make_pattern(modeSide, N, 0x4002u);
    double normL = 0.0;
    double normR = 0.0;
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

void t_kbeyond::update_modulators() {
    double rate = clampd(modrate, 0.0, 5.0);
    for (int i = 0; i < N; ++i) {
        double idx = (double)(i + 1) / (double)N;
        double warp = lerp(0.4, 1.2, idx);
        fdn_phaseInc[i] = 2.0 * M_PI * rate * warp / sr;
        fdn_phase[i] = std::fmod(fdn_phase[i], 2.0 * M_PI);
        fdn_read[i] = clampd(fdn_len[i], 8.0, (double)fdn[i].size() - 4.0);
    }
}

void t_kbeyond::update_decay() {
    double rt60 = std::max(decay, 0.0);
    double regenNorm = clampd(regen, 0.0, 0.999);
    double dampLFNorm = clampd(dampLF, 0.0, 1.0);
    double dampMFNorm = clampd(dampMF, 0.0, 1.0);
    double dampHFNorm = clampd(dampHF, 0.0, 1.0);
    dampLF_mul = lerp(1.0, 0.35, dampLFNorm);
    dampMF_mul = lerp(1.0, 0.45, dampMFNorm);
    dampHF_mul = lerp(1.0, 0.12, dampHFNorm);
    double minDecay = 0.05;
    for (int i = 0; i < N; ++i) {
        double delaySamples = fdn_len[i];
        double delaySeconds = sr > 0.0 ? delaySamples / sr : 0.0;
        double gain = 0.0;
        if (delaySeconds > 0.0) {
            if (regenNorm <= 0.0) {
                gain = 0.0;
            } else if (rt60 > 0.0) {
                double regenScale = lerp(0.05, 1.0, regenNorm);
                double tau = std::max(rt60 * regenScale, minDecay);
                double exponent = (-3.0 * delaySeconds) / tau;
                gain = std::pow(10.0, exponent);
            } else {
                gain = regenNorm;
            }
        }
        fdn_decay[i] = clampd(gain, 0.0, 0.99995);
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
    double rate = clampd(uwalkRate, 0.0, 8.0);
    double incBase = (rate > 0.0 && sr > 0.0) ? (2.0 * M_PI * rate) / sr : 0.0;
    for (int i = 0; i < N; ++i) {
        double idx = (double)(i + 1) / (double)N;
        double warp = lerp(0.35, 1.75, idx);
        uwalkPhaseInc[i] = incBase * warp;
        if (rate <= 0.0)
            uwalkPhaseInc[i] = 0.0;
        double ditherWarp = lerp(0.55, 1.45, idx);
        qditherPhaseInc[i] = incBase * ditherWarp;
        if (rate <= 0.0)
            qditherPhaseInc[i] = 0.0;
    }
}

void t_kbeyond::reset_quantum_walk() {
    const double g = 0.6180339887498948482; // golden ratio reciprocal for decorrelation
    for (int i = 0; i < N; ++i) {
        double seed = ((double)(i + 1) * g + 0.37) * 2.0 * M_PI;
        uwalkPhase[i] = std::fmod(seed, 2.0 * M_PI);
        if (uwalkPhase[i] < 0.0)
            uwalkPhase[i] += 2.0 * M_PI;
        uwalkState[i] = 0.0;
        double ditherSeed = ((double)(i + 1) * g + 0.11) * 2.0 * M_PI;
        qditherPhase[i] = std::fmod(ditherSeed, 2.0 * M_PI);
        if (qditherPhase[i] < 0.0)
            qditherPhase[i] += 2.0 * M_PI;
    }
    update_quantum_walk();
}

static inline double wrap_phase(double phase) {
    if (phase >= 2.0 * M_PI)
        phase -= 2.0 * M_PI;
    else if (phase < 0.0)
        phase += 2.0 * M_PI;
    return phase;
}

void t_kbeyond::apply_quantum_dither(std::array<double, N> &vector) {
    double coherenceAmt = clampd(coherence, 0.0, 1.0);
    double maxAngle = lerp(0.0, 0.42, coherenceAmt * coherenceAmt);
    bool advancePhaseOnly = (maxAngle <= 0.0);

    auto rotate_pair = [&](int idx, int jdx, double angle) {
        double c = std::cos(angle);
        double s = std::sin(angle);
        double a = vector[idx];
        double b = vector[jdx];
        vector[idx] = c * a - s * b;
        vector[jdx] = s * a + c * b;
    };

    auto advance_phase = [&](int idx) {
        qditherPhase[idx] = wrap_phase(qditherPhase[idx] + qditherPhaseInc[idx]);
    };

    if (!advancePhaseOnly) {
        for (int start = 0; start < 2; ++start) {
            for (int i = start; i < N; i += 2) {
                int j = (i + 1) % N;
                double phase = qditherPhase[i];
                double angle = maxAngle * std::sin(phase);
                rotate_pair(i, j, angle);
                advance_phase(i);
            }
        }
        for (int i = 0; i < N; ++i) {
            int j = (i + 5) % N;
            double phase = qditherPhase[(i + 3) % N];
            double angle = 0.5 * maxAngle * std::sin(phase * 0.5);
            rotate_pair(i, j, angle);
        }
    } else {
        for (int i = 0; i < N; ++i)
            advance_phase(i);
    }
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
        apply_walsh_hadamard16(*source, output);
        break;
    case MixMode::Hybrid:
        apply_hybrid_diffusion(u, *source, output, diffusionScratch);
        break;
    default:
        apply_householder<N>(u, *source, output);
        break;
    }
}

void t_kbeyond::apply_quantum_walk(std::array<double, N> &feedback) {
    double coherenceAmt = clampd(coherence, 0.0, 1.0);
    if (coherenceAmt <= 0.0)
        return;

    std::array<double, N> neighborMix {};
    double stateMix = lerp(0.08, 0.25, coherenceAmt);
    for (int i = 0; i < N; ++i) {
        double phase = uwalkPhase[i];
        double sinA = std::sin(phase);
        double sinB = std::sin(phase * 1.7320508075688772 + (double)i * 0.4115);
        double blend = 0.5 * (sinA + sinB);
        int a = (i + 1) % N;
        int b = (i + 5) % N;
        double neighbor = lerp(feedback[a], feedback[b], 0.5 * (blend + 1.0));
        neighborMix[i] = lerp(feedback[i], neighbor, coherenceAmt);

        double drift = neighborMix[i] * (0.6 * sinA + 0.4 * sinB);
        uwalkState[i] = lerp(uwalkState[i], drift, stateMix);

        uwalkPhase[i] += uwalkPhaseInc[i];
        if (uwalkPhase[i] > 2.0 * M_PI)
            uwalkPhase[i] -= 2.0 * M_PI;
    }

    for (int i = 0; i < N; ++i) {
        double modulation = uwalkState[i] * 0.2;
        feedback[i] = lerp(feedback[i], neighborMix[i] + modulation, coherenceAmt);
    }
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
        x->predSamps = clampd(x->predelay * x->sr, 0.0, std::max(0.0, (double)x->predL.size() - 4.0));
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

    double widthNorm = clampd(x->width, 0.0, 2.0);
    double mix = clampd(x->mix, 0.0, 1.0);
    double earlyAmt = clampd(x->early, 0.0, 1.0);
    double focusAmt = clampd(x->focus, 0.0, 1.0);
    double moddepth = clampd(x->moddepth, 0.0, 32.0);

    for (long i = 0; i < sampleframes; ++i) {
        double inL = (nin > 0 && ins[0]) ? ins[0][i] : 0.0;
        double inR = (nin > 1 && ins[1]) ? ins[1][i] : inL;

        double dryL = inL;
        double dryR = inR;

        double predLen = x->predSamps;
        double predOutL = x->predL.readFrac(predLen);
        double predOutR = x->predR.readFrac(predLen);
        x->predL.write(inL + x->tiny());
        x->predR.write(inR + x->tiny());

        double midIn = 0.5 * (predOutL + predOutR);
        double sideIn = 0.5 * (predOutL - predOutR);

        double earlyL = 0.0;
        double earlyR = 0.0;
        x->render_early(predOutL, predOutR, widthNorm, earlyAmt, focusAmt, earlyL, earlyR);

        double qTarget = x->laserExcite * clampd(x->laser, 0.0, 1.0);
        if (x->qswitchWindowSamples > 0 && x->laserDiffusion > 0.0) {
            double env = x->qswitchEnv;
            double coef = (qTarget > env) ? x->qswitchAttack : x->qswitchRelease;
            env = lerp(qTarget, env, coef);
            x->qswitchEnv = env;
            if (env > 0.01)
                x->qswitchCounter = x->qswitchWindowSamples;
        } else {
            x->qswitchCounter = 0;
            x->qswitchEnv = 0.0;
        }

        std::array<double, t_kbeyond::N> vec {};
        for (int l = 0; l < t_kbeyond::N; ++l) {
            double mod = 0.0;
            if (x->modrate > 0.0 && moddepth > 0.0) {
                mod = std::sin(x->fdn_phase[l]) * moddepth;
                x->fdn_phase[l] += x->fdn_phaseInc[l];
                if (x->fdn_phase[l] > 2.0 * M_PI)
                    x->fdn_phase[l] -= 2.0 * M_PI;
            }
            double read = clampd(x->fdn_len[l] + mod, 2.0, (double)x->fdn[l].size() - 3.0);
            x->fdn_read[l] = read;
            double sig = x->fdn[l].readFrac(read);
            sig = x->fdn_tilt[l].process(sig);
            sig = x->fdn_lp[l].process(sig);
            vec[l] = sig;
            x->fdn_out[l] = sig;
        }

        x->apply_diffusion(vec, x->fdn_fb);
        double qMix = 0.0;
        if (x->qswitchCounter > 0 && x->laserDiffusion > 0.0) {
            double norm = (double)x->qswitchCounter / (double)std::max<long>(1, x->qswitchWindowSamples);
            double envelope = std::sin(0.5 * M_PI * norm);
            qMix = envelope * x->laserDiffusion * clampd(x->qswitchEnv * 1.2, 0.0, 1.0);
            if (qMix < 1.0e-6)
                qMix = 0.0;
            --x->qswitchCounter;
        }
        if (qMix > 0.0) {
            std::array<double, t_kbeyond::N> alt {};
            apply_walsh_hadamard16(vec, alt);
            for (int l = 0; l < t_kbeyond::N; ++l)
                x->fdn_fb[l] = lerp(x->fdn_fb[l], alt[l], qMix);
        }
        x->apply_quantum_walk(x->fdn_fb);

        double tailL = 0.0;
        double tailR = 0.0;
        for (int l = 0; l < t_kbeyond::N; ++l) {
            double fb = x->fdn_fb[l];
            double low = x->fdn_low[l].process(fb);
            double band = x->fdn_high[l].process(fb);
            double high = fb - band;
            double mid = band - low;
            double damped = low * x->dampLF_mul + mid * x->dampMF_mul + high * x->dampHF_mul;
            double feedback = damped * x->fdn_decay[l];
            double injection = midIn * x->inWeights[l] + sideIn * x->outWeightsL[l] * 0.15;
            x->fdn[l].write(feedback + injection + x->tiny());
            tailL += x->fdn_out[l] * x->outWeightsL[l];
            tailR += x->fdn_out[l] * x->outWeightsR[l];
        }

        double wetL = earlyL + tailL;
        double wetR = earlyR + tailR;

        if (outL)
            outL[i] = lerp(dryL, wetL, mix);
        if (outR)
            outR[i] = lerp(dryR, wetR, mix);
    }
}

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

