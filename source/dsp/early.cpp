#include "early.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "dsp_assert.h"

namespace kbeyond::dsp {

namespace {
constexpr double kGoldenRatio = 1.6180339887498948482;
}

void EarlySection::setup(double sr,
                         double size,
                         double focus,
                         const std::vector<double> &tapPattern,
                         const std::vector<double> &laserMidPattern,
                         const std::vector<double> &laserSidePattern) {
    sampleRate_ = sr;
    const double maxMs = 120.0;
    const std::size_t len = static_cast<std::size_t>(std::ceil(ms2samp(maxMs * 1.25, sr) + 8.0));
    earlyBufMid_.setup(len);
    earlyBufSide_.setup(len);

    const double baseMs = 5.2;
    const double scale = lerp(0.7, 1.45, size);
    const int pairs = static_cast<int>(kEarlyTaps / 2);
    int patternCount = pairs;
    if (kEarlyTaps % 2 != 0)
        ++patternCount;
    if (patternCount <= 0)
        patternCount = 1;

    const auto get_pattern = [&](int index) {
        if (index < static_cast<int>(tapPattern.size()))
            return tapPattern[static_cast<std::size_t>(index)];
        if (pairs > 1)
            return static_cast<double>(index) / static_cast<double>(pairs - 1);
        return 0.5;
    };

    const auto storeTap = [&](int index, long samp, double gain, double pan) {
        const double theta = (pan + 1.0) * (0.25 * M_PI);
        earlyDel_[index] = samp;
        earlyGain_[index] = gain;
        earlyCos_[index] = std::cos(theta);
        earlySin_[index] = std::sin(theta);
    };

    dsp_assert_msg(tapPattern.empty() || std::is_sorted(tapPattern.begin(), tapPattern.end()),
                   "Tap pattern must be sorted");

    std::vector<long> uniqueTapDelays;
    uniqueTapDelays.reserve(static_cast<std::size_t>(patternCount));
    const long minTapSpacing = 2;
    const long maxSample = static_cast<long>(std::floor(static_cast<double>(earlyBufMid_.size()) - 2.0));

    const auto is_valid_delay = [&](long candidate) {
        for (long prev : uniqueTapDelays) {
            if (std::llabs(candidate - prev) < minTapSpacing)
                return false;
        }
        return true;
    };

    const auto resolve_delay = [&](double samples) {
        const long minSample = 1;
        long base = clampl(static_cast<long>(std::floor(samples)), minSample, maxSample);
        constexpr int offsets[] = {0, 1, -1, 2, -2, 3, -3, 4, -4};
        for (int offset : offsets) {
            long candidate = clampl(base + offset, minSample, maxSample);
            if (is_valid_delay(candidate)) {
                uniqueTapDelays.push_back(candidate);
                return candidate;
            }
        }
        for (long delta = minTapSpacing; delta < minTapSpacing + 32; ++delta) {
            long up = clampl(base + delta, minSample, maxSample);
            if (is_valid_delay(up)) {
                uniqueTapDelays.push_back(up);
                return up;
            }
            long down = clampl(base - delta, minSample, maxSample);
            if (is_valid_delay(down)) {
                uniqueTapDelays.push_back(down);
                return down;
            }
        }
        uniqueTapDelays.push_back(base);
        return base;
    };

    for (int p = 0; p < pairs; ++p) {
        const double idxNorm = get_pattern(p);
        double ms = baseMs * std::pow(kGoldenRatio, idxNorm * 1.2);
        ms = std::min(ms * scale, maxMs);
        const long samp = resolve_delay(ms2samp(ms, sr));
        const double gain = std::pow(0.72, static_cast<double>(p) + 1.0);
        const double panMag = 1.0 - idxNorm;
        const double panLeft = -panMag;
        const double panRight = panMag;
        storeTap(p, samp, gain, panLeft);
        storeTap(static_cast<int>(kEarlyTaps) - 1 - p, samp, gain, panRight);
    }

    if (kEarlyTaps % 2 != 0) {
        const int center = pairs;
        const double idxNorm = (center < static_cast<int>(tapPattern.size())) ? tapPattern[static_cast<std::size_t>(center)] : 0.65;
        double ms = baseMs * std::pow(kGoldenRatio, idxNorm);
        ms = std::min(ms * scale, maxMs);
        const long samp = resolve_delay(ms2samp(ms, sr));
        const double gain = std::pow(0.72, static_cast<double>(pairs) + 1.0);
        storeTap(center, samp, gain, 0.0);
    }

    if (!uniqueTapDelays.empty()) {
        auto sorted = uniqueTapDelays;
        std::sort(sorted.begin(), sorted.end());
        for (std::size_t i = 1; i < sorted.size(); ++i)
            dsp_assert_msg(sorted[i] - sorted[i - 1] >= minTapSpacing, "Tap spacing below minimum");
    }

    double energy = 0.0;
    for (int tap = 0; tap < static_cast<int>(kEarlyTaps); ++tap) {
        double wL = earlyGain_[static_cast<std::size_t>(tap)] * earlyCos_[static_cast<std::size_t>(tap)];
        double wR = earlyGain_[static_cast<std::size_t>(tap)] * earlySin_[static_cast<std::size_t>(tap)];
        energy += wL * wL + wR * wR;
    }
    dsp_assert_msg(energy > 0.0, "Early tap energy invalid");
    if (energy > 0.0) {
        const double scaleEnergy = 1.0 / std::sqrt(energy);
        for (double &gain : earlyGain_)
            gain *= scaleEnergy;
    }

    const double focusNorm = clampd(focus, 0.0, 1.0);
    const double baseLaser = lerp(6.5, 17.5, size);
    const double groupSpacing = lerp(3.0, 9.5, size);
    const double clusterSpan = lerp(2.2, 7.8, focusNorm);
    const double chirpSpan = lerp(0.6, 4.2, focusNorm);
    const double jitterAmt = lerp(0.08, 0.35, focusNorm);

    uint32_t lfsr = 0x1Fu;
    const auto next_code = [&]() {
        const int bit = ((lfsr >> 4) ^ (lfsr >> 2)) & 1;
        lfsr = ((lfsr << 1) & 0x1F) | static_cast<uint32_t>(bit);
        return (lfsr & 1U) ? 1.0 : -1.0;
    };

    const double maxDelay = static_cast<double>(earlyBufMid_.size()) - 2.0;
    for (int idx = 0; idx < static_cast<int>(kLaserTaps); ++idx) {
        const double norm = (kLaserTaps > 1) ? static_cast<double>(idx) / static_cast<double>(kLaserTaps - 1) : 0.0;
        const double groupFloat = norm * static_cast<double>(kLaserGroups);
        int group = static_cast<int>(std::floor(groupFloat));
        if (group >= static_cast<int>(kLaserGroups))
            group = static_cast<int>(kLaserGroups) - 1;
        const double local = groupFloat - static_cast<double>(group);
        const double startMs = baseLaser + groupSpacing * static_cast<double>(group);
        const double sweep = clusterSpan * (1.0 + 0.27 * static_cast<double>(group));
        const double chirp = chirpSpan * local * local;
        const double jitter = jitterAmt * next_code() * (1.0 - local);
        const double tapMs = clampd(startMs + local * sweep + chirp + jitter, 1.0, maxMs - 1.0);
        const double delaySamp = clampd(ms2samp(tapMs, sr), 1.0, maxDelay);
        laserDelay_[static_cast<std::size_t>(idx)] = delaySamp;

        const double shape = next_code();
        laserShape_[static_cast<std::size_t>(idx)] = shape;
        const double baseGain = 0.22 * std::pow(0.78, static_cast<double>(idx) * 0.55);
        const double bright = lerp(0.78, 1.32, focusNorm) * (1.0 + 0.22 * shape);
        const double midScale = (idx < static_cast<int>(laserMidPattern.size()))
                                    ? lerp(0.65, 1.45, laserMidPattern[static_cast<std::size_t>(idx)])
                                    : 1.0;
        laserMidGain_[static_cast<std::size_t>(idx)] = baseGain * bright * midScale;
        const double sideSign = next_code();
        const double sideScale = (idx < static_cast<int>(laserSidePattern.size()))
                                     ? lerp(0.65, 1.45, laserSidePattern[static_cast<std::size_t>(idx)])
                                     : 1.0;
        const double sideWeight = lerp(0.18, 0.55, focusNorm) * sideSign * sideScale;
        laserSideGain_[static_cast<std::size_t>(idx)] = baseGain * sideWeight;
        const double panBase = lerp(-0.9, 0.9, norm);
        const double panJitter = 0.25 * next_code() * (1.0 - local);
        const double pan = clampd(panBase + panJitter, -1.0, 1.0);
        const double theta = (pan + 1.0) * (0.25 * M_PI);
        laserCos_[static_cast<std::size_t>(idx)] = std::cos(theta);
        laserSin_[static_cast<std::size_t>(idx)] = std::sin(theta);
    }
}

void EarlySection::resetState() {
    laserEnv_ = 0.0;
    laserExcite_ = 0.0;
    qswitchEnv_ = 0.0;
    qswitchCounter_ = 0;
    laserPhase_ = 0.0;
}

void EarlySection::updateGate(double gate) {
    const double g = clampd(gate, 0.0, 1.0);
    laserGateScaled_ = 0.004 + 0.25 * g * g;
}

void EarlySection::updateWindow(double windowSeconds, double sr) {
    sampleRate_ = sr;
    const double winSec = clampd(windowSeconds, 0.2, 0.6);
    qswitchWindowSamples_ = static_cast<long>(std::llround(winSec * sr));
    if (qswitchWindowSamples_ < 1)
        qswitchWindowSamples_ = 1;
    if (qswitchCounter_ > qswitchWindowSamples_)
        qswitchCounter_ = qswitchWindowSamples_;
}

void EarlySection::updatePhaseIncrement(double focus, double sr) {
    sampleRate_ = sr;
    const double focusNorm = clampd(focus, 0.0, 1.0);
    const double sweepHz = lerp(0.5, 7.0, focusNorm);
    laserPhaseInc_ = 2.0 * M_PI * sweepHz / std::max(1.0, sr);
}

void EarlySection::updateEnvelopeCoefficients(double sr) {
    sampleRate_ = sr;
    const double attackTime = 0.0035;
    const double releaseTime = 0.12;
    laserEnvAttack_ = std::exp(-1.0 / std::max(1.0, sr * attackTime));
    laserEnvRelease_ = std::exp(-1.0 / std::max(1.0, sr * releaseTime));
    const double qAttackTime = 0.004;
    const double qReleaseTime = 0.24;
    qswitchAttack_ = std::exp(-1.0 / std::max(1.0, sr * qAttackTime));
    qswitchRelease_ = std::exp(-1.0 / std::max(1.0, sr * qReleaseTime));
}

void EarlySection::render(double inL,
                          double inR,
                          double widthNorm,
                          double earlyAmt,
                          double focusNorm,
                          double clusterAmt,
                          double &earlyL,
                          double &earlyR,
                          uint32_t &rng) {
    widthNorm = clampd(widthNorm, 0.0, 2.0);
    focusNorm = clampd(focusNorm, 0.0, 1.0);
    const double wMainBase = 0.5 * (1.0 + widthNorm);
    const double wCrossBase = 0.5 * (1.0 - widthNorm);
    const double wMain = lerp(1.0, wMainBase, focusNorm);
    const double wCross = lerp(0.0, wCrossBase, focusNorm);
    const double midIn = 0.5 * (inL + inR);
    const double sideIn = 0.5 * (inL - inR);
    const double detector = std::max(std::fabs(midIn), std::fabs(sideIn));
    const double envCoef = (detector > laserEnv_) ? laserEnvAttack_ : laserEnvRelease_;
    laserEnv_ = lerp(detector, laserEnv_, envCoef);
    const double gateBase = clampd((laserEnv_ - laserGateScaled_) / 0.3, 0.0, 1.0);
    const double clusterGate = std::sqrt(gateBase);
    const double focusGate = clusterGate * focusNorm;
    laserExcite_ = focusGate;

    earlyBufMid_.write(midIn + tiny_noise(rng));
    earlyBufSide_.write(sideIn + tiny_noise(rng));

    const double sideBlendBase = clampd(0.5 * widthNorm, 0.0, 1.0);
    const double sideBlend = lerp(0.0, sideBlendBase, focusNorm);

    double left = 0.0;
    double right = 0.0;
    for (int tap = 0; tap < static_cast<int>(kEarlyTaps); ++tap) {
        const double tapMid = earlyBufMid_.readInt(earlyDel_[static_cast<std::size_t>(tap)]);
        const double tapSide = earlyBufSide_.readInt(earlyDel_[static_cast<std::size_t>(tap)]);
        const double gain = earlyGain_[static_cast<std::size_t>(tap)];
        const double cosTheta = earlyCos_[static_cast<std::size_t>(tap)];
        const double sinTheta = earlySin_[static_cast<std::size_t>(tap)];
        const double baseL = tapMid * gain * cosTheta;
        const double baseR = tapMid * gain * sinTheta;
        const double sideL = tapSide * gain * cosTheta;
        const double sideR = -tapSide * gain * sinTheta;
        const double mixL = baseL + sideBlend * sideL;
        const double mixR = baseR + sideBlend * sideR;
        left += wMain * mixL + wCross * mixR;
        right += wCross * mixL + wMain * mixR;
    }

    const double clusterAmount = clampd(clusterAmt, 0.0, 1.0);
    if (clusterAmount > 0.0) {
        const double phaseNow = laserPhase_;
        laserPhase_ += laserPhaseInc_;
        if (laserPhase_ >= 2.0 * M_PI)
            laserPhase_ -= 2.0 * M_PI;
        const double swirl = std::sin(phaseNow);
        const double swirlB = std::sin(phaseNow * 0.5 + 1.0471975511965976);
        const double swirlMix = 0.5 * (swirl + swirlB);
        const double modDepth = lerp(0.25, 0.9, clusterAmount * clusterAmount);
        double gate = focusGate * focusGate;
        if (gate < 1.0e-6)
            gate = 0.0;
        if (gate > 0.0) {
            double clusterL = 0.0;
            double clusterR = 0.0;
            for (int tap = 0; tap < static_cast<int>(kLaserTaps); ++tap) {
                const double tapMid = earlyBufMid_.readFrac(laserDelay_[static_cast<std::size_t>(tap)]);
                const double tapSide = earlyBufSide_.readFrac(laserDelay_[static_cast<std::size_t>(tap)]);
                const double mod = 1.0 + modDepth * swirlMix * laserShape_[static_cast<std::size_t>(tap)];
                const double midGain = laserMidGain_[static_cast<std::size_t>(tap)] * mod;
                const double sideGain = laserSideGain_[static_cast<std::size_t>(tap)] * mod;
                const double cosTheta = laserCos_[static_cast<std::size_t>(tap)];
                const double sinTheta = laserSin_[static_cast<std::size_t>(tap)];
                const double baseL = tapMid * midGain * cosTheta;
                const double baseR = tapMid * midGain * sinTheta;
                const double sideL = tapSide * sideGain * cosTheta;
                const double sideR = -tapSide * sideGain * sinTheta;
                clusterL += baseL + sideL;
                clusterR += baseR + sideR;
            }
            left += clusterAmount * gate * clusterL;
            right += clusterAmount * gate * clusterR;
        }
    } else {
        laserPhase_ += laserPhaseInc_;
        if (laserPhase_ >= 2.0 * M_PI)
            laserPhase_ -= 2.0 * M_PI;
    }

    if (focusNorm < 1.0) {
        const double center = 0.5 * (left + right);
        const double blend = focusNorm * focusNorm;
        left = lerp(center, left, blend);
        right = lerp(center, right, blend);
    }

    earlyL = left * earlyAmt;
    earlyR = right * earlyAmt;
}

double EarlySection::computeQSwitchMix(double clusterAmt, double diffusion) {
    const double qTarget = laserExcite_ * clampd(clusterAmt, 0.0, 1.0);
    if (qswitchWindowSamples_ > 0 && diffusion > 0.0) {
        double env = qswitchEnv_;
        const double coef = (qTarget > env) ? qswitchAttack_ : qswitchRelease_;
        env = lerp(qTarget, env, coef);
        qswitchEnv_ = env;
        if (env > 0.01)
            qswitchCounter_ = qswitchWindowSamples_;
    } else {
        qswitchCounter_ = 0;
        qswitchEnv_ = 0.0;
    }

    double qMix = 0.0;
    if (qswitchCounter_ > 0 && diffusion > 0.0) {
        const double norm = static_cast<double>(qswitchCounter_) /
                            static_cast<double>(std::max<long>(1, qswitchWindowSamples_));
        const double envelope = std::sin(0.5 * M_PI * norm);
        qMix = envelope * diffusion * clampd(qswitchEnv_ * 1.2, 0.0, 1.0);
        if (qMix < 1.0e-6)
            qMix = 0.0;
        --qswitchCounter_;
    }
    return qMix;
}

} // namespace kbeyond::dsp

