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
    setup_fdn();
    refresh_filters();
    update_diffusion();
    update_output_weights();
    update_modulators();
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
    auto storeTap = [&](int index, long samp, double gain, double pan) {
        double theta = (pan + 1.0) * (0.25 * M_PI);
        earlyDel[index] = samp;
        earlyGain[index] = gain;
        earlyCos[index] = std::cos(theta);
        earlySin[index] = std::sin(theta);
    };
    for (int p = 0; p < pairs; ++p) {
        double idx = (double)p / (double)std::max(1, pairs - 1);
        double ms = baseMs * std::pow(phi, idx * 1.2);
        ms = std::min(ms * scale, maxMs);
        long samp = (long)clampd(std::floor(ms2samp(ms, sr)), 1.0, (double)earlyBufMid.size() - 2.0);
        double gain = std::pow(0.72, (double)p + 1.0);
        double panMag = 1.0 - idx;
        double panLeft = -panMag;
        double panRight = panMag;
        storeTap(p, samp, gain, panLeft);
        storeTap(kEarlyTaps - 1 - p, samp, gain, panRight);
    }
    if (kEarlyTaps % 2 != 0) {
        int center = pairs;
        double ms = baseMs * std::pow(phi, 0.65);
        ms = std::min(ms * scale, maxMs);
        long samp = (long)clampd(std::floor(ms2samp(ms, sr)), 1.0, (double)earlyBufMid.size() - 2.0);
        double gain = std::pow(0.72, (double)pairs + 1.0);
        storeTap(center, samp, gain, 0.0);
    }
}

void t_kbeyond::render_early(double inL, double inR, double widthNorm, double earlyAmt, double &earlyL, double &earlyR) {
    widthNorm = clampd(widthNorm, 0.0, 2.0);
    double wMain = 0.5 * (1.0 + widthNorm);
    double wCross = 0.5 * (1.0 - widthNorm);
    double midIn = 0.5 * (inL + inR);
    double sideIn = 0.5 * (inL - inR);
    earlyBufMid.write(midIn + tiny());
    earlyBufSide.write(sideIn + tiny());
    double sideBlend = clampd(0.5 * widthNorm, 0.0, 1.0);
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
    earlyL = left * earlyAmt;
    earlyR = right * earlyAmt;
}

void t_kbeyond::setup_fdn() {
    const double phi = 1.6180339887498948482;
    double baseMs = lerp(15.0, 55.0, size);
    double spread = lerp(0.65, 1.35, size);
    for (int i = 0; i < N; ++i) {
        double idx = (double)i / (double)std::max(1, N - 1);
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
        double angle = (2.0 * M_PI * (double)i) / (double)N;
        inWeights[i] = 0.9 / (double)N * (1.0 + 0.35 * std::sin(angle * 0.73 + 0.2));
    }
    update_decay();
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

void t_kbeyond::update_diffusion() {
    make_phi_vector<N>(u, clampd(phiweight, 0.0, 1.0));
}

void t_kbeyond::update_output_weights() {
    double widthNorm = clampd(width, 0.0, 2.0);
    double normL = 0.0;
    double normR = 0.0;
    for (int i = 0; i < N; ++i) {
        double angle = (2.0 * M_PI * (double)i) / (double)N;
        double l = std::sin(angle * 0.91 + 0.17);
        double r = std::cos(angle * 1.07 - 0.11);
        double mid = 0.5 * (l + r);
        double side = 0.5 * (l - r);
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

void t_kbeyond::apply_diffusion(const std::array<double, N> &input, std::array<double, N> &output) {
    switch (mixMode) {
    case MixMode::Householder:
        apply_householder<N>(u, input, output);
        break;
    case MixMode::WHT:
        apply_walsh_hadamard16(input, output);
        break;
    case MixMode::Hybrid:
        apply_hybrid_diffusion(u, input, output, diffusionScratch);
        break;
    default:
        apply_householder<N>(u, input, output);
        break;
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
        x->render_early(predOutL, predOutR, widthNorm, earlyAmt, earlyL, earlyR);

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

    CLASS_ATTR_SYM(c, "mode_mix", 0, t_kbeyond, modeMixSym);
    CLASS_ATTR_ACCESSORS(c, "mode_mix", kbeyond_attr_get_mode_mix, kbeyond_attr_set_mode_mix);
    CLASS_ATTR_ENUM(c, "mode_mix", 0, "householder wht hybrid");
    CLASS_ATTR_LABEL(c, "mode_mix", 0, "Diffusion Mode");

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    s_kbeyond_class = c;
}

#endif // KBEYOND_UNIT_TEST

