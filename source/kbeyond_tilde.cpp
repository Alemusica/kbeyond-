#include "kbeyond_tilde.h"
#include "householder_phi16.h"

#include <algorithm>
#include <numeric>
#include <cstdio>
#include <cstddef>

static constexpr std::size_t kAssistStringMax = 256;

static t_class* s_kbeyond_class = nullptr;

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

void t_kbeyond::setup_sr(double newsr) {
    sr = newsr > 1.0 ? newsr : 48000.0;
    vs = std::max<long>(64, vs);
    filterDynamic = clampd(filter, 0.0, 1.0);
    widthDynamic = clampd(width, 0.0, 2.0);
    moddepthDynamic = clampd(moddepth, 0.0, 32.0);
    setup_predelay();
    setup_early();
    setup_fdn();
    refresh_filters();
    update_householder();
    update_output_weights();
    update_modulators();
    setup_detector();
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
    earlyBuf.setup(len);
    const double phi = 1.6180339887498948482;
    double baseMs = 5.2;
    double scale = lerp(0.7, 1.45, size);
    for (int i = 0; i < kEarlyTaps; ++i) {
        double idx = (double)i / (double)std::max(1, kEarlyTaps - 1);
        double ms = baseMs * std::pow(phi, idx * 1.2);
        ms = std::min(ms * scale, maxMs);
        long samp = (long)clampd(std::floor(ms2samp(ms, sr)), 1.0, (double)earlyBuf.size() - 2.0);
        earlyDel[i] = samp;
        double decay = std::pow(0.72, (double)i + 1.0);
        earlyGain[i] = decay;
        earlyPan[i] = lerp(-1.0, 1.0, idx);
    }
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
        double extra = std::max(32.0, moddepth + 16.0);
        size_t bufLen = (size_t)std::ceil(len + extra);
        bufLen = std::max<size_t>(bufLen, 64);
        fdn[i].setup(bufLen);
        fdn_len[i] = clampd(len, 8.0, (double)bufLen - 4.0);
        fdn_phase[i] = std::fmod((double)(i + 1) * 1.2345, 2.0 * M_PI);
        fdn_out[i] = 0.0;
        fdn_fb[i] = 0.0;
        double angle = (2.0 * M_PI * (double)i) / (double)N;
        inWeights[i] = 0.9 / (double)N * (1.0 + 0.35 * std::sin(angle * 0.73 + 0.2));
    }
}

void t_kbeyond::refresh_filters() {
    double hf = lerp(1200.0, sr * 0.45, clampd(filterDynamic, 0.0, 1.0));
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

void t_kbeyond::update_householder() {
    make_phi_vector<N>(u, clampd(phiweight, 0.0, 1.0));
}

void t_kbeyond::update_output_weights() {
    double widthNorm = clampd(width, 0.0, 2.0);
    for (int i = 0; i < N; ++i) {
        double angle = (2.0 * M_PI * (double)i) / (double)N;
        double l = std::sin(angle * 0.91 + 0.17);
        double r = std::cos(angle * 1.07 - 0.11);
        double mid = 0.5 * (l + r);
        double side = 0.5 * (l - r);
        baseMid[i] = mid;
        baseSide[i] = side;
    }
    apply_width_runtime(widthNorm);
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

void t_kbeyond::setup_detector() {
    double srSafe = std::max(1.0, sr);
    detectorLP700.setCutoff(srSafe, 700.0);
    detectorLP3k.setCutoff(srSafe, 3000.0);
    detectorLP700.z = 0.0;
    detectorLP3k.z = 0.0;
    detectorFast = detectorSlow = 0.0;
    auto timeToCoeff = [srSafe](double seconds) {
        double t = std::max(1e-4, seconds);
        return 1.0 - std::exp(-1.0 / (t * srSafe));
    };
    detectorFastCoeff = timeToCoeff(0.03);
    detectorSlowCoeff = timeToCoeff(0.45);
    motionSlewCoeff = timeToCoeff(0.12);
    widthSlewCoeff = timeToCoeff(0.08);
    filterSlewCoeff = timeToCoeff(0.25);
    moddepthSlewCoeff = timeToCoeff(0.18);
    motionState = 0.0;
    widthDynamic = clampd(widthDynamic, 0.0, 2.0);
    filterDynamic = clampd(filterDynamic, 0.0, 1.0);
    moddepthDynamic = clampd(moddepthDynamic, 0.0, 32.0);
}

void t_kbeyond::apply_width_runtime(double widthValue) {
    widthDynamic = clampd(widthValue, 0.0, 2.0);
    for (int i = 0; i < N; ++i) {
        double mid = baseMid[i];
        double side = baseSide[i];
        outWeightsL[i] = mid + widthDynamic * side;
        outWeightsR[i] = mid - widthDynamic * side;
    }
}

void t_kbeyond::apply_filter_runtime(double normalizedFilter) {
    filterDynamic = clampd(normalizedFilter, 0.0, 1.0);
    refresh_filters();
}

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
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->regen, 0.0, 0.999, nullptr);
}

t_max_err kbeyond_attr_set_derez(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->derez, 0.0, 1.0, &t_kbeyond::refresh_filters);
}

t_max_err kbeyond_attr_set_filter(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->filter, 0.0, 1.0, nullptr);
    if (!err) {
        x->apply_filter_runtime(x->filter);
    }
    return err;
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

t_max_err kbeyond_attr_set_modrate(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->modrate, 0.0, 5.0, &t_kbeyond::update_modulators);
}

t_max_err kbeyond_attr_set_moddepth(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    t_max_err err = kbeyond_attr_set_double(x, attr, argc, argv, &x->moddepth, 0.0, 32.0, &t_kbeyond::setup_fdn);
    if (!err)
        x->update_modulators();
    if (!err)
        x->moddepthDynamic = clampd(x->moddepth, 0.0, 32.0);
    return err;
}

t_max_err kbeyond_attr_set_phiweight(t_kbeyond* x, void* attr, long argc, t_atom* argv) {
    return kbeyond_attr_set_double(x, attr, argc, argv, &x->phiweight, 0.0, 1.0, &t_kbeyond::update_householder);
}

t_max_err kbeyond_attr_set_mode_mix(t_kbeyond* x, void*, long argc, t_atom* argv) {
    if (!x)
        return MAX_ERR_GENERIC;
    t_symbol* sym = nullptr;
    if (argc > 0) {
        if (atom_gettype(argv) == A_SYM)
            sym = atom_getsym(argv);
        else if (atom_gettype(argv) == A_LONG)
            sym = atom_getlong(argv) == 0 ? gensym("householder") : gensym("wht");
    }
    if (!sym)
        sym = gensym("householder");

    x->mode_mix = sym;
    if (sym == gensym("wht") || sym == gensym("WHT")) {
        x->mixMode = t_kbeyond::MixMode::WHT;
    } else if (sym == gensym("hybrid") || sym == gensym("Hybrid")) {
        x->mixMode = t_kbeyond::MixMode::Hybrid;
    } else {
        x->mixMode = t_kbeyond::MixMode::Householder;
        x->mode_mix = gensym("householder");
    }
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

    x->mode_mix = gensym("householder");
    x->mixMode = t_kbeyond::MixMode::Householder;

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

    double mix = clampd(x->mix, 0.0, 1.0);
    double earlyAmt = clampd(x->early, 0.0, 1.0);
    double regen = clampd(x->regen, 0.0, 0.999);
    double baseWidth = clampd(x->width, 0.0, 2.0);
    double baseFilter = clampd(x->filter, 0.0, 1.0);
    double baseModdepth = clampd(x->moddepth, 0.0, 32.0);

    double widthDyn = clampd(x->widthDynamic, 0.0, 2.0);
    double filterDyn = clampd(x->filterDynamic, 0.0, 1.0);
    double moddepthDyn = clampd(x->moddepthDynamic, 0.0, 32.0);
    double motionState = clampd(x->motionState, 0.0, 1.0);

    double widthSlew = x->widthSlewCoeff;
    double filterSlew = x->filterSlewCoeff;
    double moddepthSlew = x->moddepthSlewCoeff;
    double motionSlew = x->motionSlewCoeff;
    double detectorFast = x->detectorFast;
    double detectorSlow = x->detectorSlow;
    double detectorFastCoeff = x->detectorFastCoeff;
    double detectorSlowCoeff = x->detectorSlowCoeff;

    OnePoleLP &lp700 = x->detectorLP700;
    OnePoleLP &lp3k = x->detectorLP3k;

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

        double low700 = lp700.process(midIn);
        double band = lp3k.process(midIn - low700);
        double power = band * band;
        detectorFast += detectorFastCoeff * (power - detectorFast);
        detectorSlow += detectorSlowCoeff * (power - detectorSlow);
        double fastEnergy = std::max(detectorFast, 1e-12);
        double diff = std::max(0.0, detectorFast - detectorSlow);
        double motionRaw = clampd(diff / (fastEnergy + 1e-12), 0.0, 1.0);
        motionRaw = std::sqrt(motionRaw);
        motionState += motionSlew * (motionRaw - motionState);

        double widthTarget = clampd(baseWidth + motionState * (2.0 - baseWidth) * 0.4, 0.0, 2.0);
        widthDyn += widthSlew * (widthTarget - widthDyn);
        x->apply_width_runtime(widthDyn);
        widthDyn = x->widthDynamic;
        double widthNorm = widthDyn;

        double filterTarget = clampd(baseFilter + motionState * (1.0 - baseFilter) * 0.35, 0.0, 1.0);
        filterDyn += filterSlew * (filterTarget - filterDyn);
        double filterCandidate = clampd(filterDyn, 0.0, 1.0);
        if (std::fabs(filterCandidate - x->filterDynamic) > 1e-4)
            x->apply_filter_runtime(filterCandidate);
        filterDyn = x->filterDynamic;

        double moddepthTarget = clampd(baseModdepth + motionState * (32.0 - baseModdepth) * 0.5, 0.0, 32.0);
        moddepthDyn += moddepthSlew * (moddepthTarget - moddepthDyn);
        double moddepth = clampd(moddepthDyn, 0.0, 32.0);
        x->moddepthDynamic = moddepth;

        x->earlyBuf.write(midIn + x->tiny());
        double earlyL = 0.0;
        double earlyR = 0.0;
        for (int tap = 0; tap < t_kbeyond::kEarlyTaps; ++tap) {
            double tapSample = x->earlyBuf.readInt(x->earlyDel[tap]);
            double energy = tapSample * x->earlyGain[tap];
            double pan = x->earlyPan[tap];
            double side = energy * pan * widthNorm;
            double mid = energy;
            earlyL += mid + side;
            earlyR += mid - side;
        }
        earlyL *= earlyAmt;
        earlyR *= earlyAmt;

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

        switch (x->mixMode) {
        case t_kbeyond::MixMode::Householder:
            apply_householder<t_kbeyond::N>(x->u, vec, x->fdn_fb);
            break;
        case t_kbeyond::MixMode::WHT:
            apply_walsh_hadamard<t_kbeyond::N>(vec, x->fdn_fb);
            break;
        case t_kbeyond::MixMode::Hybrid: {
            std::array<double, t_kbeyond::N> hh {};
            std::array<double, t_kbeyond::N> wht {};
            apply_householder<t_kbeyond::N>(x->u, vec, hh);
            apply_walsh_hadamard<t_kbeyond::N>(vec, wht);
            double blend = clampd(x->phiweight, 0.0, 1.0);
            for (int l = 0; l < t_kbeyond::N; ++l)
                x->fdn_fb[l] = lerp(wht[l], hh[l], blend);
            break;
        }
        }

        double tailL = 0.0;
        double tailR = 0.0;
        for (int l = 0; l < t_kbeyond::N; ++l) {
            double feedback = x->fdn_fb[l] * regen;
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

    x->widthDynamic = widthDyn;
    x->filterDynamic = filterDyn;
    x->moddepthDynamic = clampd(moddepthDyn, 0.0, 32.0);
    x->motionState = clampd(motionState, 0.0, 1.0);
    x->detectorFast = detectorFast;
    x->detectorSlow = detectorSlow;
}

extern "C" C74_EXPORT void ext_main(void *r) {
    t_class *c = class_new("kbeyond~", (method)kbeyond_new, (method)kbeyond_free, sizeof(t_kbeyond), 0L, A_GIMME, 0);

    class_addmethod(c, (method)kbeyond_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)kbeyond_assist, "assist", A_CANT, 0);

    CLASS_ATTR_DOUBLE(c, "regen", 0, t_kbeyond, regen);
    CLASS_ATTR_ACCESSORS(c, "regen", NULL, kbeyond_attr_set_regen);
    CLASS_ATTR_FILTER_CLIP(c, "regen", 0.0, 0.999);

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

    CLASS_ATTR_DOUBLE(c, "modrate", 0, t_kbeyond, modrate);
    CLASS_ATTR_ACCESSORS(c, "modrate", NULL, kbeyond_attr_set_modrate);
    CLASS_ATTR_FILTER_CLIP(c, "modrate", 0.0, 5.0);

    CLASS_ATTR_DOUBLE(c, "moddepth", 0, t_kbeyond, moddepth);
    CLASS_ATTR_ACCESSORS(c, "moddepth", NULL, kbeyond_attr_set_moddepth);
    CLASS_ATTR_FILTER_CLIP(c, "moddepth", 0.0, 32.0);

    CLASS_ATTR_DOUBLE(c, "phiweight", 0, t_kbeyond, phiweight);
    CLASS_ATTR_ACCESSORS(c, "phiweight", NULL, kbeyond_attr_set_phiweight);
    CLASS_ATTR_FILTER_CLIP(c, "phiweight", 0.0, 1.0);

    CLASS_ATTR_SYM(c, "mode_mix", 0, t_kbeyond, mode_mix);
    CLASS_ATTR_ACCESSORS(c, "mode_mix", NULL, kbeyond_attr_set_mode_mix);
    CLASS_ATTR_ENUM(c, "mode_mix", 0, "householder wht hybrid");

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    s_kbeyond_class = c;
}

