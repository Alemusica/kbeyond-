# AGENTS.md

> Istruzioni operative per agent (es. Codex) su **kbeyond-**.
> Obiettivo: qualità “flagship” in 5 fasi, senza stravolgere il core FDN.

## 0) Goal & Tier
- **Live (CPU-only)**: early bilanciata, FDN stabile, bassa latenza.
- **Massima Qualità (CPU-HQ)**: RT60 vero, damping 3-bande, mixing HH/WHT/Hybrid, evoluzione senza chorus.
- **Hyper Qualità (HYQ, GPU)**: late tail su GPU (Metal/MPSGraph), decorrelazioni “quantum/laser”.

## 1) Guard-rails
- Niente allocazioni nel thread audio; aggiornamenti per blocco; FTZ/DAZ.
- Param smoothing 5–20 ms; trasformazioni **energy-preserving** (M/S, mixing).
- Core FDN intatto; aggiunte modulari e reversibili.

## 2) Roadmap (fasi → branch/PR)
**A — Hardening (`a-hardening`)**
1. Early energy-balanced (tap specchiati, pan equal-power, width in M/S ortonormale). Test impulso mono → |L−R| < 0.1 dB.
2. RT60 vero `@decay`: per-linea `g_i = 10^(−3·d_i/@decay)`; damping LF/MF/HF (1-pole).
3. Vitalità low-CPU: 2–4 allpass con jitter lento + micro-jitter su pochi tap.

**B — Potenza (`b-power`)**
4. Mixing **WHT/Hybrid** oltre Householder(φ).
5. Detector range–doppler (RMS lento + banda 0.7–3 kHz) che modula `width/hf_damp/modDepth` con slew.

**C — Novità sonore (`c-phonon-hq`, `c-laser`)**
6. **PHONON**: scambi di energia tra micro-modi con rotazioni di Givens a blocchi (unitario). Parametri `@phonon_*`.
7. Early “RADAR-coded” (chirp/m-seq) + **Q-switch** di diffusione 200–600 ms.

**D — PrimeFear (`d-primefear`)**
8. **PrimeModes** (`prime/aureo/plastica/prime_aureo`) per ER/FDN + routing M/S. Guard-rails: co-primalità, distanza ≥ 2 smp, auto-bilanciamento.

**E — HYQ GPU (`e-gpu-latecon`)**
9. Late tail partizionata su GPU (MPSGraph FFT); prima partizione CPU/vDSP; fallback automatico.
10. PHONON/Laser in frequenza (per-bin/per-banda).

## 3) Definition of Done (per PR)
- Early |L−R| < 0.1 dB; RT60 errore ≤ ±5% (per banda); nessun ringing corto.
- Niente chorus/swirl (pitch-tracking della coda pulito).
- CPU: HQ ≈ ≤ 1–1.3%/core @48 kHz (indicativo). GPU: throughput > CPU e fallback testato.
- Docs e preset aggiornati; test QA verdi (vedi `docs/QA.md`).

## 4) Percorsi & Documenti
- Roadmap: `docs/roadmap_it.md`
- QA & metriche: `docs/QA.md`
- PHONON: `docs/PHONON.md`
- Prime modes: `docs/PRIME_MODES.md`
