# Roadmap di sviluppo — kbeyond~

> **Obiettivo:** portare `kbeyond~` a qualità “flagship” senza stravolgere il core, e rendere il codice **modulare** (no megafile) per lavorare rapidi e sicuri (umani e agent).

---

## 0) FACTORIZATION-FIRST (prima di tutto)

### 0.1 Principi
- **Equivalenza funzionale**: la fattorizzazione non deve cambiare il suono (bit-identico entro ±1 LSB).
- **Zero logica nuova**: solo split/nomi/include.
- **Boundary chiari**: header, dipendenze minime, invarianti documentati.

### 0.2 Mappa moduli (target)
source/
├─ kbeyond_wrapper_max.cpp
├─ dsp/
│ ├─ params.h, buffers.h
│ ├─ filters.{h,cpp}
│ ├─ mod.{h,cpp}
│ ├─ early.{h,cpp}
│ ├─ fdn.{h,cpp}
│ ├─ mixing.{h,cpp}
│ ├─ decay.{h,cpp}
│ ├─ phonon.{h,cpp} (stub)
│ ├─ primemodes.{h,cpp} (stub)
│ └─ detector.{h,cpp}
├─ gpu/
│ ├─ latecon.metal (stub)
│ └─ latecon.mm (stub)
└─ util/dsp_assert.h

markdown
Copy code

### 0.3 Sequenza split (PR meccaniche) — `f-factorize`
1) `params.h`, `buffers.h` • 2) `filters.*` • 3) `early.*` • 4) `fdn.*` • 5) `mixing.*` • 6) `decay.*` (RT60 pronto ma disattivato).  
**DoD F0:** A/B IR 30 s → diff RMS < −100 dBFS; CPU ~ invariata; build ok.

---

## 1) Fase A — Hardening (bugfix + controllabilità)

**A1. Early bilanciata** — tap ER specchiati L/R; pan equal-power; `width` in **M/S ortonormale**. Test: |L−R| < 0,1 dB.  
**A2. RT60 reale + Material (3 bande)** — `@decay (s)` con `g_i = 10^(−3·d_i/@decay)` (prima del mixing) + 1-pole LF/MF/HF.  
**A3. Vitalità low-CPU** — 2–4 allpass diffusori con jitter lento + micro-jitter su pochi tap.  
**Nota semplificazione** — **Laser OFF** (temporaneo, no-op sugli attributi).
  - **Stato branch `a-hardening/no-laser-2d`**: disattivazione temporanea completata e in merge-review.
  - **Prossimo passo**: aggiornare le voci **A1/A2** con lo stato corrente (tap/pan bilanciati + RT60 reale) e consolidare la documentazione QA.

**DoD A** — Balance OK; RT60 ≤ ±5%/banda; densità 50–500 ms↑; no ringing; CPU nel budget HQ.

---

## 2) Fase B — Potenza (Spazio/Dinamica)

**B1. Matrix Engine v2 (16/32/**64**)**
- **WHT Full + Φ-doping**: H64 con permutazioni/flip a passo aureo.
- **GoldenWave (time-sliced)**: 6 stadi butterfly applicati a rotazione con schedulazione aurea → “onde” che si propagano.
- **2D separabile (8×8 ⊗ 8×8)**: alterna righe/colonne per increspature anisotrope.

*Parametri*: `@matrix_dim{16|32|64}`, `@mode_mix{HH|WHT|GoldenWave|2D|Hybrid}`, `@wave_speed`, `@phi_step`, `@perm_depth`, `@coherence`, `@seed`.  
*Guard-rails*: step **unitari** (Δlivello < 0,1 dB); `diag(G)` prima del mixing; `64` riservato a HQ/HYQ.

**B2. Detector range–doppler** — *Range* (RMS lento), *Doppler* (RMS 0.7–3 kHz) → modula `width/hf_damp/modDepth` con slew.

**DoD B** — Densità ≥ HH baseline; movimento senza chorus/pumping; CPU sotto controllo (SIMD, per-blocco).

---

## 3) Fase C — Quantum

**C1. Dither unitario + quantum walk** — rotazioni di Givens (20–40 ms), alternanza U₁/U₂ (seed), energia preservata.  
**C2. Interpolazione “quantistica” tra matrici** — SLERP 2×2 per stadio (percorso unitario), opzionale geodetica approssimata U(N).

**DoD C** — Niente righe di pitch stabili; transizioni tra matrici senza bump > 0,1 dB né swirl.

---

## 4) Fase D — PrimeFear

**D1. PrimeModes** (prime/aureo/plastica/prime_aureo) per ER/FDN + routing M/S differenziato.  
**D2. Guard-rails** — co-primalità linee, min distanza tap ≥ 2 smp, auto-bilanciamento ER.

**DoD D** — Pattern non periodici, M/S musicali, zero artefatti.

---

## 5) Fase E — HYQ GPU (Metal)

**E0. Manager Metal/MPSGraph** — device/queue, buffer Shared, FFT/iFFT batched, triplo buffering, fallback vDSP.  
**E1. LateCon su GPU** — FFT→mul per partizione→iFFT→OLA (prima partizione su CPU/vDSP).  
**E2. GPU Matrix Engine v2** — butterfly 2×2 in compute (per-blocco).

**DoD E** — Throughput GPU > CPU su preset lunghi o molte istanze; equivalenza mag/fase entro tolleranza; fallback indolore.

---

## Parametri (compatti)

**Core**: `@decay @predelay @size @width @mix @color @hf_damp @lf_damp @dispersion`  
**Early**: `@early @erdensity @focus @angle @onsetspread`  
**Spazio/Dinamica**: `@matrix_dim @mode_mix(HH|WHT|GoldenWave|2D|Hybrid) @wave_speed @phi_step @perm_depth @coherence @motion`  
**Quantum**: `@coherence @uwalk_rate @quant_interp @quant_rate @quant_depth`  
**PrimeFear**: `@mode_er @mode_late @mode_mid @mode_side @seed`  
**GPU**: `@gpu(Off|Auto|Force) @fftSize @partSize @queueDepth`

---

## QA & Metriche

- **Early balance**: |L−R| < 0,1 dB.  
- **RT60**: errore ≤ ±5% vs `@decay` per banda.  
- **Densità**: hit-rate 50–500 ms in crescita; zero periodi corti.  
- **No chorus/swirl**: pitch-tracking pulito.  
- **Stereo**: IACC stabile (0–80 ms).  
- **CPU/GPU**: profili per preset; fallback GPU testato.  
- **Factorization**: bit-match ±1 LSB tra pre e post split.  
- **Matrix v2**: Δlivello < 0,1 dB su tutte le modalità/transizioni.

---

## PR plan (piccole, chiare)

`f-factorize` → split • `a-hardening` → A1–A3 • `b-power` → B1–B2 • `c-phonon-hq`/`c-quant-interp` → C1–C2 • `d-primefear` → D1–D2 • `e-gpu-latecon`/`e-gpu-matrix` → E1–E2

**Ogni PR**: changelog, preset aggiornati, test QA verdi. Evitare refactor massivi nelle PR di feature.

---
