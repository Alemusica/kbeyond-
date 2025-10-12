# docs/ROADMAP.md — kbeyond~

## Stato attuale (Ottobre 2025)

* **Fase F – Factorization**: completata; moduli DSP separati e bit‑match confermato.
* **Fase A – Hardening**:
  ‑ ✅ **A1** (base M/S ortonormale, side‑injection dedicata, width rinormalizzato) — test `side_impulse_width_balance`, `run_side_width_energy_test`.
  ‑ ✅ **A2** (RT60 reale + damping LF/MF/HF con clamp dinamico) — QA “Dynamic damping ≥ baseline” chiuso.
  ‑ ✅ **A3** **Vitalità HQ** (senza compromessi, orientata ad Apple Silicon).
* **Fasi successive**: priorità **GPU (E)** → poi **Matrix v2/Detector (B)** → quindi **Quantum (C)** e **PrimeFear (D)**.

### 0) FACTORIZATION‑FIRST

Principi: equivalenza funzionale (bit‑identico ±1 LSB), zero logica nuova, boundary chiari.
Sequenza split: `params/buffers` → `filters` → `early` → `fdn` → `mixing` → `decay` (RT60 pronto ma disattivato).
**DoD F0**: IR 30 s → diff RMS < −100 dBFS; CPU ~ invariata; build ok.

### 1) Fase A — Hardening (completata)

* **A1. Early bilanciata** — tap ER specchiati L/R; pan equal‑power; `width` M/S ortonormale.
* **A2. RT60 reale + Material (3 bande)** — `@decay (s)` con `g_i = 10^(−3·d_i/@decay)` + 1‑pole LF/MF/HF.
* **A3. Vitalità HQ** — jitter/all‑pass lenti e micro‑variazioni **senza vincoli di budget CPU**; niente swirl/chorus; mantenere trasformazioni unitarie.
  **Nota** — **Laser OFF** (temporaneo, `#if KBEYOND_ENABLE_LASER`).

**DoD A** — Balance OK; RT60 ≤ ±5%/banda; densità 50–500 ms↑; no ringing.

### 2) Fase E — HYQ GPU (Metal) — **anticipata**

* **E0. Manager Metal/MPSGraph** + fallback vDSP; buffer Shared; triplo buffering.
* **E1. LateCon su GPU** — FFT→mul partizionata→iFFT→OLA (**prima partizione su CPU/vDSP**). Attributi: `@gpu {Off|Auto|Force}`, `@fftSize`, `@partSize`, `@queueDepth`.
* **E2. GPU Matrix Engine v2** — butterfly 2×2 unitari in compute, scheduling per‑blocco (GoldenWave/2D).
* **DoD E** — Parità numerica CPU↔GPU entro soglia; throughput GPU > CPU su preset lunghi/multi‑istanza; fallback indolore.

### 3) Fase B — Potenza (CPU)

* **B1. Matrix Engine v2 (16/32/64)** — WHT Full + Φ‑doping; **GoldenWave**; **2D** separabile (8×8 ⊗ 8×8); SIMD NEON/Accelerate su Apple Silicon.
  * ✅ Modalità `wht2d` (H4 ⊗ H4) disponibile come `@mode_mix` a parità di energia.
* **B2. Detector range–doppler** — Range (RMS lento) + Doppler (0.7–3 kHz) → mod `width/hf_damp/modDepth` con slew **per blocco** (no per‑sample).

**DoD B** — Densità ≥ HH baseline; movimento senza chorus/pumping; Δlivello < 0.1 dB tra modalità.

### 4) Fase C — Quantum

* **C1. Dither unitario + quantum walk** (Givens 20–40 ms).
* **C2. Interpolazione quantistica (U(N))** — SLERP 2×2 per stadio.

### 5) Fase D — PrimeFear

* **D1. PrimeModes** per ER/FDN + routing M/S.
* **D2. Guard‑rails**: co‑primalità, min distanza tap ≥ 2 smp, auto‑bilanciamento ER.

## QA & Metriche

* **Parità CPU↔GPU**: RMS diff ≤ **−60 dBFS** statico, ≤ **−50 dBFS** in transitori (E1/E2).
* **Mixer M/S**: energia conservata; **THD side‑tone < −90 dBFS** a `width ∈ {0.5,1,2}`.
* **RT60**: errore ≤ **±5%** per banda (Schroeder + fit [−5, −35] dB).
* **Densità**: hit‑rate 50–500 ms in crescita; no periodi corti.
* **No chorus/swirl**; IACC 0–80 ms stabile.
* **GPU**: throughput > CPU; fallback GPU testato.
* **Matrix v2**: Δlivello < **0.1 dB** tra modalità/transizioni.

---
