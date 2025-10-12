# docs/ROADMAP.md — kbeyond~

## Stato attuale (Ottobre 2025)

* **Fase F – Factorization**: completata; moduli DSP separati e bit‑match confermato.
* **Fase A – Hardening**:
  ‑ ✅ **A1** (base M/S ortonormale, side‑injection dedicata, width rinormalizzato) — test `side_impulse_width_balance`, `run_side_width_energy_test`.
  ‑ ✅ **A2** (RT60 reale + damping LF/MF/HF con clamp dinamico) — QA “Dynamic damping ≥ baseline” chiuso.
  ‑ ⏳ **A3** (vitalità low‑CPU) — allpass/jitter lenti da progettare.
* **Fasi B–E**: pianificate (Matrix v2, Detector range–doppler, Quantum, PrimeFear, GPU LateCon/Matrix).

### 0) FACTORIZATION‑FIRST

Principi: equivalenza funzionale (bit‑identico ±1 LSB), zero logica nuova, boundary chiari.
Sequenza split: `params/buffers` → `filters` → `early` → `fdn` → `mixing` → `decay` (RT60 pronto ma disattivato).
**DoD F0**: IR 30 s → diff RMS < −100 dBFS; CPU ~ invariata; build ok.

### 1) Fase A — Hardening

* **A1. Early bilanciata** — tap ER specchiati L/R; pan equal‑power; `width` M/S ortonormale.
* **A2. RT60 reale + Material (3 bande)** — `@decay (s)` con `g_i = 10^(−3·d_i/@decay)` + 1‑pole LF/MF/HF.
* **A3. Vitalità low‑CPU** — 2–4 allpass con jitter lento + micro‑jitter su pochi tap.
  **Nota** — **Laser OFF** (temporaneo, `#if KBEYOND_ENABLE_LASER`).

**DoD A** — Balance OK; RT60 ≤ ±5%/banda; densità 50–500 ms↑; no ringing; CPU nel budget HQ.

### 2) Fase B — Potenza

* **B1. Matrix Engine v2 (16/32/64)** — WHT Full + Φ‑doping; **GoldenWave**; **2D** separabile (8×8 ⊗ 8×8).
* **B2. Detector range–doppler** — Range (RMS lento) + Doppler (0.7–3 kHz) → mod `width/hf_damp/modDepth`.

**DoD B** — Densità ≥ HH baseline; movimento senza chorus/pumping; CPU sotto controllo (SIMD, per‑blocco).

### 3) Fase C — Quantum

* **C1. Dither unitario + quantum walk** (Givens 20–40 ms).
* **C2. Interpolazione quantistica (U(N))** — SLERP 2×2 per stadio.

### 4) Fase D — PrimeFear

* **D1. PrimeModes** per ER/FDN + routing M/S.
* **D2. Guard‑rails**: co‑primalità, min distanza tap ≥ 2 smp, auto‑bilanciamento ER.

### 5) Fase E — HYQ GPU (Metal)

* **E0. Manager Metal/MPSGraph** + fallback vDSP.
* **E1. LateCon** su GPU (prima partizione CPU/vDSP).
* **E2. GPU Matrix Engine v2** (butterfly 2×2 compute).

## QA & Metriche

* Early balance: |L−R| < 0.1 dB.
* RT60: errore ≤ ±5% per banda (metodo Schroeder + fit [−5, −35] dB).
* Densità: hit‑rate 50–500 ms in crescita; no periodi corti.
* No chorus/swirl; IACC 0–80 ms stabile.
* CPU/GPU profili; fallback GPU testato.
* Matrix v2: Δlivello < 0.1 dB.

---

