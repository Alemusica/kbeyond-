# Roadmap di sviluppo — kbeyond~

> **Obiettivo:** portare `kbeyond~` a qualità “flagship” senza stravolgere il core, e rendere il codice **modulare** (non più un megafile), così che umani e agent (Codex) possano lavorare veloci e in sicurezza.

---

## 0) FACTORIZATION-FIRST (prima di tutto)

### 0.1 Principi
- **Equivalenza funzionale:** la fattorizzazione **non deve cambiare il suono** (bit-identico entro ±1 LSB su test di riferimento).
- **Zero logica nuova:** spostiamo codice, nomi, include; **nessuna feature** in questa fase.
- **Boundary chiari:** ogni modulo ha header, dipendenze minime e **invarianti documentati**.

### 0.2 Mappa moduli (target)
source/
├─ kbeyond_wrapper_max.cpp # Adapter Max (I/O, param gateway)
├─ dsp/
│ ├─ params.h # enum/IDs, range, smoothing
│ ├─ buffers.h # ring & delay lines (aligned)
│ ├─ filters.{h,cpp} # 1-pole HP/LP, tilt/color, allpass TPT
│ ├─ mod.{h,cpp} # LFO ricorrenti, jitter, seed
│ ├─ early.{h,cpp} # taps/cluster, pan equal-power, M/S ortonormale
│ ├─ fdn.{h,cpp} # FDN 16×16, Householder(φ), hook WHT/Hybrid
│ ├─ mixing.{h,cpp} # Householder/WHT, rotazioni di Givens (unitary)
│ ├─ decay.{h,cpp} # RT60 per-linea, damping LF/MF/HF
│ ├─ phonon.{h,cpp} # (stub) scattering tra micro-modi
│ ├─ primemodes.{h,cpp} # (stub) prime/aureo/plastica/prime_aureo
│ └─ detector.{h,cpp} # (stub) range/doppler
├─ gpu/
│ ├─ latecon.metal # (stub) kernel pointwise
│ └─ latecon.mm # (stub) MPSGraph FFT/iFFT + OLA
└─ util/
└─ dsp_assert.h # invariants, FTZ/DAZ helpers

markdown
Copy code

### 0.3 Sequenza di fattorizzazione (PR “meccaniche”) — `f-factorize`
1. Estrarre **`params.h`** (IDs, smoothing) e **`buffers.h`** (delay/ring).
2. Spostare filtri 1-pole/allpass in **`filters.{h,cpp}`**.
3. Estrarre **`early.{h,cpp}`** (generazione taps + pan equal-power; logica invariata).
4. Estrarre **`fdn.{h,cpp}`** (rete, Householder φ, loop process).
5. Creare **`mixing.{h,cpp}`** (oggi solo Householder; WHT rimane hook).
6. Creare **`decay.{h,cpp}`** (formula RT60 **pronta ma disattivata** via flag interno).

**DoD F0**
- Test A/B con IR di 30 s: differenza RMS < −100 dBFS; CPU ~ invariata; build/pack passano.  
- *Nota per Codex:* mappa i vecchi simboli → nuovi file; non cambiare i nomi pubblici dell’external Max. Inserisci `// TODO(feature)` dove metteremo le funzioni nuove (C–E).

---

## 1) Fase A — Hardening (bugfix + controllabilità)

**A1. Early bilanciata**
- Taps ER **specchiati** L/R; **pan equal-power**; `width` applicato **in M/S ortonormale**.  
- **Test:** impulso mono → sweep `@early/@width/@phiweight` → |L−R| < **0,1 dB**.

**A2. RT60 reale + Material (3 bande)**
- Nuovo `@decay (s)`: per linea `i` con ritardo `d_i` → `g_i = 10^(−3·d_i/@decay)` (**prima** del mixing).  
- Damping **LF/MF/HF** (1-pole/linea) mappato su un controllo “Material”.  
- `@regen` resta un **trim** di densità, non l’RT60.

**A3. Vitalità a bassa CPU**
- 2–4 **all-pass diffusori** con **jitter lento** (random-walk band-limited).  
- **Micro-jitter** ±1–2 campioni su pochi taps late (seed stabile).

**DoD A**
- Balance OK; errore **RT60 ≤ ±5%** (per banda); densità IR 50–500 ms in crescita; niente ringing corto; CPU nel budget HQ.

---

## 2) Fase B — Potenza (Spazio/Dinamica)

**B1. Matrix Engine v2 (scalabile 16/32/**64**)**
- Motore di mixing unitario modulare con tre modalità ortonormali (tutte **energy-preserving**):
  - **B1.a — WHT “Full + Φ-doping” (64×64)**  
    Applica `H64` e, per evitare staticità, usa **permutazioni circolanti** e **flip di segno** guidati da passo **aureo** (φ).  
    *Effetto:* densità massima con drift energetico “organico”.
  - **B1.b — WHT “GoldenWave” time-sliced**  
    Scomponi `H64` in **6 stadi butterfly** e applicali **a rotazione** (1 stadio per campione/blocco) con **schedulazione aurea**; l’energia “si propaga” a fronti come onde.  
    *Effetto:* movimento continuo e non periodico.
  - **B1.c — WHT 2D separabile (8×8 ⊗ 8×8)**  
    Mappa 64 linee su una griglia 8×8 e alterna trasformate su **righe** e **colonne** (scheduler aureo).  
    *Effetto:* increspature **anisotrope**, forte percezione spaziale “tipo pozza d’acqua”.

  **Parametri nuovi**
  - `@matrix_dim{16|32|64}`, `@mode_mix{HH|WHT|GoldenWave|2D|Hybrid}`,  
    `@wave_speed` (stadi/s o campioni/stadio), `@phi_step` (on/off), `@perm_depth`,  
    `@coherence` (deviazione angolare nei 2×2), `@seed`.

  **Guard-rails**
  - Tutti gli step sono **ortonormali** (variazione livello < 0,1 dB).  
  - `diag(G)` (RT60) resta **prima** del mixing.  
  - Per **Live** limita `@matrix_dim ≤ 32`; `64` riservato a HQ/HYQ.

**B2. Detector range–doppler**
- *Range* = RMS lento (50–150 ms); *Doppler* = RMS veloce 0.7–3 kHz (10–20 ms).  
- Con slew, modula `width`, `hf_damp`, `modDepth` (macro `@motion`).

**B3. 3D1 / FOA (opzionale, CPU)**
- Bus **FOA interno (W,X,Y[,Z])** come modalità opzionale: il Matrix Engine v2 opera per-canale FOA con **render stereo** energy-preserving.  
- Parametri: `@space_mode{Stereo|FOA}`, `@elevation`.

**DoD B**
- Mix unitarî (|ΔdB| < 0,1) su sweep `@mode_mix/@wave_speed/@phi_step`.  
- Densità 50–500 ms ≥ HH baseline; niente chorus/pumping.  
- CPU: con N=64 restare nel budget HQ (SIMD, per-blocco).

---

## 3) Fase C — Quantum / Laser

**C1. Dither unitario + Quantum walk**
- **Rotazioni di Givens** su stato (linee o micro-modi) **ogni 20–40 ms**; energia conservata.  
- Alternanza **U₁/U₂** (quantum walk) con seed; controlli: `@coherence`, `@uwalk_rate`.

**C2. Laser (coerenza HF)**
- Early **mode-locked** con **chirp LFM** o m-sequence (cluster **non densi**, “modali”).  
- **Q-switch**: aumento temporaneo di diffusione nei primi **200–600 ms** della coda.

**C3. Interpolazione “quantistica” tra matrici (unitary path)**
- Transizioni **energy-preserving** tra stati di mixing (es. HH ↔ WHT, GoldenWave ↔ 2D):  
  - **Per-stadio SLERP 2×2** (interpola gli angoli dei butterfly), sempre ortonormale.  
  - (HQ) **Geodetica approssimata** su U(N) componendo micro-rotazioni.  
- Parametri: `@quant_interp` (0–1), `@quant_rate` (Hz), `@quant_depth`.

**DoD C**
- Nessuna riga di pitch stabile; HF “setosa” opzionale; variazione spettrale non stazionaria (distanza L2 mediana fra finestre successive > soglia).  
- **C3:** crossfade tra matrici senza bump > 0,1 dB né swirl.

---

## 4) Fase D — PrimeFear (pattern matematici)

**D1. PrimeModes (nuovo layer)**
- GeneratorI per ER/FDN: **prime**, **aureo** (Beatty φ/φ²), **plastica** (ρ), **prime_aureo**.  
- Parametri: `@mode_er @mode_late @mode_mid @mode_side @seed`.  
- **Routing M/S** differenziato (Mid/Side con pattern diversi).

**D2. Guard-rails**
- **Co-primalità** delle linee FDN; distanza minima tra taps ≥ **2 smp**; **auto-bilanciamento** ER.

**DoD D**
- Pattern non periodici, combinazioni M/S musicali, zero artefatti.

---

## 5) Fase E — HYQ GPU (Metal)

**E0. Manager Metal/MPSGraph**
- Device/queue, buffer **Shared**, batch **FFT/iFFT** (MPSGraph); **triplo buffering**; fallback vDSP.

**E1. LateCon partizionata su GPU (+ PHONON/Laser in freq)**
- Pipeline: **FFT → mul per partizione → iFFT → OLA**.  
- Prima partizione piccola **CPU/vDSP** (maschera latenza).  
- **PHONON-HYQ:** micro-rotazioni **per-bin/per-banda**; **Laser** HF ±0.2–0.5 dB per 50–150 ms.

**E2. GPU Matrix Engine v2 (opzionale, HQ/HYQ)**
- Porting dei **butterfly 2×2** del Matrix Engine v2 in **compute Metal** (applicazione *per-blocco* sugli stati delle linee).  
- Mantieni path **CPU** come fallback (switch `@gpu Auto/Off/Force`).

**E3. 3D1 / FOA su GPU (opzionale)**
- Rotazioni FOA (3×3/4×4) e **beamforming** in compute, *per-blocco* (quasi gratis se E2 è attivo).

**DoD E**
- Throughput GPU > CPU su preset lunghi o molte istanze; equivalenza mag/fase entro tolleranza; fallback indolore.

---

## 6) Parametri (compatti)

- **Core:** `@decay @predelay @size @width @mix @color @hf_damp @lf_damp @dispersion`  
- **Early:** `@early @erdensity @focus @angle @onsetspread`  
- **Spazio/Dinamica:** `@matrix_dim @mode_mix(HH|WHT|GoldenWave|2D|Hybrid) @wave_speed @phi_step @perm_depth @coherence @motion`  
- **Quantum/Laser:** `@coherence @uwalk_rate @laser_mode(Off|Qswitch|ModeLocked) @quant_interp @quant_rate @quant_depth`  
- **PHONON:** `@phonon(Off|Light|HQ|HYQ) @phonon_temp @phonon_mfp @phonon_coh @phonon_disp @phonon_defects`  
- **PrimeFear:** `@mode_er @mode_late @mode_mid @mode_side @seed`  
- **GPU:** `@gpu(Off|Auto|Force) @fftSize @partSize @queueDepth`  
- **FOA (opz.):** `@space_mode(Stereo|FOA) @elevation`

---

## 7) QA & Metriche

- **Early balance:** impulso mono; sweep `@early/@width` → |L−R| < **0,1 dB**.  
- **RT60:** log-slope per banda → errore ≤ **±5%** vs `@decay`.  
- **Densità:** hit-rate 50–500 ms crescente; nessun periodo corto.  
- **No chorus/swirl:** assenza di righe stabili nel pitch-tracking della coda.  
- **Stereo:** IACC stabile (0–80 ms).  
- **CPU/GPU:** profili per preset; fallback GPU verificato.  
- **Factorization:** A/B bit-match (±1 LSB) tra pre e post `f-factorize` su IR standard.  
- **Matrix v2:** variazione livello < 0,1 dB su tutte le modalità e transizioni (C3).

---

## 8) Pianificazione PR (piccole, chiare)

- **`f-factorize`** → Split meccanico in moduli (sez. 0.3).  
- **`a-hardening`** → A1 + A2 + A3.  
- **`b-power`** → B1(a,b,c) + B2 + (opz.) B3.  
- **`c-phonon-hq`** e **`c-laser`** → C1 + C2; **`c-quant-interp`** → C3.  
- **`d-primefear`** → D1 + D2.  
- **`e-gpu-latecon`** → E0 + E1; **`e-gpu-matrix`** → E2; **`e-gpu-foa`** → E3.

**Per ogni PR:** changelog, preset aggiornati, test QA verdi. Evitare refactor massivi nelle PR di feature.

---

## 9) Istruzioni per Codex (per non “perdersi”)

**Entrypoint:**
- Max wrapper: `source/kbeyond_wrapper_max.cpp`  
- Core process: `source/dsp/fdn.{h,cpp}`, `source/dsp/early.{h,cpp}`  
- Param gateway: `source/dsp/params.h`

**Regole:**
- Non cambiare l’API dell’external; aggiunte **solo** nei moduli dedicati.  
- Aggiornare `docs/QA.md` se cambia una metrica; aggiungere preset in `presets/`.  
- Preferire **commit granulari** (≤ 200 LOC) con descrizione “cosa/perché”.

---

## 10) Preset di rotta (per ascoltare subito)

- **Hall HQ** — `@decay 6.5s`, `@hf_damp medio`, `@mode_mix Hybrid`, `@motion 0.3`, ER **plastica**, Late **prime_aureo**.  
- **Plate Dense** — `@decay 2.2s`, `@dispersion ↑`, `@focus 0.7`, ER **aureo**, Late **prime**.  
- **Chamber HYQ (GPU)** — `@decay 4.0s`, `@gpu Auto`, `@fftSize 4096`, `@coherence 0.7`, `@laser_mode Qswitch`, Mid=**prime**, Side=**aureo**.  
- **Ambient Infinite** — `@decay 13s`, `@uwalk_rate 0.5 Hz`, ER **plastica** scarna.

---

### Appendice A — Invarianti da rispettare
- Householder/WHT e Givens sono **unitari** → energia preservata.  
- `width` si applica **solo** in M/S con matrice **ortonormale**.  
- Tutti i tempi sono in **secondi** (SR-aware) → quantizzati a campioni al `prepare()`.  
- Nessuna **allocazione** nel thread audio; **FTZ/DAZ** attivi; aggiornamenti **per blocco**.

### Appendice B — Rischi & Mitigazioni
- **Fattorizzazione** che altera il suono → A/B IR + test automatici; bloccare la PR se mismatch.  
- **GPU latenza** → prima partizione CPU; fallback se `commandBuffer` sfora.  
- **Pattern Prime/Aureo** troppo densi → guard-rails (spacing ≥ 2 smp, auto-bilanciamento ER).
