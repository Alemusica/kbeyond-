# AGENTS.md

> Istruzioni operative per agent (es. Codex) su **kbeyond-**.
> Obiettivo: qualità “flagship” senza stravolgere il core FDN, codice modulare e sicuro in real-time.

## 0) Goal & Tier
- **Live (CPU-only)**: early bilanciata, FDN stabile, latenza bassa.
- **Massima Qualità (CPU-HQ)**: RT60 “vero”, damping 3-bande, mixing HH/WHT/Hybrid, evoluzione senza chorus.
- **Hyper Qualità (HYQ, GPU)**: late tail su GPU (Metal/MPSGraph), decorrelazioni “quantum”.

> **Semplificazione attuale**: percorso **Laser disabilitato** (no-op sugli attributi), focus su **matrici 2D/WHT**.

## 1) Guard-rails
- Niente allocazioni nel **thread audio**; **FTZ/DAZ**; aggiornamenti **per blocco**.
- Parametri con **smoothing** 5–20 ms; trasformazioni **energy-preserving** (M/S ortonormale, mixing unitari).
- **Core FDN intatto**; feature modulari e reversibili (nessun refactor invasivo in PR funzionali).

## 2) Roadmap (fasi → branch/PR)
**F — Factorization (`f-factorize`)**
- Split meccanico in moduli (`params.h`, `buffers.h`, `filters.*`, `early.*`, `fdn.*`, `mixing.*`, `decay.*`, `detector.*`, stub `phonon.*`, `primemodes.*`). **DoD:** bit-match entro ±1 LSB e CPU ~ invariata.

**A — Hardening (`a-hardening`)**
1. Early energy-balanced (tap specchiati, pan equal-power, width in M/S ortonormale). Test impulso mono → |L−R| < 0.1 dB.
2. RT60 vero `@decay`: per-linea `g_i = 10^(−3·d_i/@decay)`; damping LF/MF/HF (1-pole).
3. Vitalità low-CPU: 2–4 allpass con jitter lento + micro-jitter su pochi tap.

**B — Potenza (`b-power`)**
4. **Matrix Engine v2** (scalabile 16/32/**64**):  
   WHT “Full + Φ-doping”, **GoldenWave** time-sliced (stadi butterfly a schedulazione aurea), **2D separabile** (H8 ⊗ H8).
5. Detector **range–doppler** (RMS lento + banda 0.7–3 kHz) che modula `width/hf_damp/modDepth` con slew.

**C — Quantum (`c-phonon-hq`, `c-quant-interp`)**
6. **PHONON** (scambi di energia fra micro-modi con rotazioni di Givens, unitario).  
7. **Interpolazione “quantistica”** fra matrici (SLERP 2×2 per stadio, percorso unitario in U(N)).

**D — PrimeFear (`d-primefear`)**
8. **PrimeModes** (prime/aureo/plastica/prime_aureo) per ER/FDN + routing M/S. Guard-rails: co-primalità, distanza ≥ 2 smp, auto-bilanciamento.

**E — HYQ GPU (`e-gpu-latecon`, `e-gpu-matrix`)**
9. Late tail partizionata su GPU (MPSGraph FFT), prima partizione CPU/vDSP, fallback automatico.  
10. **GPU Matrix Engine v2** (butterfly 2×2 in compute, per-blocco).

## 3) Definition of Done (per PR)
- Early |L−R| < 0.1 dB; **RT60** errore ≤ ±5% per banda; nessun ringing corto.
- Niente chorus/swirl (pitch-tracking pulito).
- **CPU:** HQ ≲ 1–1.3%/core @48 kHz (indic.) • **GPU:** throughput > CPU e fallback testato.
- **Docs/preset** aggiornati; **test QA verdi** (vedi `docs/QA.md`).

## 4) Percorsi & Documenti
- Roadmap: `docs/ROADMAP.md`
- QA & metriche: `docs/QA.md`
- Prime modes: `docs/PRIME_MODES.md`
- PHONON: `docs/PHONON.md`

---

## 5) NEXT PR (istruzioni puntuali per Codex)

**Branch**: `a-hardening/no-laser-2d`  
**Titolo PR**: `FDN hardening: laser OFF (temp), 2D/WHT basis, side-injection vec, predelay min, damping clamp`

**Cosa fare**
1) **Laser OFF** (temp.)  
   - Racchiudi il codice laser in `#if KBEYOND_ENABLE_LASER` (default 0).  
   - Gli attributi laser restano **no-op** (API compatibili).
2) **Base M/S & side-injection dedicata**  
   - In `t_kbeyond`: mantieni `outBaseMid/outBaseSide` come *voicing*.  
   - Crea `outWeightsSide` **ortonormale** rispetto a Mid (Gram–Schmidt) solo per **iniezione side** in FDN.  
   - Output L/R generato da `(mid ± width*side)` e **rinormalizzato** (∑w²=1 per canale).
3) **FDN write_feedback**  
   - Firma con vettore **sideWeights** separato (non riusare outL): iniezione = `midIn*inWeights[i] + sideIn*sideWeights[i]*0.15`.
4) **Predelay minimo**  
   - `predSamps >= 1.0` (evita “range=0” nel detector a @predelay=0).
5) **Clamp damping dinamico**  
   - MF/HF **mai sotto baseline**: `dampMFValue = max(dampMFValue, baselineMF)`, idem HF.  
   - Chiude il fail “**Dynamic damping dropped below baseline**”. :contentReference[oaicite:1]{index=1}
6) **Test QA**  
   - Aggiungi **side-impulse width-mod test**: impulso puramente side, width modulata [0..2]; **energia L/R bilanciata** entro **±3%**.  
   - **Motion/width**: “settling” più lungo (es. 128 blocchi brevi) prima di verificare crescita/ritorno.

**DoD PR**
- ✅ Tutti i test verdi (inclusi nuovi)  
- ✅ Nessun bump di livello >0.1 dB fra modalità di mixing  
- ✅ Docs aggiornati (`docs/QA.md`, `docs/ROADMAP.md`)  
- ✅ Voce in `docs/BUGLOG.md` per il damping clamp

**File da toccare (indicativi)**
- `source/kbeyond_tilde.{h,cpp}` (basis M/S, width, predelay min, laser guard)  
- `source/dsp/fdn.{h,cpp}` (firma + uso sideWeights)  
- `tests/*` (nuovo test side-impulse + settling motion)  
- `docs/QA.md`, `docs/ROADMAP.md`, `docs/BUGLOG.md`

---
