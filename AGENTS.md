# AGENTS.md — kbeyond~ (Istruzioni operative per assistenti / Codex)

> Obiettivo: qualità **flagship** senza stravolgere il core FDN. Codice modulare e **real‑time safe**.

## 0) Goal & Tier

* **Live (CPU‑only)** — early bilanciata, FDN stabile, latenza bassa.
* **Massima Qualità (CPU‑HQ)** — RT60 “vero”, damping 3‑bande, mixing HH/WHT/Hybrid, evoluzione **senza chorus**.
* **Hyper Qualità (HYQ, GPU)** — late tail su GPU (Metal/MPSGraph), decorrelazioni “quantum”.

> **Semplificazione attuale**: percorso **Laser disabilitato** (no‑op sugli attributi), focus su **matrici 2D/WHT**.

### Flag & Tier mapping (default di build/esecuzione)

| Tier      | Macro/Flag             | Default                                      |
| --------- | ---------------------- | -------------------------------------------- |
| Live      | `KBEYOND_ENABLE_LASER` | **0** (attributi accettati ma no‑op)         |
|           | Mixing base            | **WHT/2D** (unitario; preferire WHT per CPU) |
|           | Motion/Detector        | **Off** di default                           |
| CPU‑HQ    | `KBEYOND_ENABLE_LASER` | 0                                            |
|           | Mixing                 | WHT/Hybrid; modulazioni moderate             |
|           | Motion/Detector        | On (slew lento; clamp Δ per blocco)          |
| HYQ (GPU) | GPU Late‑Con           | On (fallback CPU/vDSP)                       |
|           | Matrix Engine GPU      | pianificato (butterfly 2×2 compute)          |

## 1) Guard‑rails

* **Thread audio**: nessuna allocazione; **FTZ/DAZ**; aggiornamenti **per blocco**.
* **Smoothing**: parametri 5–20 ms; trasformazioni **energy‑preserving** (M/S ortonormale, mixing unitari).
* **Core FDN intatto**; feature **modulari e reversibili** (no refactor invasivi nelle PR funzionali).

## 2) Roadmap (fasi → branch/PR)

**F — Factorization (`f-factorize`)**
Split meccanico in moduli (`params.h`, `buffers.h`, `filters.*`, `early.*`, `fdn.*`, `mixing.*`, `decay.*`, `detector.*`, stub `phonon.*`, `primemodes.*`). **DoD**: bit‑match entro ±1 LSB e CPU ~ invariata.

**A — Hardening (`a-hardening`)**

1. Early energy‑balanced (tap specchiati, pan equal‑power, width in M/S ortonormale). Test impulso mono → |L−R| < **0.1 dB**.
2. RT60 vero `@decay`: per‑linea `g_i = 10^(−3·d_i/@decay)`; damping LF/MF/HF (1‑pole).
3. Vitalità low‑CPU (2–4 allpass/jitter lento) — *non bloccante*.

**B — Potenza (`b-power`)**
4. **Matrix Engine v2** (scalabile 16/32/**64**): WHT “Full + Φ‑doping”, **GoldenWave** time‑sliced (butterfly con schedulazione aurea), **2D separabile** (H8 ⊗ H8).
5. Detector **range–doppler** (RMS lento + banda 0.7–3 kHz) che modula `width/hf_damp/modDepth` con slew.

**C — Quantum (`c-phonon-hq`, `c-quant-interp`)**
6. **PHONON** (scambi di energia fra micro‑modi con rotazioni di Givens, unitario).
7. **Interpolazione “quantistica”** fra matrici (SLERP 2×2 per stadio, percorso unitario in U(N)).

**D — PrimeFear (`d-primefear`)**
8. **PrimeModes** (prime/aureo/plastica/prime_aureo) per ER/FDN + routing M/S. Guard‑rails: co‑primalità, distanza ≥ 2 smp, auto‑bilanciamento.

**E — HYQ GPU (`e-gpu-latecon`, `e-gpu-matrix`)**
9. Late tail partizionata su GPU (MPSGraph FFT), prima partizione CPU/vDSP, fallback automatico.
10. **GPU Matrix Engine v2** (butterfly 2×2 in compute, per‑blocco).

## 3) Definition of Done (per PR)

* Early |L−R| < **0.1 dB**; **RT60** errore ≤ **±5%** per banda; nessun ringing corto.
* No **chorus/swirl** (pitch‑tracking pulito).
* **CPU:** HQ ≲ 1–1.3%/core @48 kHz (indic.) • **GPU:** throughput > CPU e fallback testato.
* **Docs/preset** aggiornati; **test QA verdi** (vedi `docs/QA.md`).

## 4) Percorsi & Documenti

* Roadmap: `docs/ROADMAP.md`
* QA & metriche: `docs/QA.md`
* Prime modes: `docs/PRIME_MODES.md`
* PHONON: `docs/PHONON.md`

---

## 5) NEXT PR (istruzioni puntuali per Codex)

**Branch**: `a-hardening/no-laser-2d`
**Titolo PR**: `FDN hardening: laser OFF (temp), 2D/WHT basis, side-injection vec, predelay min, damping clamp`

### Cosa fare

1. **Laser OFF** (temp.)
   ‑ Racchiudi il codice laser in `#if KBEYOND_ENABLE_LASER` (default **0**).
   ‑ Gli attributi laser restano **no‑op** (API compatibili).
2. **Base M/S & side‑injection dedicata**
   ‑ In `t_kbeyond`: mantieni `outBaseMid/outBaseSide` come *voicing*.
   ‑ Crea `outWeightsSide` **ortonormale** rispetto a Mid (Gram–Schmidt) **solo** per **iniezione side** in FDN.
   ‑ Output L/R generato da `(mid ± width*side)` e **rinormalizzato** (∑w² = 1 per canale).
3. **FDN write_feedback**
   ‑ Firma con vettore **sideWeights** separato (non riusare outL): iniezione = `midIn*inWeights[i] + sideIn*sideWeights[i]*0.15`.
4. **Predelay minimo**
   ‑ `predSamps >= 1.0` (evita “range=0” nel detector a @predelay=0).
5. **Clamp damping dinamico**
   ‑ MF/HF **mai sotto baseline**: `dampMFValue = max(dampMFValue, baselineMF)`, idem HF.
   ‑ Chiude il fail “Dynamic damping dropped below baseline”.
6. **Test QA**
   ‑ Aggiungi **side‑impulse width‑mod test**: impulso puramente side, width modulata [0..2]; **energia L/R bilanciata** entro **±3%**.
   ‑ **Motion/width**: “settling” più lungo (p.es. **128 blocchi brevi**) prima di verificare crescita/ritorno.
   ‑ (**Opzionale** per questa PR) **Leak‑compensation** del mixer M/S **spostata a block‑rate** o disattivata; QA **THD side‑tone** < −90 dBFS.

### DoD PR

* ✅ Tutti i test verdi (inclusi nuovi)
* ✅ Nessun bump di livello > **0.1 dB** fra modalità di mixing
* ✅ Docs aggiornati (`docs/QA.md`, `docs/ROADMAP.md`)
* ✅ Voce in `docs/BUGLOG.md` per il damping clamp

### File da toccare (indicativi)

* `source/kbeyond_tilde.{h,cpp}` (basis M/S, width, predelay min, laser guard)
* `source/dsp/fdn.{h,cpp}` (firma + uso sideWeights)
* `tests/*` (nuovo test side‑impulse + settling motion + THD side‑tone)
* `docs/QA.md`, `docs/ROADMAP.md`, `docs/BUGLOG.md`

---

## 6) Invarianti & Glossario (per evitare ambiguità)

* **M/S ortonormale:** `mid = v/||v||`; `side = (w − (w·mid)mid)/||...||` — poi `L=(mid + w*side)/√(1+w²)`, `R=(mid − w*side)/√(1+w²)`.
* **Unitario/energy‑preserving:** tutte le trasformazioni (WHT/Householder/Hybrid) devono preservare l’energia entro tolleranza numerica.
* **Stabilità FDN:** feedback/gain per‑linea < 1; nessuna periodicità corta; co‑primalità delle lunghezze.
* **RT‑safe:** nessuna alloc; smoothing parametri per blocco; evitare sin/cos per‑sample costosi.

## 7) Misura RT60 (metodo di riferimento)

1. **IR** con impulso; **per‑banda** (LF/MF/HF) con filtri 1‑pole come in codice.
2. **Schroeder** reverse integration → curva in dB.
3. **Fit lineare** su [−5, −35] dB (o intervallo stabile) → RT60 atteso.
4. Errore relativo per banda ≤ **±5%**.

## 8) Mini‑checklist PR (copia/incolla)

```markdown
- [ ] RT‑safe: nessuna alloc/lock; FTZ/DAZ attivi.
- [ ] Base M/S ortonormale; ∑wL²=∑wR²=1.
- [ ] Predelay min ≥ 1 campione.
- [ ] Clamp dinamico damping ≥ baseline.
- [ ] Side‑impulse width sweep: energia L/R ±3%.
- [ ] RT60 per banda entro ±5%.
- [ ] Nessun chorus/swirl; check pitch‑tracking.
- [ ] Docs (QA/ROADMAP/BUGLOG) aggiornati.
```

---


