Roadmap di sviluppo — kbeyond~

Obiettivo: portare kbeyond~ a qualità “flagship” senza stravolgere il core, e rendere il codice modulare (non più un megafile), così che umani e agent (Codex) possano lavorare veloci e in sicurezza.

0) FACTORIZATION‑FIRST (prima di tutto)
0.1 Principi

Equivalenza funzionale: la fattorizzazione non deve cambiare il suono (bit‑identico entro ±1 LSB su test di riferimento).

Zero logica nuova: spostiamo codice, nomi, include; niente feature.

Boundary chiari: ogni modulo ha header unico, dipendenze minime, e invarianti documentati.

0.2 Mappa moduli (target)
source/
 ├─ kbeyond_wrapper_max.cpp      # Adapter Max (I/O, param gateway)
 ├─ dsp/
 │   ├─ params.h                 # enum/IDs, range, smoothing
 │   ├─ buffers.h                # ring, delay lines (aligned)
 │   ├─ filters.{h,cpp}          # 1-pole HP/LP, tilt/color, allpass TPT
 │   ├─ mod.{h,cpp}              # LFO ricorrenti, jitter, seed
 │   ├─ early.{h,cpp}            # taps/cluster, pan equal-power, M/S ortonormale
 │   ├─ fdn.{h,cpp}              # FDN 16×16, Householder(φ), WHT/Hybrid (stub)
 │   ├─ mixing.{h,cpp}           # Householder/WHT, Givens (unitary ops)
 │   ├─ decay.{h,cpp}            # RT60 per-linea, damping LF/MF/HF
 │   ├─ phonon.{h,cpp}           # (stub) scattering unitario su micro-modi
 │   ├─ primemodes.{h,cpp}       # (stub) prime/aureo/plastica/prime_aureo
 │   └─ detector.{h,cpp}         # (stub) range/doppler
 ├─ gpu/
 │   ├─ latecon.metal            # (stub) kernel pointwise
 │   └─ latecon.mm               # (stub) MPSGraph FFT/iFFT + OLA
 └─ util/
     └─ dsp_assert.h             # invariants, FTZ/DAZ helpers

0.3 Sequenza di fattorizzazione (PR “meccaniche”)

F0. Split iniziale (branch f-factorize)

Estrarre params.h (IDs e smoothing) e buffers.h (delay/ring) dal megafile.

Spostare filtri 1‑pole/allpass in filters.{h,cpp}.

Estrarre early.{h,cpp}: generazione taps + pan equal‑power (senza cambiare logica).

Estrarre fdn.{h,cpp}: rete, feedback matrix (Householder φ), ciclo process.

Creare mixing.{h,cpp} con le sole funzioni Householder (WHT rimane stub).

Creare decay.{h,cpp} con la formula RT60 ma disattivata (flag interno) → solo collegamento.

DoD F0

Test A/B con IR di 30 s: differenza RMS < −100 dBFS; CPU ~ invariata; build/pack passano.

Nota per Codex: mappa i vecchi simboli → nuovi file; non cambiare i nomi pubblici dell’external Max. Inserisci // TODO(feature) dove metteremo le funzioni nuove (C–E).

1) Fase A — Hardening (bugfix + controllabilità)

A1. Early bilanciata

Taps ER specchiati L/R; pan equal‑power; width applicato in M/S ortonormale.

Test: impulso mono → sweep @early/@width/@phiweight → |L−R| < 0,1 dB.

A2. RT60 reale + Material (3 bande)

Nuovo @decay (s): per linea i con ritardo d_i → g_i = 10^(−3·d_i/@decay) (applicato prima del mixing).

Damping LF/MF/HF: filtri 1‑pole per linea; mapping da parametri musicali.

@regen resta solo trim di densità.

A3. Vitalità a bassa CPU

2–4 allpass diffusori con jitter lento (random‑walk band‑limited).

Micro‑jitter ±1–2 campioni su pochi taps late (seed stabile).

DoD A

Balance OK; errore RT60 ≤ ±5% (per banda); densità IR in crescita 50–500 ms; no ringing corto; CPU ≤ budget HQ.

2) Fase B — Potenza (Spazio/Dinamica)

B1. Mixing WHT/Hybrid

Aggiungere Walsh–Hadamard (somme/sottrazioni) e Hybrid (HH↔WHT alternato).

Switch runtime @mode_mix = HH | WHT | Hybrid.

B2. Detector range–doppler

Range = RMS lento (50–150 ms); Doppler = RMS veloce 0.7–3 kHz (10–20 ms).

Mappa con slew → width, hf_damp, modDepth (macro @motion).

DoD B

Spazialità/dinamica percepibile senza chorus/pumping; CPU invariata ~.

3) Fase C — Quantum / Laser

C1. Dither unitario + Quantum walk

Rotazioni di Givens su vettore di stato (linee o micro‑modi) ogni 20–40 ms; energia conservata.

Alternanza U₁/U₂ (quantum walk) con seed; parametri: @coherence, @uwalk_rate.

C2. Laser (coerenza HF)

Early mode‑locked con chirp LFM o m‑sequence (cluster non densi, “modali”).

Q‑switch: aumento temporaneo di diffusione nei primi 200–600 ms della coda.

DoD C

Niente linee di pitch stabili; “setosità” HF attivabile; variazione spettrale temporalmente non stazionaria (distanza L2 mediana fra finestre successive > soglia).

4) Fase D — PrimeFear (pattern matematici)

D1. PrimeModes (nuovo layer)

GeneratorI per ER/FDN: prime, aureo (Beatty φ/φ²), plastica (ρ), prime_aureo.

Parametri: @mode_er @mode_late @mode_mid @mode_side @seed.

Routing M/S differenziato (Mid/Side con pattern diversi).

D2. Guard‑rails

Co‑primalità linee FDN; distanza minima tra taps ≥ 2 smp; auto‑bilanciamento ER.

DoD D

Pattern non periodici, combinazioni M/S musicali, zero artefatti.

5) Fase E — HYQ GPU (Metal)

E0. Manager Metal/MPSGraph

Device/queue, buffer Shared, batch FFT/iFFT (MPSGraph); triplo buffering; fallback vDSP.

E1. LateCon partizionata su GPU (+ PHONON/Laser in freq)

Pipeline: FFT → mul per partizione → iFFT → OLA. Prima partizione piccola CPU/vDSP (maschera latenza).

PHONON‑HYQ: micro‑rotazioni per‑bin/per‑banda; Laser HF ±0.2–0.5 dB per 50–150 ms.

DoD E

Throughput GPU > CPU su preset lunghi; equivalenza mag/fase entro tolleranza; fallback indolore in caso di lag GPU.

6) Parametri (compatti)

Core: @decay @predelay @size @width @mix @color @hf_damp @lf_damp @dispersion
Early: @early @erdensity @focus @angle @onsetspread
Spazio/Dinamica: @mode_mix(HH/WHT/Hybrid) @motion
Quantum/Laser: @coherence @uwalk_rate @laser_mode(Off/Qswitch/ModeLocked)
PHONON: @phonon(Off/Light/HQ/HYQ) @phonon_temp @phonon_mfp @phonon_coh @phonon_disp @phonon_defects
PrimeFear: @mode_er @mode_late @mode_mid @mode_side @seed
GPU: @gpu(Off/Auto/Force) @fftSize @partSize @queueDepth

7) QA & Metriche

Early balance: impulso mono; sweep @early/@width → |L−R| < 0,1 dB.

RT60: log‑slope per banda → errore ≤ ±5% vs @decay.

Densità: hit‑rate 50–500 ms crescente; nessun periodo corto.

No chorus/swirl: assenza di righe stabili nel pitch‑tracking della coda.

Stereo: IACC stabile (0–80 ms).

CPU/GPU: profili per preset; fallback GPU verificato.

Factorization: A/B bit‑match (entro ±1 LSB) tra pre e post F0 su IR standard.

8) Pianificazione PR (piccole, chiare)

f-factorize: Split meccanico in moduli (0.3).

a-hardening: A1+A2+A3.

b-power: B1+B2.

c-phonon-hq e c-laser: C1+C2.

d-primefear: D1+D2.

e-gpu-latecon: E0+E1.

Ogni PR: changelog, preset aggiornati, test QA verdi. Niente refactor massivi dentro PR di feature.

9) Istruzioni per Codex (per non “perdersi”)

Entry points:

Max wrapper: source/kbeyond_wrapper_max.cpp

Core process: source/dsp/fdn.{h,cpp} e source/dsp/early.{h,cpp}

Param gateway: source/dsp/params.h

Regole:

Non cambiare l’API dell’external; inserire feature solo nei moduli dedicati.

Aggiornare docs/QA.md se cambia una metrica; aggiungere preset in presets/.

Preferire commit granulari (≤200 LOC) con descrizione “cosa/perché”.

10) Preset di rotta (per ascoltare subito)

Hall HQ: @decay 6.5s, @hf_damp medio, @mode_mix Hybrid, @motion 0.3, ER plastica, Late prime_aureo.

Plate Dense: @decay 2.2s, @dispersion ↑, @focus 0.7, ER aureo, Late prime.

Chamber HYQ (GPU): @decay 4.0s, @gpu Auto, @fftSize 4096, @coherence 0.7, @laser_mode Qswitch, Mid=prime, Side=aureo.

Ambient Infinite: @decay 13s, @uwalk_rate 0.5 Hz, ER plastica scarna.

Appendice A — Invarianti da rispettare

Householder/WHT e Givens sono unitari → energia preservata.

width si applica solo in M/S con matrice ortonormale.

Tutti i tempi sono in secondi (SR‑aware) → quantizzati a campioni al prepare().

Nessuna allocazione nel thread audio; FTZ/DAZ attivi; aggiornamenti per blocco.

Appendice B — Rischi & Mitigazioni

Fattorizzazione che altera il suono → A/B IR + test automatici; bloccare la PR se mismatch.

GPU latenza → prima partizione CPU; fallback se commandBuffer sfora.

Pattern Prime/Aureo troppo densi → guard‑rails (spacing ≥2 smp, auto‑bilanciamento ER).
