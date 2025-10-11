# PHONON (fononi)

**Obiettivo**: coda viva senza chorus tramite scambi di energia fra micro-modi (unitario).

## Idea
- K micro-modi (banda ERB) + rotazioni di Givens a blocchi (20–40 ms), energia costante.
- Parametri: `@phonon(Off/Light/HQ/HYQ)`, `@phonon_temp`, `@phonon_mfp`, `@phonon_coh`, `@phonon_disp`, `@phonon_defects`.

## Posizionamento
- CPU-HQ: vicino all’FDN (late nativa).
- HYQ: per-bin/per-banda dopo la convoluzione su GPU.

## Scheduler consigliato
- Light: m = 2–4, ε ≤ 0.02
- HQ:    m = 6–10, ε ≤ 0.04
- HYQ:   m = 12–16, ε ≤ 0.06
