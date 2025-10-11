# QA & Metriche

## Test rapidi
- **Early balance**: impulso mono; sweep `@early/@width` → |L−R| < 0.1 dB.
- **RT60**: misura log-slope per banda → errore ≤ ±5% su `@decay`.
- **Densità**: hit-rate 50–500 ms crescente; niente periodicità corte (co-primalità linee).
- **No chorus/swirl**: assenza di righe stabili nel pitch-tracking della coda.
- **Stereo**: IACC stabile (0–80 ms).
- **CPU/GPU**: profilo per preset; fallback GPU indolore.

## Preset di prova
- Hall HQ
- Plate Dense
- Chamber HYQ
- Ambient Infinite (in `presets/`).

## Definition of Done (per PR)
- Tutti i test sopra OK, docs/preset aggiornati, nessuna regressione.
