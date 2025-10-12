# QA & Metriche

## Test rapidi (accettazione)

- **Early balance (mono impulso)**  
  |L−R| < **0.1 dB** su sweep `@early/@width`.

- **RT60 reale**  
  Misura log-slope per banda → errore ≤ **±5%** rispetto a `@decay`.

- **Densità**  
  Hit-rate 50–500 ms crescente; nessuna periodicità corta (co-primalità linee OK).

- **No chorus/swirl**  
  Nessuna riga stabile nel pitch-tracking della coda.

- **Stereo**  
  IACC stabile 0–80 ms.

- **CPU/GPU**  
  Profilo per preset; nessun glitch; **fallback GPU** verificato.

- **Factorization**  
  A/B bit-match (±1 LSB) tra pre e post `f-factorize` su IR standard.

---

## Test specifici (Hardening & 2D/WHT)

- **Clamp damping dinamico ≥ baseline**  
  Con motion attivo, i valori dinamici MF/HF **non scendono mai** sotto `decayState` (baseline).  
  _Motivo: aveva causato il fail “Dynamic damping dropped below baseline” in QA._ :contentReference[oaicite:2]{index=2}

- **Side-impulse width-modulation**
  Impulso **puro side** (L=+1, R=−1), width modulata lentamente in [0..2] (`test_side_impulse_width_balance`).
  **Energia L/R bilanciata** entro **±3%** sull’intera finestra (base M/S ortonormale + side-injection dedicata).

- **Motion/width (settling)**
  Per verificare crescita/ritorno della larghezza, usa **128 blocchi** di assestamento a ciascun livello di ampiezza (range/doppler con slew).

- **Unità dei pesi di uscita**
  `∑wL² = ∑wR² = 1` dopo `apply_width`; tolleranza numerica ≤ 1e−9.

- **Diffusione WHT2D**
  `test_mix_modes` copre `householder/wht/wht2d/hybrid`: tutte le modalità preservano energia (±1e−9) e forniscono risposte distinte.

- **Mix @100% senza collasso**
  `test_wet_tail_makeup_balance` esegue un impulso lato wet con mix ∈ {0.25, 0.5, 0.75, 1.0} e verifica che la coda RMS a mix=1.0 non cada oltre **−1 dB** rispetto a 0.75.

---

## Preset di prova
- **Hall HQ** • **Plate Dense** • **Chamber HYQ** • **Ambient Infinite** (in `presets/`).

---

## Definition of Done (per PR)
- Tutti i test di cui sopra **OK**.
- `docs/BUGLOG.md` aggiornato per ogni bugfix (Sintomo → Root cause → Fix → Commit SHA → Verifica).
- Documentazione/preset aggiornati; nessuna regressione CPU.
