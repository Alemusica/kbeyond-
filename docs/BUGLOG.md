# BUGLOG

> Registro ufficiale dei bugfix e delle soluzioni. Una voce per **ogni** issue/PR chiusa.
> Ordine: dal più recente al più vecchio. Mantieni le voci sintetiche e verificabili.

## Template (copia e incolla)
```markdown
## [YYYY-MM-DD] #<issue> — Titolo breve
- **Sintomi**: cosa succedeva (riassunto + link all’issue).
- **Root cause**: causa individuata (file/linea, logica).
- **Soluzione**: cosa è stato cambiato e perché funziona.
- **File toccati**: `source/...`, `docs/...` (lista corta).
- **Commit**: `<sha8>` (branch) — PR #<n>
- **Verifica**: test effettuati (es. impulso mono |L−R|<0.1 dB; RT60 ±5%).
- **Note/rollback**: eventuali rischi o come tornare indietro.
```

---

## [2024-05-24] #N/A — Dynamic damping dropped below baseline
- **Sintomi**: i test QA registravano il fallimento “Dynamic damping dropped below baseline” durante cicli di motion con width modulata.
- **Root cause**: mancavano i clamp per MF/HF damping dinamici, permettendo ai valori calcolati di scendere sotto il baseline derivato da `decayState`.
- **Soluzione**: introdotta macro di disattivazione laser + clamp `max(dynamics, baseline)` su MF/HF; predelay forzato ≥1 sample per stabilizzare il detector.
- **File toccati**: `source/kbeyond_tilde.cpp`, `source/dsp/fdn.cpp`, `tests/qa_side_impulse_width.cpp`, `tests/qa_motion_settling.cpp`, `docs/QA.md`, `docs/ROADMAP.md`.
- **Commit**: pending (branch `a-hardening/no-laser-2d`) — PR #TBD.
- **Verifica**: suite QA con test `qa_side_impulse_width` e `qa_motion_settling` verde; regressione CPU negativa.
## [2025-10-12] #000 — Wet mix makeup @100%
- **Sintomi**: spingendo `@mix` verso 1.0 il livello uscita crollava di parecchi dB, rendendo inascoltabile il full-wet e rompendo i preset wide.
- **Root cause**: a `mix=1.0` il codice azzerava la compensazione (`wetMakeup`), lasciando la coda senza make-up mentre il dry veniva rimosso.
- **Soluzione**: applicato `wetMakeup` su tutto il range (0<@mix≤1) con clamp ≤8× e QA che misura la coda RMS fino a full-wet.
- **File toccati**: `source/kbeyond_tilde.cpp`, `tests/early_reflections_test.cpp`, `docs/QA.md`, `docs/ROADMAP.md`.
- **Commit**: (pending, branch `a-hardening/no-laser-2d`) — PR #TBD
- **Verifica**: `test_wet_tail_makeup_balance` (KBEYOND_BUILD_TESTS=ON, richiede `MAX_SDK_ROOT`).
- **Note/rollback**: ripristinare il commit precedente se servisse tornare alla curva originale del mix; nessun preset da migrare.
