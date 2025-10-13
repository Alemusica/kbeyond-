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

## [2025-10-14] #002 — Mix L/R statico su base ortonormale
- **Sintomi**: la coda wet sfarfallava sui preset wide; `mix_mid_side_to_lr` modulava a frequenza audio (ratio mid/side) con collassi istantanei e differenze marcate rispetto a kBeyond/Airwindows.
- **Root cause**: i coefficienti L/R venivano ricombinati per-sample tramite `ratio = absSide/absMid`, ignorando i pesi normalizzati generati in `apply_width`; quando `mid` si avvicinava a zero il mix diventava instabile.
- **Soluzione**: proiezione L/R statica: `apply_width` usa `(mid ± width*side)` rinormalizzato, salva i dot product Mid/Side↦L/R e `mix_mid_side_to_lr` li riutilizza senza dipendenze da ampiezze istantanee.
- **File toccati**: `source/kbeyond_tilde.cpp`, `tests/motion_detector_test.cpp`, `docs/QA.md`, `docs/ROADMAP.md`.
- **Commit**: (pending, branch `a-hardening/no-laser-2d`) — PR #TBD
- **Verifica**: `run_mid_side_mix_normalization`, `test_side_impulse_width_balance`, `run_side_width_energy_test`.
- **Note/rollback**: reintrodurre il mix dinamico ripristina il glitch stereo e il collapse della coda.

## [2025-10-13] #001 — Dynamic damping ≥ baseline
- **Sintomi**: QA "Dynamic damping dropped below baseline" falliva quando il detector range/doppler riduceva l'energia wet, causando valori MF/HF < baseline.
- **Root cause**: il calcolo dinamico di `dampMFValue/dampHFValue` non bloccava esplicitamente il minimo rispetto ai coefficienti calcolati in `update_decay`.
- **Soluzione**: applicato clamp con `std::max` verso `decayState.dampMF/dampHF` direttamente nel perform, assicurando che l'automation dinamica non scenda mai sotto la baseline.
- **File toccati**: `source/kbeyond_tilde.cpp`, `tests/early_reflections_test.cpp`, `docs/QA.md`, `docs/ROADMAP.md`.
- **Commit**: (pending, branch `a-hardening/no-laser-2d`) — PR #TBD
- **Verifica**: `motion_tests::run_motion_moddepth_response` (QA) + `test_side_impulse_width_balance` (width 0..2, energia bilanciata ±3%).
- **Note/rollback**: per ripristinare il comportamento precedente rimuovere il clamp nel perform, ma riemerge il fail QA.

## [2025-10-12] #000 — Wet mix makeup @100%
- **Sintomi**: spingendo `@mix` verso 1.0 il livello uscita crollava di parecchi dB, rendendo inascoltabile il full-wet e rompendo i preset wide.
- **Root cause**: a `mix=1.0` il codice azzerava la compensazione (`wetMakeup`), lasciando la coda senza make-up mentre il dry veniva rimosso.
- **Soluzione**: applicato `wetMakeup` su tutto il range (0<@mix≤1) con clamp ≤8× e QA che misura la coda RMS fino a full-wet.
- **File toccati**: `source/kbeyond_tilde.cpp`, `tests/early_reflections_test.cpp`, `docs/QA.md`, `docs/ROADMAP.md`.
- **Commit**: (pending, branch `a-hardening/no-laser-2d`) — PR #TBD
- **Verifica**: `test_wet_tail_makeup_balance` (KBEYOND_BUILD_TESTS=ON, richiede `MAX_SDK_ROOT`).
- **Note/rollback**: ripristinare il commit precedente se servisse tornare alla curva originale del mix; nessun preset da migrare.
