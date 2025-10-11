# CONTRIBUTING

Grazie per contribuire a **kbeyond~**. Obiettivo: qualità “flagship” con codice modulare e sicuro in real-time.

## Branch & fasi
- **Branch per fase**: `f-factorize`, `a-hardening`, `b-power`, `c-phonon-hq`, `c-laser`, `c-quant-interp`, `d-primefear`, `e-gpu-latecon`, `e-gpu-matrix`, `e-gpu-foa`.
- Riferimento: `docs/ROADMAP.md` (Definition of Done per ogni fase).

## PR: dimensione & stile
- **Commit perfetti** con messaggi chiari (“cosa/perché”).
- **Nessun refactor massivo se non necessario** dentro PR di feature.
- Usa i **template issue/PR** in `.github/ISSUE_TEMPLATE/` e `PULL_REQUEST_TEMPLATE.md`.
- Collega sempre un’issue: `Closes #<n>`.

## Bugfix policy (obbligatoria)
- Apri un’**issue** con template *Bug*.
- La PR di fix deve aggiornare **`docs/BUGLOG.md`** con:
  - Sintomi, root cause, soluzione, file toccati, commit SHA, verifica/metriche.
- Spunta “**BUGLOG aggiornato**” nella checklist PR.

## Fattorizzazione (factorization-first)
- PR **meccaniche** di split: nessuna logica nuova.
- **Bit-match** su IR di riferimento: differenza RMS < −100 dBFS.
- Mantieni API pubbliche dell’external invariate.

## Coding guidelines
- C++17, **zero alloc** nel thread audio, **FTZ/DAZ** attivi.
- Aggiornamenti **per blocco** (64–128 frame); LFO per ricorrenza (no sin/cos per-sample).
- **Energy-preserving**: M/S ortonormale, Householder/WHT/Givens unitarî.
- **SR-aware**: tempi in secondi, quantizzati a campioni nel `prepare()`.
- **SIMD** dove utile; buffer allineati.

## QA minimo per ogni PR
- Early |L−R| < 0.1 dB (se ER).
- RT60 per banda errore ≤ ±5% (se coda).
- Densità 50–500 ms senza periodi corti.
- No chorus/swirl (pitch-tracking pulito).
- IACC 0–80 ms stabile (stereo).
- Profilo **CPU/GPU** e **fallback GPU** (se rilevante).
- **Docs/preset** aggiornati (`docs/QA.md`, `docs/ROADMAP.md` se cambiano metriche, `presets/`).

## Worklog (facoltativo ma consigliato)
- Mantieni **`docs/WORKLOG.md`** per note di sessione (branch, commits, prossimi passi).
- Script utili: `scripts/session_resume.sh` per riassunto commit/PR/issue.

## CI
- La CI macOS deve compilare; allega artefatti (`.mxo`) quando possibile.
- Non introdurre dipendenze non documentate; aggiorna `README` se servono.

Grazie! Ogni PR che rispetta **QA + BUGLOG** accelera review e merge.
