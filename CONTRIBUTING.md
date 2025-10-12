# CONTRIBUTING.md — kbeyond~

Grazie per contribuire a **kbeyond~**. Obiettivo: qualità “flagship” con codice modulare e sicuro in real‑time.

## Branch & fasi

* Branch per fase: `f-factorize`, `a-hardening`, `b-power`, `c-phonon-hq`, `c-laser`, `c-quant-interp`, `d-primefear`, `e-gpu-latecon`, `e-gpu-matrix`, `e-gpu-foa`.
* Riferimento: `docs/ROADMAP.md` (Definition of Done per ogni fase).

## PR: dimensione & stile

* Commit piccoli, messaggi chiari (**cosa/perché**).
* Evita refactor massivi nelle PR di feature.
* Usa i template Issue/PR in `.github/ISSUE_TEMPLATE/` e `PULL_REQUEST_TEMPLATE.md`.
* Collega sempre un’issue: `Closes #<n>`.

## Bugfix policy (obbligatoria)

* Apri un’Issue con template **Bug**.
* La PR di fix deve aggiornare `docs/BUGLOG.md` con: **Sintomi → Root cause → Soluzione → File → Commit SHA → Verifica/metriche**.
* Spunta “BUGLOG aggiornato” nella checklist PR.

## Factorization‑first

* PR **meccaniche** di split: nessuna logica nuova.
* Bit‑match su IR di riferimento: **diff RMS < −100 dBFS**.
* Mantieni API pubbliche dell’external invariate.

## Coding guidelines

* **C++17**, zero alloc nel thread audio, **FTZ/DAZ** attivi.
* Aggiornamenti **per blocco** (64–128 frame); LFO di servizio (no sin/cos per‑sample).
* **Energy‑preserving**: M/S ortonormale, Householder/WHT/Givens **unitari**.
* **SR‑aware**: tempi in secondi, quantizzati a campioni nel `prepare()`.
* SIMD dove utile; buffer allineati.

## QA minimo per ogni PR

* Early |L−R| < **0.1 dB** (se ER).
* **RT60** per banda errore ≤ **±5%** (se coda).
* Densità 50–500 ms senza periodi corti.
* **No chorus/swirl** (pitch‑tracking pulito).
* IACC 0–80 ms stabile (stereo).
* Profilo **CPU/GPU** e fallback GPU (se rilevante).
* Docs/preset aggiornati (`docs/QA.md`, `docs/ROADMAP.md` se cambiano metriche, `presets/`).

## Worklog (consigliato)

* Mantieni `docs/WORKLOG.md` per note di sessione (branch, commits, prossimi passi).
* Script utili: `scripts/session_resume.sh` (riassunto commit/PR/issue).

## CI

* La CI macOS deve compilare; allega artefatti (`.mxo`) quando possibile.
* Non introdurre dipendenze non documentate; aggiorna README se servono.

---

