# CONTRIBUTING

- **Branch per fase**: `a-hardening`, `b-power`, `c-phonon-hq`, `c-laser`, `d-primefear`, `e-gpu-latecon`.
- **Commit piccoli** con descrizioni chiare; niente refactor massivi nelle PR di feature.
- **PR checklist**: build ok, test QA ok, docs/preset aggiornati, no regressioni CPU.
- **Coding**: C++17, niente allocazioni nel thread audio, SIMD dove utile, buffer allineati.
- **Stile**: SR-aware (tempi in secondi), smoothing params 5–20 ms, trasformazioni energy-preserving.
