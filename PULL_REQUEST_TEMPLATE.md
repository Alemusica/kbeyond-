## Cosa fa questa PR
<!-- Breve riassunto della modifica e del beneficio sonoro/tecnico -->

## Perché
<!-- Motivazione: bug fix, hardening, performance, nuova feature coerente con ROADMAP -->

## Come è implementata
- Moduli toccati:
- Note su performance/latency:
- Compatibilità preset/params: (breaking changes? sì/no)

## QA (spuntare)
- [ ] **Early balance** |L−R| < 0.1 dB (se toccata ER)
- [ ] **RT60** errore ≤ ±5% (se toccata la coda)
- [ ] **No chorus/swirl** (pitch-tracking coda pulito)
- [ ] **CPU** in target / **GPU** fallback testato (se rilevante)
- [ ] **Docs/preset** aggiornati (`docs/QA.md`, `presets/` se serve)
- [ ] **BUGLOG** aggiornato (`docs/BUGLOG.md`) *se la PR chiude un bug*
- [ ] **Bit-match** post-fattorizzazione (se PR “meccanica” di split)

## Test eseguiti
<!-- Impulso mono, misure RT60 per banda, hit-rate 50–500ms, IACC, profilo CPU/GPU, ecc. -->

## Relazionato a
Closes #...

## Note per il reviewer
<!-- Rischi, rollback, flag sperimentali, TODO follow-up -->
