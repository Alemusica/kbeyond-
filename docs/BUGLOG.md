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

## Esempio (sostituisci con dati reali)
- **Sintomi**: fallimento test “Dynamic damping dropped below baseline” durante motion.
- **Root cause**: i valori dinamici MF/HF potevano scendere sotto i baseline calcolati in `decayState`.
- **Soluzione**: clamp dei valori dinamici a `max(dynamics, baseline)`; predelay minimo a 1 sample; base M/S ortonormale con vettore `side` dedicato all’iniezione.
- **File toccati**: `source/kbeyond_tilde.cpp`, `source/dsp/fdn.cpp`, `tests/...`
- **Commit**: `<sha8>` (branch `a-hardening/no-laser-2d`) — PR #<n>
- **Verifica**: side-impulse width-mod test verde; motion settling esteso; nessuna regressione CPU.
