# PrimeModes (prime/aureo/plastica/prime_aureo)

## Scopo
Generatori di pattern per ER e FDN + routing M/S differenziato:
- **prime**: ritardi/tap basati su numeri primi (co-primalità, LCM grande).
- **aureo**: Beatty(φ) / Beatty(φ²) per partizionare senza sovrapposizioni.
- **plastica**: costante plastica ρ (Padovan/Perrin) per distribuzioni a bassa coincidenza.
- **prime_aureo**: intersezione/adiacenza fra prime e Beatty(φ).

## Parametri
`@mode_er @mode_late @mode_mid @mode_side @seed`

## Guard-rails
- Distanza minima tra tap ≥ 2 campioni.
- Auto-bilanciamento stereo su ER (equal-power + specchio).
