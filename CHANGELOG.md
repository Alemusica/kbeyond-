# Changelog

## [Unreleased]
- Riassegnato `@early` alla diffusione stereo e `@focus` al livello delle riflessioni precoci, mantenendo il cluster laser legato al focus.
- Bilanciato il mix wet/dry con crossfade equal-power e compensazione automatica del livello wet per evitare cali a 100% mix.
- Ricalibrato il mapping di `@regen`/`@decay` per preservare un floor udibile e una crescita progressiva del feedback.
- Aggiunto test DSP che misura `fdn_decay` e l'energia d'uscita agli step 0.1/0.5/0.9.
- Introdotto detector range–doppler con attributo `@motion` per modulare `@width`, damping MF/HF e `@moddepth` mantenendo trasformazioni energy-preserving.
- Estesa la suite di test con verifiche su width dinamica, mod depth e normalizzazione dell'energia.

# [0.2.0] - 2024-07-01
- Aggiunto attributo `@mode_mix` con modalità Householder, Walsh-Hadamard e Hybrid.
- Implementato diffusore selezionabile con helper dedicati e normalizzazione dell'uscita.
- Aggiornati help patcher, preset e documentazione per descrivere le nuove modalità.
- Estesa la suite di test con verifiche sull'energia delle modalità di mix.

## [0.1.0] - 2024-06-01
- Prima pubblicazione dell'external `kbeyond~` con motore FDN 16×16 e Householder(φ).
- Aggiunti file di build (CMake, Makefile, script macOS/Windows) e package metadata.
- Inclusi help patcher, preset di esempio e documentazione.
