# README.md — kbeyond~ (Max/MSP external)

`kbeyond~` è un riverbero **FDN 16×16** con diffusori **Householder/WHT/Hybrid**. Progetto C++17 per Max 8/9, con help patcher, preset e metadati package.

## Caratteristiche

* 16 linee di ritardo; predelay SR‑aware; early φ‑spaced con controllo stereo **M/S**.
* Detector **range–doppler** (RMS lento + 0.7–3 kHz) per modulare `width/damping/moddepth` con slew controllato da `@motion`.
* Modulatori per‑linea (`@modrate`, `@moddepth`).
* Controlli: `@regen @derez @filter @early @focus @predelay @mix @width @size @color @modrate @moddepth @motion @phiweight @mode_mix`.

> **Laser**: presente come attributi ma **disabilitato di default** (`KBEYOND_ENABLE_LASER=0`). Gli slider hanno effetto **no‑op** finché non abilitato a build‑time.

## Requisiti

* **Max SDK** 8+ (`MAX_SDK_ROOT` impostata).
* **CMake** ≥ 3.18; compilatore **C++17** (clang++ su macOS, MSVC su Windows).

## Build (CMake)

```bash
export MAX_SDK_ROOT=/percorso/al/max-sdk
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Se `ext.h` / `z_dsp.h` non si trovano, verificare il contenuto del Max SDK (vedi `scripts/verify_max_sdk.sh`).

## Test (facoltativi)

```bash
cmake -S . -B build -DKBEYOND_BUILD_TESTS=ON -DMAX_SDK_ROOT=$MAX_SDK_ROOT
cmake --build build -j
ctest --test-dir build -V
```

## Package Max

* Copia gli externals in `externals/` (non versionato).
* `help/kbeyond~.maxhelp` e `presets/` restano nelle relative cartelle.
* Aggiorna `package-info.json` se cambi versione/compatibilità.

## Note DSP

* `@predelay` in secondi (0–0.5 s), quantizzato a campioni; **min 1 sample** per evitare “range=0”.
* `@width` agisce su ER e tail M/S; proiezione L/R **ortonormale** (∑wL² = ∑wR² = 1).
* Mapping `@regen/@decay` con floor sui coefficienti per evitare collassi.

---
