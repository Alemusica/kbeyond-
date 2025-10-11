# kbeyond~ Max/MSP external

`kbeyond~` è un riverbero FDN 16×16 che utilizza una trasformata di Householder pesata con φ per creare una coda spaziosa e modulata. Il progetto include tutto il necessario per compilare l'external con il Max SDK su macOS e Windows, insieme a help patcher, preset e metadati di package.

## Caratteristiche principali
- 16 linee di ritardo con matrice Householder(φ) calcolata in tempo reale.
- Predelay sensibile al sample-rate, riflessioni precoci φ-spaced con controllo stereo M/S.
- Modulatori indipendenti per ciascuna linea FDN con parametri `@modrate` e `@moddepth`.
- Controlli completi: `@regen @derez @filter @early @predelay @mix @width @size @color @modrate @moddepth @phiweight`.
- Implementazione C++17, 64-bit only, compatibile con Max 8/9.

## Struttura del repository
```
source/                  # Sorgenti C++ dell'external
help/kbeyond~.maxhelp     # Patcher di esempio con controlli mix/width
presets/kbeyond_presets.json
scripts/                  # Script di build automatica per macOS e Windows
Makefile.macos            # Build rapida clang++ -> bundle .mxo
CMakeLists.txt            # Configurazione CMake cross-platform
package-info.json         # Metadati del Max package
README.md / CHANGELOG.md / LICENSE / .gitignore
```

## Requisiti
- Max SDK 8 o successivo. Impostare la variabile d'ambiente `MAX_SDK_ROOT` sul percorso root del Max SDK (ad esempio `/Applications/Max.app/Contents/Resources/C74/packages/max-sdk-8.5.5`).
- CMake ≥ 3.18, un compilatore C++17 (clang++ su macOS, MSVC su Windows).

## Build con CMake
```bash
export MAX_SDK_ROOT=/percorso/al/max-sdk   # oppure passare -DMAX_SDK_ROOT=/percorso/al/max-sdk a cmake
./scripts/verify_max_sdk.sh               # controllo opzionale ma consigliato
cmake -S . -B build/Release -DCMAKE_BUILD_TYPE=Release
cmake --build build/Release --config Release
```
Se `cmake` termina immediatamente con un errore che menziona `ext.h` o `z_dsp.h`, significa che il percorso fornito in `MAX_SDK_ROOT` (o passato tramite `-DMAX_SDK_ROOT`) non contiene un Max SDK valido (ad esempio la cartella `source/c74support` manca oppure è una versione incompleta). Assicurarsi di avere scaricato e scompattato il Max SDK ufficiale e riprovare. Lo script `scripts/verify_max_sdk.sh` permette di controllare rapidamente che i file richiesti (`ext.h`, `z_dsp.h`) siano presenti prima di configurare la build.

Sono supportati sia l'archivio `.zip` distribuito da Cycling '74 (layout `max-sdk-base/c74support`) sia il checkout GitHub del repository del Max SDK (layout `source/c74support`). Il rilevamento automatico userà il percorso corretto per gli include e le framework path.
Gli artefatti prodotti (`kbeyond~.mxo` su macOS, `kbeyond~.mxe64` su Windows) verranno creati in `build/Release`. Copiare il file risultante nella cartella `externals/` del package Max.

## Script automatizzati
- macOS: `scripts/build_macos.sh`
- Windows (PowerShell): `scripts/build_windows.ps1`
Entrambi richiedono `MAX_SDK_ROOT` impostata e generano gli artefatti nelle sottocartelle `build/macos-cmake` o `build/windows-cmake`.

## Makefile macOS
Per una build rapida del bundle `.mxo` direttamente con clang++:
```bash
export MAX_SDK_ROOT=/percorso/al/max-sdk
make -f Makefile.macos
```
Il bundle è scritto in `build/macos/kbeyond~.mxo/Contents/MacOS/kbeyond~`.

## Package Max
- Copiare gli externals compilati in `externals/` (non versionato in repo).
- Lasciare i file di help e preset nelle rispettive cartelle.
- Aggiornare `package-info.json` se si modifica la versione o le compatibilità.

## Documentazione
- `help/kbeyond~.maxhelp` mostra una catena minima `noise~ -> kbeyond~ -> ezdac~` con controllo degli attributi `@mix` e `@width`.
- `presets/kbeyond_presets.json` offre tre preset (Hall Larga, Plate Densa, Room Intima) per un richiamo rapido.

## Note di implementazione DSP
- Il predelay (`@predelay`) è espresso in secondi (0–0.5 s) e ricampionato al cambio di sample-rate.
- La larghezza stereo (`@width`) agisce sia sulle riflessioni precoci (M/S) sia sul mix di uscita FDN.
- Le linee FDN utilizzano `householder_phi16.h` per generare il vettore normalizzato e applicare la trasformata senza creare la matrice completa.
- Non sono presenti riferimenti a vecchi componenti Audio Unit; il target espone esclusivamente `ext_main`, `dsp64`, `perform64` per Max.

## Licenza
MIT — vedere `LICENSE`.
