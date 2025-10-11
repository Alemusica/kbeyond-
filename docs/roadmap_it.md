# Roadmap di sviluppo kbeyond~

Questo documento riassume le fasi evolutive pianificate per l'external `kbeyond~`, dalla stabilizzazione iniziale fino all'ottimizzazione GPU.

## Fase A — Hardening

1. **A1. Early bilanciata**
   - Implementare tap delle riflessioni precoci specchiati tra i canali.
   - Utilizzare un bilanciamento in potenza (`pan` equal-power) con controllo della larghezza basato su M/S ortonormale.
   - Validare con un test d'impulso mono assicurandosi che lo sbilanciamento |Δ| rimanga inferiore a 0,1 dB.
2. **A2. RT60 reale + Materiali a 3 bande**
   - Modellare il decadimento con un parametro `@decay` espresso in secondi.
   - Applicare per ogni linea il guadagno `g_i = 10^{(-3 \cdot d_i / @decay)}`.
   - Integrare filtri mono polo separati per basse, medie e alte frequenze.
3. **A3. Vitalità a bassa CPU**
   - Inserire 2–4 diffusori all-pass con jitter lento per aumentare la densità.
   - Aggiungere micro-jitter di ±1–2 campioni su un sottoinsieme di tap per evitare staticità.

## Fase B — Potenza (Spazio/Dinamica)

1. **B1. Mixing WHT/Hybrid**
   - Estendere la matrice di diffusione oltre la trasformata di Householder pesata con φ adottando una combinazione Walsh-Hadamard/ibrida.
   - Massimizzare la densità temporale tramite somme e sottrazioni bilanciate.
2. **B2. Detector range-doppler**
   - Implementare un rivelatore range-doppler sui segnali interni.
   - Mappare l'uscita del detector sui parametri `width`, `hf_damp` e `modDepth` con slew controllato.

## Fase C — Quantum / Laser

1. **C1. Dither unitario + Quantum walk**
   - Applicare un dither controllato unitario.
   - Introdurre un processo di quantum walk regolato dai parametri `@coherence` e `@uwalk_rate`.
2. **C2. Laser**
   - Creare cluster di riflessioni precoci "mode-locked" utilizzando chirp o sequenze m.
   - Aggiungere un Q-switch di diffusione nei primi 200–600 ms della coda per ottenere transizioni dinamiche.

## Fase D — PrimeFear

1. **D1. PrimeModes (nuovo layer)**
   - Definire modalità basate su numeri primi, sezione aurea e proporzioni plastiche sia per la ER sia per la FDN.
   - Offrire routing differenziato M/S tra le modalità.
2. **D2. Guard-rails**
   - Assicurare la coprimalità delle lunghezze delle linee e una distanza minima tra tap di almeno 2 campioni.
   - Implementare un sistema di auto-bilanciamento stereo.

## Fase E — HYQ GPU (Metal)

1. **E0. Manager Metal/MPSGraph**
   - Gestire device e queue Metal con buffer `Shared`.
   - Utilizzare MPSGraph per FFT/iFFT batched con scheduling efficiente.
   - Integrare strumenti Apple Developer per profiling e ottimizzazione.
2. **E1. LateCon partizionata su GPU**
   - Spostare la convoluzione tardiva partizionata sulla GPU con fallback vDSP.
   - Applicare interventi "quantum/laser" su fase e ampiezza per ogni bin di frequenza.

## Criteri di accettazione

- Ogni fase deve essere completata senza regressioni nei test DSP esistenti.
- Aggiornare changelog e documentazione utente dopo ogni milestone.
- Validare le nuove funzionalità con patcher Max dedicati e misurazioni oggettive.

