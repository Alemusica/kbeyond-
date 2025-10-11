#!/usr/bin/env bash
# Richiede GitHub CLI: gh auth login

set -euo pipefail

REPO="${1:-Alemusica/kbeyond-}"

# Milestones
gh api repos/$REPO/milestones -X POST -f title='A Hardening' -f state='open' || true
gh api repos/$REPO/milestones -X POST -f title='B Power' -f state='open' || true
gh api repos/$REPO/milestones -X POST -f title='C Quantum/Laser' -f state='open' || true
gh api repos/$REPO/milestones -X POST -f title='D PrimeFear' -f state='open' || true
gh api repos/$REPO/milestones -X POST -f title='E GPU' -f state='open' || true

issue() {
  local title="$1"; local body="$2"; local labels="$3"; local milestone="$4"
  gh issue create --repo "$REPO" --title "$title" --body "$body" --label $labels --milestone "$milestone"
}

issue "A1: Early energy-balanced"         "* Tap ER specchiati, pan equal-power, width in M/S ortonormale
* Test impulso mono: |L−R| < 0.1 dB
* Update docs/preset"         "task,dsp,phase:A" "A Hardening"

issue "A2: RT60 vero + damping 3-bande"         "* Param @decay (s) con g_i = 10^(−3·d_i/@decay)
* 1-pole LF/MF/HF per linea
* QA RT60 per banda ≤ ±5%"         "task,dsp,phase:A" "A Hardening"

issue "A3: Diffusori + micro-jitter"         "* 2–4 allpass con jitter lento
* Micro-jitter ±1–2 smp su pochi tap"         "task,dsp,phase:A" "A Hardening"

issue "B1: Mixing WHT/Hybrid"         "* Aggiungere modalità WHT/Hybrid oltre Householder(φ)"         "task,dsp,phase:B" "B Power"

issue "B2: Detector range–doppler"         "* RMS lento + banda 0.7–3 kHz
* Modula width/hf_damp/modDepth con slew"         "task,dsp,phase:B" "B Power"

issue "C1: PHONON (CPU-HQ)"         "* Givens a blocchi tra micro-modi
* Parametri @phonon_*"         "task,dsp,phase:C" "C Quantum/Laser"

issue "C2: Early radar + Q-switch"         "* Cluster chirp/m-seq
* Q-switch 200–600 ms"         "task,dsp,phase:C" "C Quantum/Laser"

issue "D1: PrimeModes + routing M/S"         "* prime/aureo/plastica/prime_aureo
* @mode_er/@mode_late/@mode_mid/@mode_side/@seed"         "task,dsp,phase:D" "D PrimeFear"

issue "D2: Guard-rails PrimeFear"         "* Co-primalità linee; min spacing 2 smp; auto-bilanciamento stereo"         "task,dsp,phase:D" "D PrimeFear"

issue "E1: GPU LateCon + fallback"         "* FFT partizionata MPSGraph; prima partizione CPU/vDSP; fallback"         "task,gpu,phase:E" "E GPU"
