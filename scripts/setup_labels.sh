#!/usr/bin/env bash
# Crea etichette utili. Richiede gh CLI.

REPO="${1:-Alemusica/kbeyond-}"

labels=(
  "task:#0366d6"
  "bug:#d73a4a"
  "docs:#0e8a16"
  "dsp:#5319e7"
  "gpu:#1d76db"
  "phase:A:#B60205"
  "phase:B:#D93F0B"
  "phase:C:#FBCA04"
  "phase:D:#0E8A16"
  "phase:E:#0052CC"
)

for l in "${labels[@]}"; do
  name="${l%%:*}"
  color="${l##*:}"
  gh label create "$name" --repo "$REPO" --color "$color" --force
done
