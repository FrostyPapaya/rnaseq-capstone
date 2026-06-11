#!/usr/bin/env bash
set -euo pipefail

#force the scripts to use the english UTF-8 settings
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

mkdir -p data/processed/salmon/
mkdir -p logs/
SRR_IDS=(
  "SRR5364318"
  "SRR5364317"
  "SRR5364316"
  "SRR5364315"
  "SRR5364337"
  "SRR5364336"
  "SRR5364335"
  "SRR5364334"
)
for srr in "${SRR_IDS[@]}"; do
  echo "Quantifying $srr..."
  salmon quant \
    -i salmon_index/ \
    -l A \
    -1 data/raw/${srr}_1.fastq.gz \
    -2 data/raw/${srr}_2.fastq.gz \
    -p 4 \
    --validateMappings \
    --gcBias \
    -o data/processed/salmon/${srr} \
    2> logs/salmon_${srr}.log
  echo "$srr done!"
done
echo "All samples quantified!"
