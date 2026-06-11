#!/usr/bin/env bash
set -euo pipefail

# Force the script to use the newly generated English UTF-8 settings
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

echo "Cleaning old index..."
# Remove the old directory entirely. Salmon will create it automatically.
rm -rf salmon_index

echo "Building Salmon index..."

salmon index \
  -t data/reference/Mus_musculus.GRCm38.cdna.all.fa.gz  \
  -i salmon_index \
  -k 31 \
  -p 1 \
  2>&1 | tee salmon_index.log

echo "Done"
