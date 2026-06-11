#!/usr/bin/env bash
set -euo pipefail
# Ensembl release 102, GRCm38 (mouse mm10)
# Change to GRCh38 if human
cd data/reference/
echo "Downloading mouse cDNA transcriptome..."
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
echo "Downloading genome for decoys..."
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
echo "Building decoy list..."
grep "^>" <(gunzip -c Mus_musculus.GRCm38.dna.toplevel.fa.gz) | \
  cut -d " " -f 1 | \
  sed 's/>//' > decoys.txt
echo "Concatenating transcriptome + genome for selective alignment..."
cat Mus_musculus.GRCm38.cdna.all.fa.gz \
    Mus_musculus.GRCm38.dna.toplevel.fa.gz > gentrome.fa.gz
echo "Reference files ready!"
