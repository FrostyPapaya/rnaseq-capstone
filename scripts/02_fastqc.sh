#!/usr/bin/env bash
set -euo pipefail

# Run FastQC on all FASTQ files
mkdir -p results/qc/

fastqc data/raw/*.fastq.gz -o results/qc/ -t 4

echo "FastQC complete. Running MultiQC..."
multiqc results/qc/ -o results/qc/multiqc/

echo "Done! Open results/qc/multiqc/multiqc_report.html in a browser."