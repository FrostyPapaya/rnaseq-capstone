   #!/usr/bin/env bash
set -euo pipefail

   # Download raw FASTQ files from SRA for GSE96870
   # Dataset: Mouse brain RNA-seq (Cerebellum vs Cortex)

cd data/raw/

# 4 from spinal cord, 4 from cerebellum
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
echo "Downloading $srr..."
prefetch "$srr" --max-size 100G
fasterq-dump --split-files "$srr" -O .
echo "Compressing $srr..."
pigz -p 4 "${srr}"*.fastq
done

echo "Recording checksums..."
md5sum *.fastq.gz > md5sums.txt

echo "Done! All files downloaded and checksummed."

