#!/bin/bash

# Count Ribo-seq RPF reads in candidate PC5UTR regions
# Large BAM files are not included in this repository.

set -euo pipefail

BED="projects/02_HeLa_RiboSeq_PC5UTR/data/candidate_PC5UTR_regions.bed"

BAM1="/path/to/SRR25602014.sorted.bam"
BAM2="/path/to/SRR25602013.sorted.bam"

OUT="projects/02_HeLa_RiboSeq_PC5UTR/results/PC5UTR_RPF_raw_counts.tsv"

bedtools multicov \
  -bams "$BAM1" "$BAM2" \
  -bed "$BED" \
  > "$OUT"
