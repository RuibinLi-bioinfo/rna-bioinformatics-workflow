# Workflow Notes

## 1. Data Download

Public sequencing data can be downloaded using tools such as Kingfisher or SRA Toolkit.

## 2. Quality Control

Raw FASTQ files are checked using FastQC. This step helps evaluate base quality, adapter contamination, sequence duplication, and other sequencing quality metrics.

## 3. Adapter Trimming

Trimmomatic is used to remove sequencing adapters. In this workflow, adapter trimming is performed before duplicate removal.

## 4. Duplicate Removal

PCR duplicates are removed before polyG/polyN trimming. This order is important because trimming can change read lengths and affect duplicate detection.

## 5. polyG / polyN Trimming

polyG and polyN tails are removed using cutadapt. Reads shorter than 36 nt are filtered out.

## 6. rRNA Filtering

Bowtie2 is used to align reads to a 45S pre-rRNA reference. Unaligned read pairs are retained for downstream analysis.

## 7. STAR Chimeric Alignment

STAR is used for chimeric alignment. This is important for RIC-seq analysis because reads may contain gapped or chimeric RNA interaction signals.

## 8. Separate R1 and R2 Alignment

For RIC-seq, R1 and R2 can also be aligned separately to avoid losing reads due to forced paired-end alignment constraints.
