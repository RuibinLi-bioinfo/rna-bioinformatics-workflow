#!/bin/bash

# RIC-seq preprocessing and STAR chimeric alignment workflow
# Author: Ruibin Li
# Purpose: Demonstration of a reproducible RNA bioinformatics workflow

set -euo pipefail

############################
# 0. Basic configuration
############################

BASEpath="/home/bioinfo/07_people/rbli1/RICseq"
SRRid="SRR8632820"

rRNAgenomPATH="/home/bioinfo/07_people/rbli1/Reference_Genome/rRNA"
genomPATH="/home/bioinfo/07_people/rbli1/Reference_Genome/Human"

mkdir -p "$BASEpath/Fastqc"
mkdir -p "$BASEpath/CleanData"
mkdir -p "$BASEpath/Logs"
mkdir -p "$BASEpath/STARoutput"

############################
# 1. Download public sequencing data
############################

# Activate kingfisher environment before running this step:
# mamba activate kingfisher

kingfisher get -r "$SRRid" -m ena-ascp \
  --ascp-ssh-key ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  --output-directory /home/bioinfo/07_people/rbli1 \
  --download-threads 8

############################
# 2. Initial FastQC
############################

# Activate QC environment before running this step:
# conda activate qc_align

fastqc -o "$BASEpath/Fastqc" \
  --noextract \
  --threads 4 \
  "$BASEpath/FastqFiles/${SRRid}_1.fastq.gz" \
  "$BASEpath/FastqFiles/${SRRid}_2.fastq.gz"

############################
# 3. Adapter trimming using Trimmomatic
############################

ADAPTER_PATH="/mnt/data/teacher/miniforge3/envs/qc_align/share/trimmomatic/adapters/TruSeq3-PE.fa"

trimmomatic PE -threads 4 \
  "$BASEpath/FastqFiles/${SRRid}_1.fastq.gz" \
  "$BASEpath/FastqFiles/${SRRid}_2.fastq.gz" \
  "$BASEpath/CleanData/${SRRid}_trim_R1.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_trim_U1.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_trim_R2.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_trim_U2.fq.gz" \
  ILLUMINACLIP:${ADAPTER_PATH}:2:30:10

############################
# 4. PCR duplicate removal
############################

/usr/bin/time -v perl "$BASEpath/remove_duplicated_reads_ricPaper.pl" \
  "$BASEpath/CleanData/${SRRid}_trim_R1.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_trim_R2.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_dedup_R1.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_dedup_R2.fq.gz" \
  2> "$BASEpath/Logs/${SRRid}_dedup_time.log"

############################
# 5. polyG / polyN trimming and length filtering
############################

# Note:
# The recommended order is:
# adapter trimming -> duplicate removal -> polyG/polyN trimming.
# polyN trimming can change read length, which may affect duplicate detection.

cutadapt \
  -a "G{10}" -A "G{10}" \
  -a "N{10}" -A "N{10}" \
  -e 0.1 \
  --minimum-length 36 \
  -o "$BASEpath/CleanData/${SRRid}_1_clean.fq.gz" \
  -p "$BASEpath/CleanData/${SRRid}_2_clean.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_dedup_R1.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_dedup_R2.fq.gz"

############################
# 6. FastQC after cleaning
############################

fastqc -o "$BASEpath/Fastqc" \
  --noextract \
  --threads 4 \
  "$BASEpath/CleanData/${SRRid}_1_clean.fq.gz" \
  "$BASEpath/CleanData/${SRRid}_2_clean.fq.gz"

############################
# 7. Build Bowtie2 index for 45S pre-rRNA
############################

mkdir -p "$rRNAgenomPATH/human_45S_pre_rRNA_index_bowtie2"

cd "$rRNAgenomPATH/human_45S_pre_rRNA_index_bowtie2"

bowtie2-build \
  "$rRNAgenomPATH/human_45S_pre_ribosomal_N5.fasta" \
  human_45S_pre_rRNA_index

############################
# 8. Remove 45S pre-rRNA reads using Bowtie2
############################

bowtie2 \
  -x "$rRNAgenomPATH/human_45S_pre_rRNA_index_bowtie2/human_45S_pre_rRNA_index" \
  -1 "$BASEpath/CleanData/${SRRid}_1_clean.fq.gz" \
  -2 "$BASEpath/CleanData/${SRRid}_2_clean.fq.gz" \
  --un-conc-gz "$BASEpath/CleanData/${SRRid}_NOrRNA" \
  --al-conc-gz "$BASEpath/CleanData/${SRRid}_rRNA" \
  --no-unal \
  -p 8 \
  -S /dev/null \
  2> "$BASEpath/Logs/${SRRid}_bowtie2.log"

############################
# 9. Build STAR genome index
############################

mkdir -p "$genomPATH/hg38_index_star2711"

STAR --runMode genomeGenerate \
  --genomeDir "$genomPATH/hg38_index_star2711" \
  --genomeFastaFiles "$genomPATH/hg38.fa" \
  --sjdbGTFfile "$genomPATH/hg38.gtf" \
  --runThreadN 8 \
  --genomeSAindexNbases 14

############################
# 10. STAR paired-end chimeric alignment
############################

STAR --runMode alignReads \
  --genomeDir "$genomPATH/hg38_index_star2711" \
  --readFilesIn "$BASEpath/CleanData/${SRRid}_NOrRNA.1" "$BASEpath/CleanData/${SRRid}_NOrRNA.2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$BASEpath/STARoutput/${SRRid}_" \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 100 \
  --alignIntronMin 1 \
  --scoreGapNoncan -4 \
  --scoreGapATAC -4 \
  --chimSegmentMin 15 \
  --chimJunctionOverhangMin 15 \
  --alignSJoverhangMin 15 \
  --alignSJDBoverhangMin 10 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimOutType Junctions WithinBAM SoftClip \
  --runThreadN 16

############################
# 11. STAR single-end chimeric alignment for R1
############################

STAR --runMode alignReads \
  --genomeDir "$genomPATH/hg38_index_star2711" \
  --readFilesIn "$BASEpath/CleanData/${SRRid}_NOrRNA.1" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$BASEpath/STARoutput/${SRRid}_R1_" \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 100 \
  --alignIntronMin 1 \
  --scoreGapNoncan -4 \
  --scoreGapATAC -4 \
  --chimSegmentMin 15 \
  --chimJunctionOverhangMin 15 \
  --alignSJoverhangMin 15 \
  --alignSJDBoverhangMin 10 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimOutType Junctions SeparateSAMold \
  --runThreadN 16

############################
# 12. STAR single-end chimeric alignment for R2
############################

STAR --runMode alignReads \
  --genomeDir "$genomPATH/hg38_index_star2711" \
  --readFilesIn "$BASEpath/CleanData/${SRRid}_NOrRNA.2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$BASEpath/STARoutput/${SRRid}_R2_" \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 100 \
  --alignIntronMin 1 \
  --scoreGapNoncan -4 \
  --scoreGapATAC -4 \
  --chimSegmentMin 15 \
  --chimJunctionOverhangMin 15 \
  --alignSJoverhangMin 15 \
  --alignSJDBoverhangMin 10 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimOutType Junctions SeparateSAMold \
  --runThreadN 16
