#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

###############################################################################
# RIC-seq preprocessing + rRNA filtering + STAR chimeric alignment pipeline
#
# Recommended usage:
#   conda activate qc_align
#   bash RIC-seq_sub_pipeline_rb_modified.sh
#
# To override paths without editing this script:
#   BASE_DIR=/home/bioinfo/07_people/rbli1/RICseq \
#   SRR_ID=SRR8632820 \
#   ADAPTER_PATH=/path/to/TruSeq3-PE.fa \
#   bash RIC-seq_sub_pipeline_rb_modified.sh
###############################################################################

############################
# 0. User-editable settings
############################

SRR_ID="${SRR_ID:-SRR8632820}"

BASE_DIR="${BASE_DIR:-/home/bioinfo/07_people/rbli1/RICseq}"
FASTQ_DIR="${FASTQ_DIR:-${BASE_DIR}/FastqFiles}"
CLEAN_DIR="${CLEAN_DIR:-${BASE_DIR}/CleanData}"
FASTQC_DIR="${FASTQC_DIR:-${BASE_DIR}/Fastqc}"
LOG_DIR="${LOG_DIR:-${BASE_DIR}/Logs}"
STAR_OUT_DIR="${STAR_OUT_DIR:-${BASE_DIR}/STARoutput}"

R1_RAW="${R1_RAW:-${FASTQ_DIR}/${SRR_ID}_1.fastq.gz}"
R2_RAW="${R2_RAW:-${FASTQ_DIR}/${SRR_ID}_2.fastq.gz}"

ADAPTER_PATH="${ADAPTER_PATH:-/mnt/data/teacher/miniforge3/envs/qc_align/share/trimmomatic/adapters/TruSeq3-PE.fa}"

DEDUP_SCRIPT="${DEDUP_SCRIPT:-${BASE_DIR}/remove_duplicated_reads_ricPaper.pl}"

RRNA_DIR="${RRNA_DIR:-/home/bioinfo/07_people/rbli1/Reference_Genome/rRNA}"
RRNA_FASTA="${RRNA_FASTA:-${RRNA_DIR}/human_45S_pre_ribosomal_N5.fasta}"
RRNA_INDEX_DIR="${RRNA_INDEX_DIR:-${RRNA_DIR}/human_45S_pre_rRNA_index_bowtie2}"
RRNA_INDEX_PREFIX="${RRNA_INDEX_PREFIX:-${RRNA_INDEX_DIR}/human_45S_pre_rRNA_index}"

GENOME_DIR="${GENOME_DIR:-/home/bioinfo/07_people/rbli1/Reference_Genome/Human}"
GENOME_FASTA="${GENOME_FASTA:-${GENOME_DIR}/hg38.fa}"
GENOME_GTF="${GENOME_GTF:-${GENOME_DIR}/hg38.gtf}"
STAR_INDEX_DIR="${STAR_INDEX_DIR:-${GENOME_DIR}/hg38_index_star2711}"

THREADS_QC="${THREADS_QC:-4}"
THREADS_DOWNLOAD="${THREADS_DOWNLOAD:-8}"
THREADS_RRNA="${THREADS_RRNA:-8}"
THREADS_STAR="${THREADS_STAR:-16}"

MIN_LEN="${MIN_LEN:-36}"

RUN_DOWNLOAD="${RUN_DOWNLOAD:-false}"
RUN_RAW_FASTQC="${RUN_RAW_FASTQC:-true}"
RUN_PREPROCESS="${RUN_PREPROCESS:-true}"
RUN_BUILD_RRNA_INDEX="${RUN_BUILD_RRNA_INDEX:-false}"
RUN_FILTER_RRNA="${RUN_FILTER_RRNA:-true}"
RUN_BUILD_STAR_INDEX="${RUN_BUILD_STAR_INDEX:-false}"
RUN_STAR_SINGLE_END="${RUN_STAR_SINGLE_END:-true}"
RUN_STAR_PAIRED_END="${RUN_STAR_PAIRED_END:-false}"

############################
# 1. Helper functions
############################

log() {
    echo "[$(date +'%F %T')] $*" >&2
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"
}

require_file() {
    [[ -s "$1" ]] || die "File not found or empty: $1"
}

require_bowtie2_index() {
    [[ -s "${RRNA_INDEX_PREFIX}.1.bt2" || -s "${RRNA_INDEX_PREFIX}.1.bt2l" ]] || \
        die "Bowtie2 index not found: ${RRNA_INDEX_PREFIX}. Build it first or set RUN_BUILD_RRNA_INDEX=true."
}

require_star_index() {
    [[ -s "${STAR_INDEX_DIR}/Genome" ]] || \
        die "STAR index not found: ${STAR_INDEX_DIR}. Build it first or set RUN_BUILD_STAR_INDEX=true."
}

mkdir -p "$FASTQ_DIR" "$CLEAN_DIR" "$FASTQC_DIR" "$LOG_DIR" "$STAR_OUT_DIR"

############################
# 2. Optional download
############################

if [[ "$RUN_DOWNLOAD" == "true" ]]; then
    log "Downloading ${SRR_ID} with kingfisher"

    require_command kingfisher

    kingfisher get -r "$SRR_ID" -m ena-ascp \
        --ascp-ssh-key "${ASCP_SSH_KEY:-${HOME}/.aspera/connect/etc/asperaweb_id_dsa.openssh}" \
        --output-directory "$FASTQ_DIR" \
        --download-threads "$THREADS_DOWNLOAD"
fi

############################
# 3. Raw-data FastQC
############################

if [[ "$RUN_RAW_FASTQC" == "true" ]]; then
    log "Running FastQC on raw FASTQ files"

    require_command fastqc
    require_file "$R1_RAW"
    require_file "$R2_RAW"

    fastqc -o "$FASTQC_DIR" \
        --noextract \
        --threads "$THREADS_QC" \
        "$R1_RAW" "$R2_RAW"
fi

############################
# 4. Adapter trimming -> deduplication -> polyG/N trimming
############################

TRIM_R1="${CLEAN_DIR}/${SRR_ID}_trim_R1.fq.gz"
TRIM_R2="${CLEAN_DIR}/${SRR_ID}_trim_R2.fq.gz"
TRIM_U1="${CLEAN_DIR}/${SRR_ID}_trim_U1.fq.gz"
TRIM_U2="${CLEAN_DIR}/${SRR_ID}_trim_U2.fq.gz"

DEDUP_R1="${CLEAN_DIR}/${SRR_ID}_dedup_R1.fq.gz"
DEDUP_R2="${CLEAN_DIR}/${SRR_ID}_dedup_R2.fq.gz"

CLEAN_R1="${CLEAN_DIR}/${SRR_ID}_1_clean.fq.gz"
CLEAN_R2="${CLEAN_DIR}/${SRR_ID}_2_clean.fq.gz"

if [[ "$RUN_PREPROCESS" == "true" ]]; then
    log "Step 1/3: Trimmomatic adapter trimming"

    require_command trimmomatic
    require_file "$R1_RAW"
    require_file "$R2_RAW"
    require_file "$ADAPTER_PATH"

    trimmomatic PE -threads "$THREADS_QC" \
        "$R1_RAW" \
        "$R2_RAW" \
        "$TRIM_R1" \
        "$TRIM_U1" \
        "$TRIM_R2" \
        "$TRIM_U2" \
        "ILLUMINACLIP:${ADAPTER_PATH}:2:30:10"

    log "Step 2/3: Removing PCR duplicates before terminal-base trimming"

    require_command perl
    require_file "$DEDUP_SCRIPT"

    /usr/bin/time -v perl "$DEDUP_SCRIPT" \
        "$TRIM_R1" \
        "$TRIM_R2" \
        "$DEDUP_R1" \
        "$DEDUP_R2" \
        2> "${LOG_DIR}/${SRR_ID}_dedup_time.log"

    log "Step 3/3: Cutadapt polyG trimming, terminal N trimming, and length filtering"

    require_command cutadapt

    cutadapt \
        -a "G{10}$" \
        -A "G{10}$" \
        --trim-n \
        -e 0.1 \
        --minimum-length "$MIN_LEN" \
        -o "$CLEAN_R1" \
        -p "$CLEAN_R2" \
        "$DEDUP_R1" \
        "$DEDUP_R2" \
        > "${LOG_DIR}/${SRR_ID}_cutadapt.log"

    log "Running FastQC on final cleaned FASTQ files"

    fastqc -o "$FASTQC_DIR" \
        --noextract \
        --threads "$THREADS_QC" \
        "$CLEAN_R1" "$CLEAN_R2"
fi

############################
# 5. Optional Bowtie2 rRNA index construction
############################

if [[ "$RUN_BUILD_RRNA_INDEX" == "true" ]]; then
    log "Building Bowtie2 index for 45S pre-rRNA"

    require_command bowtie2-build
    require_file "$RRNA_FASTA"

    mkdir -p "$RRNA_INDEX_DIR"

    bowtie2-build "$RRNA_FASTA" "$RRNA_INDEX_PREFIX"
fi

############################
# 6. rRNA filtering
############################

NORRNA_R1="${CLEAN_DIR}/${SRR_ID}_NOrRNA_R1.fq.gz"
NORRNA_R2="${CLEAN_DIR}/${SRR_ID}_NOrRNA_R2.fq.gz"

RRNA_R1="${CLEAN_DIR}/${SRR_ID}_rRNA_R1.fq.gz"
RRNA_R2="${CLEAN_DIR}/${SRR_ID}_rRNA_R2.fq.gz"

if [[ "$RUN_FILTER_RRNA" == "true" ]]; then
    log "Filtering 45S pre-rRNA reads with Bowtie2"

    require_command bowtie2
    require_bowtie2_index
    require_file "$CLEAN_R1"
    require_file "$CLEAN_R2"

    bowtie2 -x "$RRNA_INDEX_PREFIX" \
        -1 "$CLEAN_R1" \
        -2 "$CLEAN_R2" \
        --un-conc-gz "${CLEAN_DIR}/${SRR_ID}_NOrRNA_R%.fq.gz" \
        --al-conc-gz "${CLEAN_DIR}/${SRR_ID}_rRNA_R%.fq.gz" \
        --no-unal \
        -p "$THREADS_RRNA" \
        -S /dev/null \
        2> "${LOG_DIR}/${SRR_ID}_bowtie2_rRNA.log"
fi

############################
# 7. Optional STAR genome index construction
############################

if [[ "$RUN_BUILD_STAR_INDEX" == "true" ]]; then
    log "Building STAR index"

    require_command STAR
    require_file "$GENOME_FASTA"
    require_file "$GENOME_GTF"

    mkdir -p "$STAR_INDEX_DIR"

    STAR --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_DIR" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$GENOME_GTF" \
        --runThreadN "$THREADS_STAR" \
        --genomeSAindexNbases 14
fi

############################
# 8. STAR chimeric alignment
############################

STAR_COMMON_ARGS=(
    --runMode alignReads
    --genomeDir "$STAR_INDEX_DIR"
    --readFilesCommand zcat
    --outSAMattributes All
    --outSAMtype BAM SortedByCoordinate
    --outFilterMultimapNmax 100
    --alignIntronMin 1
    --scoreGapNoncan -4
    --scoreGapATAC -4
    --chimSegmentMin 15
    --chimJunctionOverhangMin 15
    --alignSJoverhangMin 15
    --alignSJDBoverhangMin 10
    --alignSJstitchMismatchNmax 5 -1 5 5
    --runThreadN "$THREADS_STAR"
)

if [[ "$RUN_STAR_PAIRED_END" == "true" ]]; then
    log "Running optional paired-end STAR alignment"

    require_command STAR
    require_star_index
    require_file "$NORRNA_R1"
    require_file "$NORRNA_R2"

    STAR "${STAR_COMMON_ARGS[@]}" \
        --readFilesIn "$NORRNA_R1" "$NORRNA_R2" \
        --outFileNamePrefix "${STAR_OUT_DIR}/${SRR_ID}_PE_" \
        --chimOutType Junctions WithinBAM SoftClip
fi

if [[ "$RUN_STAR_SINGLE_END" == "true" ]]; then
    log "Running STAR alignment for R1 as single-end reads"

    require_command STAR
    require_star_index
    require_file "$NORRNA_R1"

    STAR "${STAR_COMMON_ARGS[@]}" \
        --readFilesIn "$NORRNA_R1" \
        --outFileNamePrefix "${STAR_OUT_DIR}/${SRR_ID}_R1_" \
        --chimOutType Junctions SeparateSAMold

    log "Running STAR alignment for R2 as single-end reads"

    require_file "$NORRNA_R2"

    STAR "${STAR_COMMON_ARGS[@]}" \
        --readFilesIn "$NORRNA_R2" \
        --outFileNamePrefix "${STAR_OUT_DIR}/${SRR_ID}_R2_" \
        --chimOutType Junctions SeparateSAMold
fi

log "Pipeline finished for ${SRR_ID}"
log "Clean FASTQ: ${CLEAN_R1}, ${CLEAN_R2}"
log "Non-rRNA FASTQ: ${NORRNA_R1}, ${NORRNA_R2}"
log "STAR output directory: ${STAR_OUT_DIR}"
