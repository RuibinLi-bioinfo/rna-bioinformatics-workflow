# RIC-seq and RNA Bioinformatics Workflow Demo

This repository documents a reproducible RNA bioinformatics workflow demo for RIC-seq preprocessing, RNA interaction analysis, rRNA filtering, and STAR-based chimeric alignment.

## About Me

I am Ruibin Li, an undergraduate student majoring in Biotechnology at Shandong University, with a second major in Financial Mathematics and Financial Engineering.

My research interests include RNA bioinformatics, computational biology, biostatistics, RNA interaction analysis, circRNA translation, and regulatory genomics.

I am currently receiving bioinformatics research training related to super enhancers, circRNA translation, RIC-seq, and RNA interaction analysis.

## Workflow Overview

This workflow includes:

1. Public sequencing data download
2. Raw FASTQ quality control
3. Adapter trimming
4. PCR duplicate removal
5. polyG / polyN trimming and length filtering
6. rRNA filtering using Bowtie2
7. STAR genome index construction
8. STAR-based chimeric alignment
9. Separate R1 and R2 alignment for RIC-seq analysis
10. Workflow documentation for reproducible research

## Technical Skills Demonstrated

- Linux command line
- Shell scripting
- Miniconda / conda environment management
- Public sequencing data download
- FASTQ quality control
- Trimmomatic-based adapter trimming
- PCR duplicate removal
- cutadapt / fastp-based read filtering
- Bowtie2-based rRNA filtering
- STAR genome index construction
- STAR chimeric alignment
- RIC-seq preprocessing logic
- Reproducible workflow documentation

## Notes

This repository is a workflow demonstration based on public sequencing data analysis logic. Large raw sequencing files such as FASTQ, BAM, SAM, SRA files, reference genomes, and genome index files are not included.

## Contact

Ruibin Li  
Shandong University  
Email: 3652812409@qq.com

## Repository Structure

```text
rna-bioinformatics-workflow/
│
├── README.md
├── environment.yml
├── sample_metadata.csv
│
├── scripts/
│   └── RIC-seq_sub_pipeline_rb.sh
│
├── docs/
│   ├── workflow_notes.md
│   └── research_summary.md
│
└── results/
    └── example_qc_summary.csv
