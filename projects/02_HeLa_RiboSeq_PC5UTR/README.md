# HeLa Ribo-seq Analysis of Candidate PC5UTR Regions

This project summarizes an ongoing analysis of ribosome-protected fragment signals in candidate PC5UTR / 5′UTR regions using public HeLa Ribo-seq datasets.

## Background

5′UTR regions may contain upstream open reading frames, regulatory elements, or translation-related signals. In this project, I analyzed candidate PC5UTR regions to evaluate whether they show ribosome-protected fragment coverage in public HeLa Ribo-seq data.

## Research Question

Do candidate PC5UTR / 5′UTR regions show detectable RPF coverage in uninduced HeLa Ribo-seq samples?

## Data

The analysis uses two public uninduced HeLa Ribo-seq samples:

```text
SRR25602014
SRR25602013
```

Large FASTQ, BAM, and reference genome files are not included in this repository.

## Workflow Overview

1. Organize candidate PC5UTR regions in BED format
2. Use processed Ribo-seq alignment files or public Ribo-seq data
3. Count RPF reads overlapping each candidate PC5UTR region
4. Calculate region length
5. Normalize RPF counts using CPM and RPKM
6. Compare RPF signal across two biological replicates
7. Classify regions based on signal strength and replicate consistency

## Preliminary Findings

The current analysis suggests that several candidate PC5UTR regions show reproducible RPF coverage in two uninduced HeLa Ribo-seq replicates.

Regions with stronger RPF signals include:

* SRSF2_TypeB_PC5UTR
* VIM_TypeB_PC5UTR
* TOP2B_TypeA_PC5UTR
* SMC4_TypeB_PC5UTR
* YY1AP1_TypeB_PC5UTR
* ASH1L_TypeA_PC5UTR

Regions with weaker or unclear RPF signals include:

* INPP5A_TypeA_PC5UTR
* JRK_TypeB_PC5UTR
* USP24_TypeA_PC5UTR

Regions without clear RPF signal in the current analysis include:

* ICE1_TypeB_PC5UTR
* ZNF473_TypeA_PC5UTR
* PPP4R1_TypeA_PC5UTR

## Interpretation

RPF coverage suggests ribosome occupancy in candidate PC5UTR regions. However, this result alone does not prove active translation. Further analysis is needed, including IGV visualization, P-site offset analysis, 3-nt periodicity, ORF annotation, and comparison with matched RNA-seq expression.

## Next Steps

* Visualize high-signal regions in IGV
* Examine P-site distribution and 3-nt periodicity
* Compare candidate PC5UTR regions with background 5′UTR regions from GENCODE annotation
* Integrate RNA-seq expression data to estimate translation efficiency
* Expand from 12 candidate regions to a larger matched control analysis
