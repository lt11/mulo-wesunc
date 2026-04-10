# Somatic Variant Detection with "Ensable" Calling

## Overview

This pipeline performs ensemble variant calling using GATK Mutect2 and VarDict, two somatic variant callers that support both matched (tumour–normal) and unmatched (tumour-only) modes.

This workflow is designed for tumour-only analyses.

The corresponding matched-normal pipeline "mulo-wesml" is available here: https://github.com/lt11/mulo-wesml.

## Installation

First, users should install SnpEff and SnpSift, and set their paths to in "mulo-wesml.sh" (using the variables "snpeff_path" and "snpsift_path").

Then, install the other dependences:

- fastqc
- bwa
- samtools
- bcftools
- GATK
- VarDict and its helper scripts (var2vcf_paired.pl, var2vcf_valid.pl, teststrandbias.R, and testsomatic.R)
- vt
- bgzip
- tabix

---

## Input Data

### FASTQ Files

The input FASTQ file must be placed in the "exp" folder. Their IDs must be reported in the "popu_samp" variable of "mulo-wesml.sh". E.g. given a tumour-only sample "smp1t", the user should set:
```
popu_samp="smp1t"
```
The "exp" folder should then contain the following FASTQ files:
```
smp1t-R1.fq.gz
smp1t-R2.fq.gz
```

### Reference Genome and SnpEff database

The reference genome, in FASTA format, should be placed in the "rep" folder. The user should verify that the reference genome matches the ID of the reference genome reported in "anno-snpeff-unc.sh". The "ref_gen_ver" variable in the script should be edited accordingly.

The path to the databases used by SnpEff must be reported in "anno-snpeff-unc.sh". These include:
- dbNSFP (variable: "dbnsfp_path")
- cosmic coding variants (variable: "cosmcod_path")
- cosmic non-coding variants (variable: "cosmnonc_path")
- dbSNP (variable: "dbsnp_path")

### Target Regions

The regions of the capture kit used to perform the whole-exome sequencing experiments must be reported, in BED format, in the "cpk" folder.

---

## Usage

Once the "mulo-wesunc.sh" script has been edited as reported above, running the pipeline is as simple as:
```
bash mulo-wesml.sh
```

---

## Outputs

The variants called are reported in the VCF files in the "var-calls" folder.

---

## Citation

If you use this code, please cite:

Tattini, L., Yan, Y., Chaturvedi, N., & Appuswamy, R. (2025). Accurate Variant Classification in Tumour-Only Genomic Data Using Interpretable Tabular Models. bioRxiv, 2025-12.
doi: https://doi.org/10.64898/2025.12.09.693348
