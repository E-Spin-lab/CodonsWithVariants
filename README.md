# CodonsWithVariants
This repository includes CodonsWithVariants, a Python script designed to analyze variants exported from the CLC Genomics Workbench. The script assesses whether single nucleotide polymorphisms (SNPs) are situated within the mRNA transcript and the coding sequence (CDS) of a gene. Should a SNP be located within the CDS, the corresponding codon number will be identified. For this analysis, it is essential to include a the reference genome's gff3 file.

## Background
Script was created at the request of a coworker for the analysis of low-frequency mutations that emerged during a viral challenge in animals. Variant tables were generated using CLC Genomics Workbench by aligning processed Illumina paired-end read files with the reference genome.

## Required Libraries
pandas

## Usage
This script was designed for a particular experiment and not for deployment, and therefore, command-line arguments were not incorporated. Future applications involving different reference genomes and experiments would require modifications to lines 4-6.
