---
title: Working with the reference genotype data
author: Jianqiao Wang
date: '2021-11-27'
slug: working-with-the-reference-genotype-data
categories: []
tags: []
subtitle: ''
summary: ''
authors: []
lastmod: '2021-11-27T17:57:24-05:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---

# Introduction

The reference genotype database has been widely used for GWA studieds. It provides more information on the genetic variants when the individual-level genotype data is not available or the genotypes data has low sequencing depth. The reference also has been used for imputing the missing genotypes to boost the power of GWAS.


Common reference panels include:
- HapMap3 reference panels: Early reference panels 
- 1000 Genomes Phase3  (1KGP, widely used), 2506 individuals from 26 populations, 49,143,605 Sites. 504 individuals for European ancestry.
- The Haplotype Reference Consortium (HRC), the first release consists of 64,976 haplotypes at 39,235,157 SNPs. Mainly european ancestry but also includes 1KGP
- Topmed reference panel. Version r2 of the panel includes 97,256 reference samples and 308,107,085 genetic variants distributed across the 22 autosomes and the X chromosome.


Not all reference panels are public. Large reference needs some applications . Here we focus on working with 1KGP reference panel.

# Work with 1KGP

## access to the 1KGP project data 

- Raw genotype varaint call files (vcf) could be downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ at each chromosome. However, we remind the loci of variants were determined for version GRCh37. 
  - Variants on GRCh38 should check the information about the data collection at : https://www.internationalgenome.org/data-portal/data-collection/30x-grch38. It contains 2504 unrelated individuals and additional 698 related samples with high coverage. 
  
To convert them to the plink file, we need to run following commands:

```shell
for j in {1..22}; do
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$j.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
./plink --make-bed --out chr$j --maf 0.01 --vcf ALL.chr$j.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done
```
  
- The plink file format of 1KGP could be found in https://www.coggenomics.org/plink/2.0/resources#1kg_phase3 in all_phase3.pgen, .psam, .pvar format


```bash
plink2 --pgen all_phase3.pgen --pgen-info # check the information 
```







```bas
--pgen-info on all_phase3.pgen:
  Variants: 84805772
  Samples: 2504
  REF alleles are all known
  Maximum allele count for a single variant: >2, not explicitly stored
  Explicitly phased hardcalls present
  No dosages present
```



```batch
plink2 --pfile $refdir/all_phase3 vzs\
--max-alleles 2
--make-bed \
--out $refdir/all_phase3
mv $refdir/all_phase3.log $refdir/log
```


- A helpful tool for analyzing the GWAS data: plinkQC (Discussed later)
  






