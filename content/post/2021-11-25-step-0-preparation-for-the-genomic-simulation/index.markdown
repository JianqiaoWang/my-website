---
title: 'Step 0: Preparation for the genomic simulation'
author: Jianqiao Wang
date: '2021-11-25'
slug: step-0-preparation-for-the-genomic-simulation
categories: []
tags: []
subtitle: ''
summary: ''
authors: []
lastmod: '2021-11-25T16:17:55-05:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---


To do simulation with genetic data, we need to make some preparion. First, we detemine the snplist and individuals in the genotype data:

snplist.txt

Prepare data for package bigsnpr and bigstatsr


```r
library(bigsnpr)
library(bigstatsr)
PrepareX = FALSE

#-------------------------------------------------------
if(PrepareX){
rds = bigsnpr::snp_readBed2(bedfile = "../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01.bed",
                     backingfile = "cric_geno",ind.col= ind.column)
geno <- snp_attach(rds)
G <- geno$genotypes
CHR <- geno$map$chromosome
infos <- snp_fastImputeSimple(G, method = "mean2")
# To make this permanent, you need to save (modify) the file on disk
geno$genotypes <- infos
geno <- snp_save(geno)
G <- geno$genotypes
ind.keep <- snp_clumping(G, infos.chr = geno$map$chromosome,
                        infos.pos = geno$map$physical.pos,
                        thr.r2 = 0.01)
saveRDS(ind.keep, file ="ind_keep_clump.rds")
data.table::fwrite( as.data.frame(geno$map$marker.ID), file = "snplist.txt")
}
```

Prepare for the GCTA method 


```r
ref.bed = "../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01"
MakeGRM = function(ref.bed, keep.snp, keep.ind, output.dir = "./temp/"){

  grm = paste0(output.dir, "geno")

  fwrite(data.frame(keep.snp), file = paste0(output.dir,"snplist.txt"), quote = F, sep = "\t", col.names=T, na = "NA")

  system(paste0("plink --bfile ", ref.bed, " --extract ", paste0(output.dir, "snplist.txt"), " --make-grm-bin --out ", grm ) )
}
# alternative, slower
#system(paste0("./gcta64 --bfile ", ref.bed, " --extract ", paste0(output.dir, "snplist.txt"), " --make-grm-bin --out ", grm #) )
```

Prepare for the LD Score regression : Calculate the ld score


```r
System("ldsc.py --bfile ../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01 --extract ./temp/snplist.txt --ld-wind-kb 1000 --out ./temp/geno")
```
