---
title: Simulation and analysis with R, GCTA, and LDSC
author: Jianqiao Wang
date: '2021-11-23'
slug: simulation-and-analysis-with-r-and-ldsc
categories: []
tags: []
subtitle: ''
summary: 'Organize simulation preparation and provide a reusable R workflow for GCTA and LDSC.'
authors: []
lastmod: '2021-11-23T11:50:44-05:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---

My work often compares results from GCTA-GREML and LDSC. This note summarizes the end-to-end workflow—from data preparation to calling GCTA/LDSC from R—in a reusable, extensible form.

## 0. Preparation

Prepare the SNP list, sample list, and LD scores for simulation.

### 0.1 Prepare data for bigsnpr/bigstatsr

```r
library(bigsnpr)
library(bigstatsr)
PrepareX = FALSE

if(PrepareX){
  rds = bigsnpr::snp_readBed2(bedfile = "../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01.bed",
                       backingfile = "cric_geno",ind.col= ind.column)
  geno <- snp_attach(rds)
  G <- geno$genotypes
  CHR <- geno$map$chromosome
  infos <- snp_fastImputeSimple(G, method = "mean2")
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

### 0.2 Prepare for GCTA (GRM)

```r
ref.bed = "../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01"
MakeGRM = function(ref.bed, keep.snp, keep.ind, output.dir = "./temp/"){

  grm = paste0(output.dir, "geno")

  fwrite(data.frame(keep.snp), file = paste0(output.dir,"snplist.txt"), quote = F, sep = "\t", col.names=T, na = "NA")

  system(paste0("plink --bfile ", ref.bed, " --extract ", paste0(output.dir, "snplist.txt"), " --make-grm-bin --out ", grm ) )
}
```

### 0.3 Prepare LD score for LDSC

```bash
ldsc.py --bfile ../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01 --extract ./temp/snplist.txt --ld-wind-kb 1000 --out ./temp/geno
```

## 1. GCTA

GCTA is widely used for GWAS analysis, but it can be awkward to integrate into an R-based simulation pipeline. Below is a minimal workflow for running bivariate REML from R.

**Inputs**
- phenotype file (Y1, Y2)
- genotype file (geno.bed)

First, generate the GRM from the genotype file. This step only needs to be done once.

```r
MakeGRM = function(ref.bed, keep.snp, keep.ind, output.dir = "./temp/"){

  grm = paste0(output.dir, "geno")

  fwrite(data.frame(keep.snp), file = paste0(output.dir,"snplist.txt"), quote = F, sep = "\t", col.names=T, na = "NA")

  system(paste0("plink --bfile ", ref.bed, " --extract ", paste0(output.dir, "snplist.txt"), " --make-grm-bin --out ", grm ) )
}
```

Next, write the phenotype file using FID and IID.

```r
WritePheno = function( Y1, Y2, FID, IID, pheno.file  ){

  pheno = data.frame(FID = FID, IID = IID,
                     Y1 = Y1, Y2 = Y2)
  
  fwrite(pheno, file = pheno.file,
         quote = F, sep = "\t", col.names=T, na = "NA")
}
```

Then run GCTA.


```r
SimuGCTA = function(geno.file, pheno.file, output.dir){

  output.file = paste0(output.dir, "GCTA.result")

  system(paste0("./gcta64 --reml-bivar --grm ",
                geno.file," --pheno ", pheno.file,
                " --thread-num 20 --out ", output.file ))
}
```

If simulated errors have no covariance, use `--reml-bivar-nocove` for a fair comparison. To test non-zero genetic covariance, use `--reml-bivar-lrt-rg 0`.


```r
SimuGCTA.test = function(geno.file, pheno.file, output.dir){

  output.file = paste0(output.dir, "GCTA.result")
  
  system(paste0("./gcta64 --reml-bivar-nocove --reml-maxit 100 --grm ", 
                geno.file," --pheno ",
                pheno.file," --reml-bivar-lrt-rg 0 --out ",output.file))
}
```

Read the GCTA output file.


```r
read.GCTA = function(output.dir, LRT = T){
  
   filename = paste0(output.dir, "GCTA.result.hsq")
    if(file.exists(filename)){
      hsq = data.table::fread(filename, fill = T)
      Est = as.vector(t( hsq[1:3, 2:3])) %>% 
    setNames(c("Heri.1","Heri.1.Se","Heri.2","Heri.2.Se", "GeCv","GeCv.Se" ))
      if(LRT){
      pval = as.numeric(sapply(strsplit(as.character(hsq[hsq$Source == "Pval",2]), "\\("), "[[", 1))
      }
    }
   return(list(Est = Est, pval = pval))
}
```


## 2. LDSC

First, compute LD scores from the genotype file (done once).

```bash
module load ldsc
source activate ldsc-1.0.1
ldsc.py --bfile geno --ld-wind-cm 1 --out chr1_select
ldsc.py --bfile ../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01 --extract ./temp/snplist.txt --ld-wind-kb 1000 --out ./temp/geno
```

Implement LDSC regression.

```r
SimuLDSC = function(Z1, N1, Z2, N2, Nc=0, weight=T, CSTR = T,
                           info = NA, LDSCORE, corenum = 1, output.dir){
  
  if(length(Z1) != length(Z2) ){stop("Z scores should have the same length.")}

  if(!is.na(info)){
    print("update info")
    Z1.data = data.frame(SNP = info$SNP , Z = Z1, N =  N1,
                         A1 =  "A", A2 = "G", stringsAsFactors = F)
    colnames(Z1.data) = c("SNP",  "Z", "N", "A1", "A2")
    fwrite(Z1.data, file = paste0(output.dir, "Z1.sumstats") , quote = F, sep = "\t", col.names=T, na = "NA" )
    Z2.data = data.frame(SNP = info$SNP , Z = Z2, N =  N2,
                         A1 = "A",  A2 = "G", stringsAsFactors = F)
    colnames(Z2.data) = c("SNP", "Z", "N", "A1", "A2")
    fwrite(Z2.data, file = paste0(output.dir, "Z2.sumstats"), quote = F, sep = "\t", col.names=T, na = "NA" )
  }
  if(weight == T){
    name = paste0(output.dir, "LDSC")
    system(paste0("ldsc.py --rg ",paste0(J, "Z1.sumstats,"), paste0(J, "Z2.sumstats"),
                  " --ref-ld ", LDSCORE," --w-ld ",LDSCORE," --intercept-h2 1,1 --out ",name))
  }else{

    stop()

  } 
  }
```


Read the LDSC output.


```r
read.LDSC = function(output.dir){
  H = fread(paste0(output.dir, "LDSC.log"), fill = T)
  result.all = rep(NA, 8)
  result.gecr = unlist(H[(stringr::str_detect(unlist(H), "Genetic Correlation:"))]) %>% str_split(, pattern = " ") %>% unlist()
  result.heri = unlist(H[(stringr::str_detect(unlist(H), "Total Observed scale h2:"))]) %>%  str_split(, pattern = " ") %>% unlist()
  result.gecov = unlist(H[(stringr::str_detect(unlist(H), "Total Observed scale gencov:"))]) %>% str_split(, pattern = " ") %>% unlist()
  result.P = unlist(H[(stringr::str_detect(unlist(H), "P:"))]) %>%
    str_split(, pattern = " ") %>% unlist()

  if(!is.null(result.gecr)){
    result =  as.numeric(unlist(regmatches(result.gecr,gregexpr("-?\\ *[0-9]+\\.?[0-9]*(?:[Ee]\\ *-?\\ *[0-9]+)?",result.gecr, perl=TRUE))))
    result.all[1] = result[1]
    result.all[2] = result[2]
  }

  if(!is.null(result.heri)){
    result =  as.numeric(unlist(regmatches(result.heri,
                                           gregexpr("-?\\ *[0-9]+\\.?[0-9]*(?:[Ee]\\ *-?\\ *[0-9]+)?",
                                                    result.heri, perl=TRUE))))
    result.all[3] = result[2]
    result.all[4] = result[3]
    result.all[5] = result[5]
    result.all[6] = result[6]
  }

  if(!is.null(result.gecov)){

    result =  as.numeric(unlist(regmatches(result.gecov,
                                           gregexpr("-?\\ *[0-9]+\\.?[0-9]*(?:[Ee]\\ *-?\\ *[0-9]+)?",
                                                    result.gecov, perl=TRUE))))
    result.all[7] = result[1]
    result.all[8] = result[2]
    zscore = result[1]/result[2]
    pvalue = 2* (1 - pnorm(abs(zscore)))
  }else{
    pvalue = NA
  }
  return( c(GeCv.pv = pvalue, GeCr.pv = as.numeric(result.P[2])) )
}
```

