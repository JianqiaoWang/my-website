---
title: Simulation and analysis with R and GCTA
author: ''
date: '2021-11-22'
slug: simulation-and-analysis-with-r-and-gcta
categories: []
tags: [software]
subtitle: ''
summary: ''
authors: []
lastmod: '2021-11-22T22:45:32-05:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---

My dissertation research always needs to do comparison with the GCTA-GREML method. GCTA is a wonderful software designed for the dealing with the GWAS data. On the other hand, it does not always fit in the common  statistical simulation pipeline, particualy for the software R. This post is used to document how to do statistical simulation with R and GCTA based on a give genotype file.

What we need:
- phenotype file Y1, Y2
- Genotype file  geno.bed

First, we generate the GRM file based on the genotype file. This function only need be estimated once 

```r
MakeGRM = function(ref.bed, keep.snp, keep.ind, output.dir = "./temp/"){

  grm = paste0(output.dir, "geno")

  fwrite(data.frame(keep.snp), file = paste0(output.dir,"snplist.txt"), quote = F, sep = "\t", col.names=T, na = "NA")

  system(paste0("plink --bfile ", ref.bed, " --extract ", paste0(output.dir, "snplist.txt"), " --make-grm-bin --out ", grm ) )
}
```

In the second step, we need to write out the phenotype file based on the FID and IID 

```r
WritePheno = function( Y1, Y2, FID, IID, pheno.file  ){

  pheno = data.frame(FID = FID, IID = IID,
                     Y1 = Y1, Y2 = Y2)
  
  fwrite(pheno, file = pheno.file,
         quote = F, sep = "\t", col.names=T, na = "NA")
}
```

Then we implement the GCTA, 


```r
SimuGCTA = function(geno.file, pheno.file, output.dir){

  output.file = paste0(output.dir, "GCTA.result")

  system(paste0("./gcta64 --reml-bivar --grm ",
                geno.file," --pheno ", pheno.file,
                " --thread-num 20 --out ", output.file ))
}
```

When the simulated errors have no covariance, we could specify --reml-bivar-nocove for fair comparison. To test for the non-zero genetic covariance, ** --reml-bivar-lrt-rg 0 ** could be used.


```r
SimuGCTA.test = function(geno.file, pheno.file, output.dir){

  output.file = paste0(output.dir, "GCTA.result")
  
  system(paste0("./gcta64 --reml-bivar-nocove --reml-maxit 100 --grm ", 
                geno.file," --pheno ",
                pheno.file," --reml-bivar-lrt-rg 0 --out ",output.file))
}
```

Read GCTA output file 


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
