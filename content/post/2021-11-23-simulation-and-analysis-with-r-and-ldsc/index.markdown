---
title: Simulation and analysis with R and LDSC
author: Jianqiao Wang
date: '2021-11-23'
slug: simulation-and-analysis-with-r-and-ldsc
categories: []
tags: []
subtitle: ''
summary: ''
authors: []
lastmod: '2021-11-23T11:50:44-05:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---



First, we generate the ld score based on the genotype file. This function only need be estimated once 

module load ldsc
source activate ldsc-1.0.1
ldsc.py --bfile geno  --ld-wind-cm 1 --out chr1_select
ldsc.py --bfile ../../Cric_HF_RealData/rawdata/cric.filtered.maf0.01 --extract ./temp/snplist.txt --ld-wind-kb 1000 --out ./temp/geno


Implement the LDSC regression 

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


read the LDSC output 


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

