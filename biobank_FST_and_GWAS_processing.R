## Packages####
library(ggplot2)
library(gridExtra)
library(reshape2)
library(data.table)
library(CMplot)
library(Hmisc)
library(cowplot)
library(doParallel)

## Functions
nei_fst_func <- function(pm,pf){((pm-pf)^2)/(4*((pm+pf)/2)*(1-((pm+pf)/2)))}

geno_to_afreq_func <- function(hom_ref_ct,het_ct,hom_alt_ct){
  ((het_ct*0.5)+hom_alt_ct)/(hom_ref_ct+het_ct+hom_alt_ct)
}
geno_to_afreq_haploid_func <- function(ref_ct,alt_ct){
  alt_ct/(ref_ct+alt_ct)
}

fis_func_one_sex_only <- function(p,pAa){ ( pAa/(2*p*(1-p)) ) - 1 }
fis_func <- function(pm,pf,pAa){ ( pAa/(2*((pm+pf)/2)*(1-((pm+pf)/2))) ) - 1 }

bt <- function(x,n,p) {binom.test(x,n,p, alternative = c("less"))$p.value}

t.test.p_func <- function(m1,m2,se1,se2,n1,n2,m0=0,rho)
{
  #se_denom <- sqrt(  se1^2 + se2^2 )
  se_dep_denom <- sqrt( se1^2 + se2^2 - ( 2 * rho * se1 * se2 ) )
  
  #welch-satterthwaite df
  df <- ( (se1^2 + se2^2)^2 )/( (se1^2)/(n1-1) + (se2^2)/(n2-1) )
     
  #t <- (m1-m2-m0)/se_denom
  t_dep <- (m1-m2-m0)/se_dep_denom
  
  #2*pt(abs(t),df,lower.tail = F)
  
  dat <- c(2*pt(abs(t_dep),df,lower.tail = F))   
  names(dat) <- c("T_DEP_P")
  return(dat) 
}

InvTable = function(tb){
  output = rep(names(tb), tb)
  return(output)
}


dir <- "~/Dropbox/dropbox_work/data/biobank/"


## S-S VIABILITY / REPROD / TOTAL####


##Pm, pf, pm', pf' (observed) ####


## Genotype data, autosomes
for (j in 1:22){
  print(paste(j,"initiated"))
  ## Frequencies, females
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",j,"_sampleqc1_snpqc1b.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c(1:4)]
  #Total genotype counts for adult females (with 0 #offspring)
  f.frq$CT_TOT <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for females (with 0 #offspring)
  f.frq$HET_F <- f0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for females (with 0 #offspring)
  f.frq$HOM_ALT_F <- f0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for females (with 0 #offspring)
  f.frq$M_HOM_REF <- f0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for females (with 0 #offspring)
  f.frq$M_HET <- f0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for females (with 0 #offspring)
  f.frq$M_HOM_ALT <- f0.frq$TWO_ALT_GENO_CTS*0
  rm(f0.frq)
  
  for (i in c(1:14,19)){
    f1.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",j,"_sampleqc1_snpqc1b.Female",i,".gcount"))
    #Total genotype counts for adult females (with i #offspring), added to adult females with 0 #offspring
    f.frq$CT_TOT <- f.frq$CT_TOT+(f1.frq$HOM_REF_CT+f1.frq$HET_REF_ALT_CTS+f1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HET_F <- f.frq$HET_F+(f1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HOM_ALT_F <- f.frq$HOM_ALT_F+(f1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_REF <- f.frq$M_HOM_REF+(f1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HET <- f.frq$M_HET+(f1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_ALT <- f.frq$M_HOM_ALT+(f1.frq$TWO_ALT_GENO_CTS*i)
  }

  ## Frequencies, males
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",j,"_sampleqc1_snpqc1b.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c(1:4)]
  #Total genotype counts for adult males (with 0 #offspring)
  m.frq$CT_TOT <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for males (with 0 #offspring)
  m.frq$HET_M <- m0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for males (with 0 #offspring)
  m.frq$HOM_ALT_M <- m0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for males (with 0 #offspring)
  m.frq$M_HOM_REF <- m0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for males (with 0 #offspring)
  m.frq$M_HET <- m0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for males (with 0 #offspring)
  m.frq$M_HOM_ALT <- m0.frq$TWO_ALT_GENO_CTS*0
  rm(m0.frq)
  
  for (i in c(1:15,17)){
    m1.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",j,"_sampleqc1_snpqc1b.Male",i,".gcount"))
    #Total genotype counts for adult males (with i #offspring), added to adult males with 0 #offspring
    m.frq$CT_TOT <- m.frq$CT_TOT+(m1.frq$HOM_REF_CT+m1.frq$HET_REF_ALT_CTS+m1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HET_M <- m.frq$HET_M+(m1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HOM_ALT_M <- m.frq$HOM_ALT_M+(m1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_REF <- m.frq$M_HOM_REF+(m1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HET <- m.frq$M_HET+(m1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_ALT <- m.frq$M_HOM_ALT+(m1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  #Combine total allele frequency counts and 
  mf.frq <- cbind(f.frq,m.frq[,c(5:10)])
  names(mf.frq)[c(5:16)] <-c("CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M")
  #Reference allele frequency in adult females (Pf)
  mf.frq$p_F_VIABILITY <- 1-(mf.frq$HOM_ALT_F+(0.5*mf.frq$HET_F))/(mf.frq$CT_TOT_F)
  #Reference allele frequency in adult males (Pm)
  mf.frq$p_M_VIABILITY <- 1-(mf.frq$HOM_ALT_M+(0.5*mf.frq$HET_M))/(mf.frq$CT_TOT_M)
  #Adult Fst
  mf.frq$ADULT_FST <- nei_fst_func(mf.frq$p_F_VIABILITY,mf.frq$p_M_VIABILITY)
  #Projected reference allele frequency in adult females (Pf')
  mf.frq$p_F <- (mf.frq$M_HOM_REF_F+(0.5*mf.frq$M_HET_F))/(mf.frq$M_HOM_REF_F+mf.frq$M_HET_F+mf.frq$M_HOM_ALT_F) 
  #Projected reference allele frequency in adult females (Pm')
  mf.frq$p_M <- (mf.frq$M_HOM_REF_M+(0.5*mf.frq$M_HET_M))/(mf.frq$M_HOM_REF_M+mf.frq$M_HET_M+mf.frq$M_HOM_ALT_M)
  #Gametic Fst
  mf.frq$GAMETIC_FST <- nei_fst_func(mf.frq$p_F,mf.frq$p_M)
  
  write.table(mf.frq,paste0(dir,"ukb_fst_v3/ukb52049_chr",j,"_sampleqc1_snpqc1b_no_further_maf_filter.Fst"),quote=F,row.names=F,sep="\t",col.names=F)
  
  print(paste(j,"completed"))
}


## Imputed data, autosomes
for (j in 1:22){
  print(paste(j,"initiated"))
  ## Frequencies, females
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c(1:4)]
  #Total genotype counts for adult females (with 0 #offspring)
  f.frq$CT_TOT <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for females (with 0 #offspring)
  f.frq$HET_F <- f0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for females (with 0 #offspring)
  f.frq$HOM_ALT_F <- f0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for females (with 0 #offspring)
  f.frq$M_HOM_REF <- f0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for females (with 0 #offspring)
  f.frq$M_HET <- f0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for females (with 0 #offspring)
  f.frq$M_HOM_ALT <- f0.frq$TWO_ALT_GENO_CTS*0
  rm(f0.frq)
  
  for (i in c(1:14,19)){
    f1.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b.Female",i,".gcount"))
    #Total genotype counts for adult females (with i #offspring), added to adult females with 0 #offspring
    f.frq$CT_TOT <- f.frq$CT_TOT+(f1.frq$HOM_REF_CT+f1.frq$HET_REF_ALT_CTS+f1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HET_F <- f.frq$HET_F+(f1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HOM_ALT_F <- f.frq$HOM_ALT_F+(f1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_REF <- f.frq$M_HOM_REF+(f1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HET <- f.frq$M_HET+(f1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_ALT <- f.frq$M_HOM_ALT+(f1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  ## Frequencies, males
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c(1:4)]
  #Total genotype counts for adult males (with 0 #offspring)
  m.frq$CT_TOT <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for males (with 0 #offspring)
  m.frq$HET_M <- m0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for males (with 0 #offspring)
  m.frq$HOM_ALT_M <- m0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for males (with 0 #offspring)
  m.frq$M_HOM_REF <- m0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for males (with 0 #offspring)
  m.frq$M_HET <- m0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for males (with 0 #offspring)
  m.frq$M_HOM_ALT <- m0.frq$TWO_ALT_GENO_CTS*0
  rm(m0.frq)
  
  for (i in c(1:15,17)){
    m1.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b.Male",i,".gcount"))
    #Total genotype counts for adult males (with i #offspring), added to adult males with 0 #offspring
    m.frq$CT_TOT <- m.frq$CT_TOT+(m1.frq$HOM_REF_CT+m1.frq$HET_REF_ALT_CTS+m1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HET_M <- m.frq$HET_M+(m1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HOM_ALT_M <- m.frq$HOM_ALT_M+(m1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_REF <- m.frq$M_HOM_REF+(m1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HET <- m.frq$M_HET+(m1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_ALT <- m.frq$M_HOM_ALT+(m1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  #Combine total allele frequency counts and 
  mf.frq <- cbind(f.frq,m.frq[,c(5:10)])
  names(mf.frq)[c(5:16)] <-c("CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M")
  #Reference allele frequency in adult females (Pf)
  mf.frq$p_F_VIABILITY <- 1-(mf.frq$HOM_ALT_F+(0.5*mf.frq$HET_F))/(mf.frq$CT_TOT_F)
  #Reference allele frequency in adult males (Pm)
  mf.frq$p_M_VIABILITY <- 1-(mf.frq$HOM_ALT_M+(0.5*mf.frq$HET_M))/(mf.frq$CT_TOT_M)
  #Adult Fst
  mf.frq$ADULT_FST <- nei_fst_func(mf.frq$p_F_VIABILITY,mf.frq$p_M_VIABILITY)
  #Projected reference allele frequency in adult females (Pf')
  mf.frq$p_F <- (mf.frq$M_HOM_REF_F+(0.5*mf.frq$M_HET_F))/(mf.frq$M_HOM_REF_F+mf.frq$M_HET_F+mf.frq$M_HOM_ALT_F) 
  #Projected reference allele frequency in adult females (Pm')
  mf.frq$p_M <- (mf.frq$M_HOM_REF_M+(0.5*mf.frq$M_HET_M))/(mf.frq$M_HOM_REF_M+mf.frq$M_HET_M+mf.frq$M_HOM_ALT_M)
  #Gametic Fst
  mf.frq$GAMETIC_FST <- nei_fst_func(mf.frq$p_F,mf.frq$p_M)
  
  
  write.table(mf.frq,paste0(dir,"ukb_fst_v3/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_no_info_filter_no_further_maf_filter.Fst"),quote=F,row.names=F,sep="\t",col.names=F)
  print(paste(j,"completed"))
}

#Mean/var LRS (observed)####

##Genotyped
#Female LRS
f.frq.all.chrom <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c("CHROM","ID","REF","ALT")]
  f.frq$NOT_MISSING_CT_0 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female1.gcount"))
  f.frq$NOT_MISSING_CT_1 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female2.gcount"))
  f.frq$NOT_MISSING_CT_2 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female3.gcount"))
  f.frq$NOT_MISSING_CT_3 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female4.gcount"))
  f.frq$NOT_MISSING_CT_4 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female5.gcount"))
  f.frq$NOT_MISSING_CT_5 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female6.gcount"))
  f.frq$NOT_MISSING_CT_6 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female7.gcount"))
  f.frq$NOT_MISSING_CT_7 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female8.gcount"))
  f.frq$NOT_MISSING_CT_8 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female9.gcount"))
  f.frq$NOT_MISSING_CT_9 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female10.gcount"))
  f.frq$NOT_MISSING_CT_10 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female11.gcount"))
  f.frq$NOT_MISSING_CT_11 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female12.gcount"))
  f.frq$NOT_MISSING_CT_12 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female13.gcount"))
  f.frq$NOT_MISSING_CT_13 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female14.gcount"))
  f.frq$NOT_MISSING_CT_14 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Female19.gcount"))
  f.frq$NOT_MISSING_CT_19 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f.frq.all.chrom[[i]] <- f.frq
  print(i)
}
f.frq.all.chrom <- do.call(rbind,f.frq.all.chrom)

#Males
m.frq.all.chrom <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c("CHROM","ID","REF","ALT")]
  m.frq$NOT_MISSING_CT_0 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male1.gcount"))
  m.frq$NOT_MISSING_CT_1 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male2.gcount"))
  m.frq$NOT_MISSING_CT_2 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male3.gcount"))
  m.frq$NOT_MISSING_CT_3 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male4.gcount"))
  m.frq$NOT_MISSING_CT_4 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male5.gcount"))
  m.frq$NOT_MISSING_CT_5 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male6.gcount"))
  m.frq$NOT_MISSING_CT_6 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male7.gcount"))
  m.frq$NOT_MISSING_CT_7 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male8.gcount"))
  m.frq$NOT_MISSING_CT_8 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male9.gcount"))
  m.frq$NOT_MISSING_CT_9 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male10.gcount"))
  m.frq$NOT_MISSING_CT_10 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male11.gcount"))
  m.frq$NOT_MISSING_CT_11 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male12.gcount"))
  m.frq$NOT_MISSING_CT_12 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male13.gcount"))
  m.frq$NOT_MISSING_CT_13 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male14.gcount"))
  m.frq$NOT_MISSING_CT_14 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male15.gcount"))
  m.frq$NOT_MISSING_CT_15 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr",i,"_sampleqc1_snpqc1b.Male17.gcount"))
  m.frq$NOT_MISSING_CT_17 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m.frq.all.chrom[[i]] <- m.frq
  print(i)
}
m.frq.all.chrom <- do.call(rbind,m.frq.all.chrom)

#Variance and mean^2 fitness for each locus, females
#Column names = LRS
names(f.frq.all.chrom)[5:20] <- as.character(c(0:14,19))
f.frq.all.chrom2 <- f.frq.all.chrom[,5:20]
#Histogram bins
bins <- as.numeric(names(f.frq.all.chrom2))
#Total number of LRS values across bins (n)
n <- rowSums(f.frq.all.chrom2)
#Mean = sum(bin*LRS-per-bin)/n
meanf.geno <- vector()
meanf.geno <- foreach (i=1:nrow(f.frq.all.chrom2), .combine='c') %dopar% {
  (sum(bins*f.frq.all.chrom2[i])/n[i])
}
#Var = sum((bin-mean)^2*LRS-per-bin)/n
varf.geno <- vector()
varf.geno <- foreach (i=1:nrow(f.frq.all.chrom2), .combine='c') %dopar% {
  sum((bins-meanf.geno[i])^2*f.frq.all.chrom2[i])/n[i]
}
#Mean^2 
meansqf.geno <- meanf.geno^2

#Variance and mean^2 fitness for each locus, males
#Column names = LRS
names(m.frq.all.chrom)[5:21] <- as.character(c(0:15,17))
m.frq.all.chrom2 <- m.frq.all.chrom[,5:21]
#Histogram bins
bins <- as.numeric(names(m.frq.all.chrom2))
#Total number of LRS values across bins (n)
n <- rowSums(m.frq.all.chrom2)
#Mean = sum(bin*LRS-per-bin)/n
meanm.geno <- vector()
for(i in 1:nrow(m.frq.all.chrom2)){
  meanm.geno[i] <- (sum(bins*m.frq.all.chrom2[i])/n[i])
}
#Var = sum((bin-mean)^2*LRS-per-bin)/n
varm.geno <- vector()
for(i in 1:nrow(m.frq.all.chrom2)){
  varm.geno[i] <- sum((bins-meanm.geno[i])^2*m.frq.all.chrom2[i])/n[i]
}
#Mean^2 
meansqm.geno <- meanm.geno^2


mf.frq.all.chrom <- cbind(f.frq.all.chrom[,1:4],varf.geno,meansqf.geno,varm.geno,meansqm.geno)
names(mf.frq.all.chrom)[5:8] <- c("VAR_LRS_F","MEANSQ_LRS_F","VAR_LRS_M","MEANSQ_LRS_M")

#write.table(mf.frq.all.chrom,paste0(dir,"/ukb_locus_by_locus_fitness_var_and_mean_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b.var_and_mean.txt"),quote=F,row.names=F,sep="\t")



##Imputed

#CHROM/ID/REF/ALT data for MAF + artefact-filtered sites
good.ids <- fread(paste0(dir,"/ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.4columns"))

registerDoParallel(cores = 4)
#Female LRS
f.frq.chrom.pos.imp <- list()
meansqf.imp <- list()
varf.imp <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c("CHROM","ID","REF","ALT")]
  f.frq$NOT_MISSING_CT_0 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female1.gcount"))
  f.frq$NOT_MISSING_CT_1 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female2.gcount"))
  f.frq$NOT_MISSING_CT_2 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female3.gcount"))
  f.frq$NOT_MISSING_CT_3 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female4.gcount"))
  f.frq$NOT_MISSING_CT_4 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female5.gcount"))
  f.frq$NOT_MISSING_CT_5 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female6.gcount"))
  f.frq$NOT_MISSING_CT_6 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female7.gcount"))
  f.frq$NOT_MISSING_CT_7 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female8.gcount"))
  f.frq$NOT_MISSING_CT_8 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female9.gcount"))
  f.frq$NOT_MISSING_CT_9 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female10.gcount"))
  f.frq$NOT_MISSING_CT_10 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female11.gcount"))
  f.frq$NOT_MISSING_CT_11 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female12.gcount"))
  f.frq$NOT_MISSING_CT_12 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female13.gcount"))
  f.frq$NOT_MISSING_CT_13 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female14.gcount"))
  f.frq$NOT_MISSING_CT_14 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Female19.gcount"))
  f.frq$NOT_MISSING_CT_19 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  print(i)
  
  #Keep only MAF/artefact-filtered sites
  f.frq <- subset(f.frq,ID %in% good.ids$ID)
  #Variance and mean^2 fitness for each locus
  #Column names = LRS
  names(f.frq)[5:20] <- as.character(c(0:14,19))
  f.frq2 <- f.frq[,5:20]
  #Histogram bins
  bins <- as.numeric(names(f.frq2))
  #Total number of LRS values across bins (n)
  n <- rowSums(f.frq2)
  #Mean = sum(bin*LRS-per-bin)/n
  meansqf.imp[[i]] <- foreach (j=1:nrow(f.frq2), .combine='rbind') %dopar% {
    (sum(bins*f.frq2[j])/n[j])^2
  }
  #Var = sum((bin-mean)^2*LRS-per-bin)/n
  varf.imp[[i]] <- foreach (j=1:nrow(f.frq2), .combine='rbind') %dopar% {
    sum((bins-sqrt(meansqf.imp[[i]][j]))^2*f.frq2[j])/n[j]
  }
  #Positional information
  f.frq.chrom.pos.imp[[i]] <- f.frq[,1:4]
}

meansqf.imp <- unlist(meansqf.imp)
varf.imp <- unlist(varf.imp)
f.frq.chrom.pos.imp <- do.call(rbind,f.frq.chrom.pos.imp)


#Males
meansqm.imp <- list()
varm.imp <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c("CHROM","ID","REF","ALT")]
  m.frq$NOT_MISSING_CT_0 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male1.gcount"))
  m.frq$NOT_MISSING_CT_1 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male2.gcount"))
  m.frq$NOT_MISSING_CT_2 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male3.gcount"))
  m.frq$NOT_MISSING_CT_3 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male4.gcount"))
  m.frq$NOT_MISSING_CT_4 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male5.gcount"))
  m.frq$NOT_MISSING_CT_5 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male6.gcount"))
  m.frq$NOT_MISSING_CT_6 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male7.gcount"))
  m.frq$NOT_MISSING_CT_7 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male8.gcount"))
  m.frq$NOT_MISSING_CT_8 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male9.gcount"))
  m.frq$NOT_MISSING_CT_9 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male10.gcount"))
  m.frq$NOT_MISSING_CT_10 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male11.gcount"))
  m.frq$NOT_MISSING_CT_11 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male12.gcount"))
  m.frq$NOT_MISSING_CT_12 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male13.gcount"))
  m.frq$NOT_MISSING_CT_13 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male14.gcount"))
  m.frq$NOT_MISSING_CT_14 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male15.gcount"))
  m.frq$NOT_MISSING_CT_15 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b.Male17.gcount"))
  m.frq$NOT_MISSING_CT_17 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  print(i)
  
  #Keep only MAF/artefact-filtered sites
  m.frq <- subset(m.frq,ID %in% good.ids$ID)
  
  #Variance and mean^2 fitness for each locus
  #Column names = LRS
  names(m.frq)[5:21] <- as.character(c(0:15,17))
  m.frq2 <- m.frq[,5:21]
  #Histogram bins
  bins <- as.numeric(names(m.frq2))
  #Total number of LRS values across bins (n)
  n <- rowSums(m.frq2)
  #Mean = sum(bin*LRS-per-bin)/n
  meansqm.imp[[i]] <- foreach (j=1:nrow(m.frq2), .combine='rbind') %dopar% {
    (sum(bins*m.frq2[j])/n[j])^2
  }
  #Var = sum((bin-mean)^2*LRS-per-bin)/n
  varm.imp[[i]] <- foreach (j=1:nrow(m.frq2), .combine='rbind') %dopar% {
    sum((bins-sqrt(meansqm.imp[[i]][j]))^2*m.frq2[j])/n[j]
  }
}
meansqm.imp <- unlist(meansqm.imp)
varm.imp <- unlist(varm.imp)


mf.frq.all.chrom <- cbind(f.frq.chrom.pos.imp[,1:4],varf.imp,meansqf.imp,varm.imp,meansqm.imp)
names(mf.frq.all.chrom)[5:8] <- c("VAR_LRS_F","MEANSQ_LRS_F","VAR_LRS_M","MEANSQ_LRS_M")

#write.table(mf.frq.all.chrom,paste0(dir,"/ukb_locus_by_locus_fitness_var_and_mean_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.var_and_mean.txt"),quote=F,row.names=F,sep="\t")



## Filter INFO score (imputed data)####

#Create lists with INFO score > 0.8, for each chromosome
chromosome_vector <- c(1:22,"X","XY")
for (i in chromosome_vector){
  imp.info <- fread(paste0(dir,"ukb_imp_mfi/ukb_mfi_chr",i,"_v3.txt"))
  imp.goodinfo <- subset(imp.info,V8>0.8)
  write.table(imp.goodinfo$V2,paste0(dir,"ukb_imp_mfi/ukb_imp_chr",i,"_v3_info_above_0.8.ids"),quote=F,col.names=F,row.names=F)
  print(i)
}


## Filter MAF####

#Only needed for observed (Unix code will be used to filter MAF in permuted data)

## Genotype data, observed
mf.frq <- fread(paste0(dir,"ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_no_further_maf_filter.Fst")) 
names(mf.frq) <- c("CHROM","ID","REF","ALT","CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M","p_F_VIABILITY","p_M_VIABILITY","ADULT_FST","p_F","p_M","GAMETIC_FST")
mf.frq$MAF <- with(mf.frq,ifelse(((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F)>0.5,1-((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F),((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F))) 
write.table(subset(mf.frq,MAF>0.01),paste0(dir,"ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01.Fst"),row.names=F,quote=F,sep="\t")
write.table(subset(mf.frq,MAF>0.01)[,2],paste0(dir,"ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01.ids"),row.names=F,quote=F,sep="\t")


## Imputed data, observed 
mf.frq <- fread(paste0(dir,"ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_no_further_maf_filter.Fst")) 
names(mf.frq) <- c("CHROM","ID","REF","ALT","CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M","p_F_VIABILITY","p_M_VIABILITY","ADULT_FST","p_F","p_M","GAMETIC_FST")
mf.frq$MAF <- with(mf.frq,ifelse(((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F)>0.5,1-((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F),((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F))) 
write.table(subset(mf.frq,MAF>0.01),paste0(dir,"ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01.Fst"),row.names=F,quote=F,sep="\t")
write.table(subset(mf.frq,MAF>0.01)[,2],paste0(dir,"ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01.ids"),row.names=F,quote=F,sep="\t")



## Filter artefacts####

#Only needed for observed (Unix code will be used to filter artefacts in permuted)


## Genotype data, autosomes, MAF>0.01

mf.frq <- fread(paste0(dir,"ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01.Fst")) 

#Necessary columns for subsequent functions
mf.frq$ALT_FREQ <- 1-with(mf.frq,((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F))
mf.frq$ALT_FREQ_F <- 1-mf.frq$p_F_VIABILITY
mf.frq$ALT_FREQ_M <- 1-mf.frq$p_M_VIABILITY
mf.frq$HET <- mf.frq$HET_M+mf.frq$HET_F
mf.frq$HOM_ALT <- mf.frq$HOM_ALT_M+mf.frq$HOM_ALT_F
mf.frq$CT_TOT <- mf.frq$CT_TOT_M+mf.frq$CT_TOT_F
mf.frq$MAF_F <- with(mf.frq,ifelse(((p_F_VIABILITY*CT_TOT_F))/(CT_TOT_F)>0.5,1-((p_F_VIABILITY*CT_TOT_F))/(CT_TOT_F),((p_F_VIABILITY*CT_TOT_F))/(CT_TOT_F)))
mf.frq$MAF_M <- with(mf.frq,ifelse(((p_M_VIABILITY*CT_TOT_M))/(CT_TOT_M)>0.5,1-((p_M_VIABILITY*CT_TOT_M))/(CT_TOT_M),((p_M_VIABILITY*CT_TOT_M))/(CT_TOT_M)))
mf.frq$MISSING_M <- 115531-mf.frq$CT_TOT_M
mf.frq$MISSING_F <- 133490-mf.frq$CT_TOT_F 

#1. Filter for excess heterozygotes (elevated Fis), assuming s=0.2 and at equilibrium
mf.frq$FIS <- fis_func(mf.frq$ALT_FREQ_M,mf.frq$ALT_FREQ_F,(mf.frq$HET)/(mf.frq$CT_TOT))
mf.frq$FIS_M <- fis_func_one_sex_only(mf.frq$ALT_FREQ_M,(mf.frq$HET_M)/(mf.frq$CT_TOT_M))
mf.frq$FIS_F <- fis_func_one_sex_only(mf.frq$ALT_FREQ_F,(mf.frq$HET_F)/(mf.frq$CT_TOT_F))
#Excess heterozygosity when s=0.2 is as follows (previous code)
#s <- 0.2
#excess_het <- ( ( (2-s)/2 * sqrt(1+(s/(2-s))^2) ) - (2-s)/2 )^2
#Based on this, and given that Mean Fis under the null is 1/2n, where n is number of individuals, and Variance Fis under the null is 1/n, we obtain q-values for Fis in both sexes 
#mf.frq$FIS_P <- p.adjust(pnorm(mf.frq$FIS,mean=(1/(2*(mf.frq$CT_TOT_M+mf.frq$CT_TOT_F)))+excess_het,sd=sqrt(1/(mf.frq$CT_TOT_M+mf.frq$CT_TOT_F)),lower.tail = F),"BH")
#males
#mf.frq$FIS_M_P <- p.adjust(pnorm(mf.frq$FIS_M,mean=(1/(2*(mf.frq$CT_TOT_M)))+excess_het,sd=sqrt(1/(mf.frq$CT_TOT_M)),lower.tail = F),"BH")
#females
#mf.frq$FIS_F_P <- p.adjust(pnorm(mf.frq$FIS_F,mean=(1/(2*(mf.frq$CT_TOT_F)))+excess_het,sd=sqrt(1/(mf.frq$CT_TOT_F)),lower.tail = F),"BH")

#Updated check for excess heterozygosity
s <- 0.2
mf.frq$FIS_P <- p.adjust(pnorm(mf.frq$FIS,mean=(1/(2*mf.frq$CT_TOT)) + (mf.frq$MAF*(1-mf.frq$MAF)/4)*(s/(1 - mf.frq$MAF*s))^2,sd=sqrt(1/(mf.frq$CT_TOT)),lower.tail = F),"BH")
mf.frq$FIS_M_P <- p.adjust(pnorm(mf.frq$FIS_M,mean=(1/(2*mf.frq$CT_TOT_M)) + (mf.frq$MAF_M*(1-mf.frq$MAF_M)/4*(s/(1 - mf.frq$MAF_M*s))^2),sd=sqrt(1/(mf.frq$CT_TOT_M)),lower.tail = F),"BH")
mf.frq$FIS_F_P <- p.adjust(pnorm(mf.frq$FIS_F,mean=(1/(2*mf.frq$CT_TOT_F)) + (mf.frq$MAF_F*(1-mf.frq$MAF_F)/4*(s/(1 - mf.frq$MAF_F*s))^2),sd=sqrt(1/(mf.frq$CT_TOT_F)),lower.tail = F),"BH")

#2. Filter based on deficit of minor allele homozygotes (as per Kasimatis et al. 2020)
#bt function is a one-sided binomial test test that asks if the number of homozygotes is less than expected

#females: If alternative allele is minor, binomial test for deficit of alternative allele homozygotes
mf.frq$MAHOM_F_P[mf.frq$ALT_FREQ_F<0.5] <- mapply(bt, mf.frq$HOM_ALT_F[mf.frq$ALT_FREQ_F<0.5],mf.frq$CT_TOT_F[mf.frq$ALT_FREQ_F<0.5],mf.frq$MAF_F[mf.frq$ALT_FREQ_F<0.5]^2)
#females: If alternative allele is major, binomial test for deficit of reference allele homozygotes 
mf.frq$MAHOM_F_P[mf.frq$ALT_FREQ_F>=0.5] <- mapply(bt, mf.frq$CT_TOT_F[mf.frq$ALT_FREQ_F>=0.5]-(mf.frq$HOM_ALT_F[mf.frq$ALT_FREQ_F>=0.5]+mf.frq$HET_F[mf.frq$ALT_FREQ_F>=0.5]),mf.frq$CT_TOT_F[mf.frq$ALT_FREQ_F>=0.5],mf.frq$MAF_F[mf.frq$ALT_FREQ_F>=0.5]^2)
mf.frq$MAHOM_F_P <- p.adjust(mf.frq$MAHOM_F_P,"BH")

#males: If alternative allele is minor, binomial test for deficit of alternative allele homozygotes
mf.frq$MAHOM_M_P[mf.frq$ALT_FREQ_M<0.5] <- mapply(bt, mf.frq$HOM_ALT_M[mf.frq$ALT_FREQ_M<0.5],mf.frq$CT_TOT_M[mf.frq$ALT_FREQ_M<0.5],mf.frq$MAF_M[mf.frq$ALT_FREQ_M<0.5]^2)
#males: If alternative allele is major, binomial test for deficit of reference allele homozygotes 
mf.frq$MAHOM_M_P[mf.frq$ALT_FREQ_M>=0.5] <- mapply(bt, mf.frq$CT_TOT_M[mf.frq$ALT_FREQ_M>=0.5]-(mf.frq$HOM_ALT_M[mf.frq$ALT_FREQ_M>=0.5]+mf.frq$HET_M[mf.frq$ALT_FREQ_M>=0.5]),mf.frq$CT_TOT_M[mf.frq$ALT_FREQ_M>=0.5],mf.frq$MAF_M[mf.frq$ALT_FREQ_M>=0.5]^2)
mf.frq$MAHOM_M_P <- p.adjust(mf.frq$MAHOM_M_P,"BH")

#both sexes: If alternative allele is minor, binomial test for deficit of alternative allele homozygotes
mf.frq$MAHOM_P[mf.frq$ALT_FREQ<0.5] <- mapply(bt, mf.frq$HOM_ALT[mf.frq$ALT_FREQ<0.5],mf.frq$CT_TOT[mf.frq$ALT_FREQ<0.5],mf.frq$MAF[mf.frq$ALT_FREQ<0.5]^2)
#both sexes: If alternative allele is major, binomial test for deficit of reference allele homozygotes 
mf.frq$MAHOM_P[mf.frq$ALT_FREQ>=0.5] <- mapply(bt, mf.frq$CT_TOT[mf.frq$ALT_FREQ>=0.5]-(mf.frq$HOM_ALT[mf.frq$ALT_FREQ>=0.5]+mf.frq$HET[mf.frq$ALT_FREQ>=0.5]),mf.frq$CT_TOT[mf.frq$ALT_FREQ>=0.5],mf.frq$MAF[mf.frq$ALT_FREQ>=0.5]^2)
mf.frq$MAHOM_P <- p.adjust(mf.frq$MAHOM_P,"BH")

#3. Filter based on significant difference in missing rate between the sexes (as per Kasimatis et al. 2020) 
mf.frq$MISSING_DIFF_P <- p.adjust(apply(mf.frq[,c("MISSING_M","CT_TOT_M","MISSING_F","CT_TOT_F")],1,function(x) chisq.test(matrix(x,nrow = 2))$p.value),"BH")
mf.frq <- mf.frq[,c(1:23,37:43)]

## Write out artefact-filtered files
write.table(subset(mf.frq,FIS_P>=0.05 & FIS_M_P>=0.05 & FIS_F_P>=0.05 & MAHOM_P>=0.05 & MAHOM_F_P>=0.05 & MAHOM_M_P>=0.05 & MISSING_DIFF_P>=0.05),paste0(dir,"ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.Fst"),row.names=F,quote=F,sep="\t")
write.table(subset(mf.frq,FIS_P>=0.05 & FIS_M_P>=0.05 & FIS_F_P>=0.05 & MAHOM_P>=0.05 & MAHOM_F_P>=0.05 & MAHOM_M_P>=0.05 & MISSING_DIFF_P>=0.05)[,c(1:25)],paste0(dir,"ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_small_v3.Fst"),row.names=F,quote=F,sep="\t")
write.table(subset(mf.frq,FIS_P>=0.05 & FIS_M_P>=0.05 & FIS_F_P>=0.05 & MAHOM_P>=0.05 & MAHOM_F_P>=0.05 & MAHOM_M_P>=0.05 & MISSING_DIFF_P>=0.05)[,c(2)],paste0(dir,"ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids"),row.names=F,quote=F,sep="\t")



## Imputed data, autosomes, MAF>0.01
mf.frq <- fread(paste0(dir,"ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01.Fst")) 

#Necessary columns for subsequent functions
mf.frq$ALT_FREQ <- 1-with(mf.frq,((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F))
mf.frq$ALT_FREQ_F <- 1-mf.frq$p_F_VIABILITY
mf.frq$ALT_FREQ_M <- 1-mf.frq$p_M_VIABILITY
mf.frq$HET <- mf.frq$HET_M+mf.frq$HET_F
mf.frq$HOM_ALT <- mf.frq$HOM_ALT_M+mf.frq$HOM_ALT_F
mf.frq$CT_TOT <- mf.frq$CT_TOT_M+mf.frq$CT_TOT_F
mf.frq$MAF_F <- with(mf.frq,ifelse(((p_F_VIABILITY*CT_TOT_F))/(CT_TOT_F)>0.5,1-((p_F_VIABILITY*CT_TOT_F))/(CT_TOT_F),((p_F_VIABILITY*CT_TOT_F))/(CT_TOT_F)))
mf.frq$MAF_M <- with(mf.frq,ifelse(((p_M_VIABILITY*CT_TOT_M))/(CT_TOT_M)>0.5,1-((p_M_VIABILITY*CT_TOT_M))/(CT_TOT_M),((p_M_VIABILITY*CT_TOT_M))/(CT_TOT_M)))
mf.frq$MISSING_M <- 115531-mf.frq$CT_TOT_M
mf.frq$MISSING_F <- 133490-mf.frq$CT_TOT_F 

#1. Filter for excess heterozygotes (elevated Fis), assuming s=0.2 and at equilibrium
mf.frq$FIS <- fis_func(mf.frq$ALT_FREQ_M,mf.frq$ALT_FREQ_F,(mf.frq$HET)/(mf.frq$CT_TOT))
mf.frq$FIS_M <- fis_func_one_sex_only(mf.frq$ALT_FREQ_M,(mf.frq$HET_M)/(mf.frq$CT_TOT_M))
mf.frq$FIS_F <- fis_func_one_sex_only(mf.frq$ALT_FREQ_F,(mf.frq$HET_F)/(mf.frq$CT_TOT_F))

#Updated check for excess heterozygosity
s <- 0.2
mf.frq$FIS_P <- p.adjust(pnorm(mf.frq$FIS,mean=(1/(2*mf.frq$CT_TOT)) + (mf.frq$MAF*(1-mf.frq$MAF)/4)*(s/(1 - mf.frq$MAF*s))^2,sd=sqrt(1/(mf.frq$CT_TOT)),lower.tail = F),"BH")
mf.frq$FIS_M_P <- p.adjust(pnorm(mf.frq$FIS_M,mean=(1/(2*mf.frq$CT_TOT_M)) + (mf.frq$MAF_M*(1-mf.frq$MAF_M)/4*(s/(1 - mf.frq$MAF_M*s))^2),sd=sqrt(1/(mf.frq$CT_TOT_M)),lower.tail = F),"BH")
mf.frq$FIS_F_P <- p.adjust(pnorm(mf.frq$FIS_F,mean=(1/(2*mf.frq$CT_TOT_F)) + (mf.frq$MAF_F*(1-mf.frq$MAF_F)/4*(s/(1 - mf.frq$MAF_F*s))^2),sd=sqrt(1/(mf.frq$CT_TOT_F)),lower.tail = F),"BH")

#2. Filter based on deficit of minor allele homozygotes (as per Kasimatis et al. 2020)
#bt function is a one-sided binomial test test that asks if the number of homozygotes is less than expected

#females: If alternative allele is minor, binomial test for deficit of alternative allele homozygotes
mf.frq$MAHOM_F_P[mf.frq$ALT_FREQ_F<0.5] <- mapply(bt, mf.frq$HOM_ALT_F[mf.frq$ALT_FREQ_F<0.5],mf.frq$CT_TOT_F[mf.frq$ALT_FREQ_F<0.5],mf.frq$MAF_F[mf.frq$ALT_FREQ_F<0.5]^2)
#females: If alternative allele is major, binomial test for deficit of reference allele homozygotes 
mf.frq$MAHOM_F_P[mf.frq$ALT_FREQ_F>=0.5] <- mapply(bt, mf.frq$CT_TOT_F[mf.frq$ALT_FREQ_F>=0.5]-(mf.frq$HOM_ALT_F[mf.frq$ALT_FREQ_F>=0.5]+mf.frq$HET_F[mf.frq$ALT_FREQ_F>=0.5]),mf.frq$CT_TOT_F[mf.frq$ALT_FREQ_F>=0.5],mf.frq$MAF_F[mf.frq$ALT_FREQ_F>=0.5]^2)
mf.frq$MAHOM_F_P <- p.adjust(mf.frq$MAHOM_F_P,"BH")

#males: If alternative allele is minor, binomial test for deficit of alternative allele homozygotes
mf.frq$MAHOM_M_P[mf.frq$ALT_FREQ_M<0.5] <- mapply(bt, mf.frq$HOM_ALT_M[mf.frq$ALT_FREQ_M<0.5],mf.frq$CT_TOT_M[mf.frq$ALT_FREQ_M<0.5],mf.frq$MAF_M[mf.frq$ALT_FREQ_M<0.5]^2)
#males: If alternative allele is major, binomial test for deficit of reference allele homozygotes 
mf.frq$MAHOM_M_P[mf.frq$ALT_FREQ_M>=0.5] <- mapply(bt, mf.frq$CT_TOT_M[mf.frq$ALT_FREQ_M>=0.5]-(mf.frq$HOM_ALT_M[mf.frq$ALT_FREQ_M>=0.5]+mf.frq$HET_M[mf.frq$ALT_FREQ_M>=0.5]),mf.frq$CT_TOT_M[mf.frq$ALT_FREQ_M>=0.5],mf.frq$MAF_M[mf.frq$ALT_FREQ_M>=0.5]^2)
mf.frq$MAHOM_M_P <- p.adjust(mf.frq$MAHOM_M_P,"BH")

#both sexes: If alternative allele is minor, binomial test for deficit of alternative allele homozygotes
mf.frq$MAHOM_P[mf.frq$ALT_FREQ<0.5] <- mapply(bt, mf.frq$HOM_ALT[mf.frq$ALT_FREQ<0.5],mf.frq$CT_TOT[mf.frq$ALT_FREQ<0.5],mf.frq$MAF[mf.frq$ALT_FREQ<0.5]^2)
#both sexes: If alternative allele is major, binomial test for deficit of reference allele homozygotes 
mf.frq$MAHOM_P[mf.frq$ALT_FREQ>=0.5] <- mapply(bt, mf.frq$CT_TOT[mf.frq$ALT_FREQ>=0.5]-(mf.frq$HOM_ALT[mf.frq$ALT_FREQ>=0.5]+mf.frq$HET[mf.frq$ALT_FREQ>=0.5]),mf.frq$CT_TOT[mf.frq$ALT_FREQ>=0.5],mf.frq$MAF[mf.frq$ALT_FREQ>=0.5]^2)
mf.frq$MAHOM_P <- p.adjust(mf.frq$MAHOM_P,"BH")

#3. Filter based on significant difference in missing rate between the sexes (as per Kasimatis et al. 2020) 
mf.frq$MISSING_DIFF_P <- p.adjust(apply(mf.frq[,c("MISSING_M","CT_TOT_M","MISSING_F","CT_TOT_F")],1,function(x) chisq.test(matrix(x,nrow = 2))$p.value),"BH")
mf.frq <- mf.frq[,c(1:23,37:43)]

## Write it out
write.table(subset(mf.frq,FIS_P>=0.05 & FIS_M_P>=0.05 & FIS_F_P>=0.05 & MAHOM_P>=0.05 & MAHOM_F_P>=0.05 & MAHOM_M_P>=0.05 & MISSING_DIFF_P>=0.05),paste0(dir,"ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst"),row.names=F,quote=F,sep="\t")
write.table(subset(mf.frq,FIS_P>=0.05 & FIS_M_P>=0.05 & FIS_F_P>=0.05 & MAHOM_P>=0.05 & MAHOM_F_P>=0.05 & MAHOM_M_P>=0.05 & MISSING_DIFF_P>=0.05)[,c(1:25)],paste0(dir,"ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_small_v3.Fst"),row.names=F,quote=F,sep="\t")
write.table(subset(auto.frq,FIS_P>=0.05 & FIS_M_P>=0.05 & FIS_F_P>=0.05 & MAHOM_P>=0.05 & MAHOM_F_P>=0.05 & MAHOM_M_P>=0.05 & MISSING_DIFF_P>=0.05)[,c(2)],paste0(dir,"ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids"),row.names=F,quote=F,sep="\t")

## Pm, pf, pm', pf', permuted1 (nchildren) ####


## Genotype data, autosomes, permuted1 (nchildren)
for (j in 1:22){
  print(paste(j,"initiated"))
  ## Frequencies, females
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm1.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c(1:4)]
  #Total genotype counts for adult females (with 0 #offspring)
  f.frq$CT_TOT <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for females (with 0 #offspring)
  f.frq$HET_F <- f0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for females (with 0 #offspring)
  f.frq$HOM_ALT_F <- f0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for females (with 0 #offspring)
  f.frq$M_HOM_REF <- f0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for females (with 0 #offspring)
  f.frq$M_HET <- f0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for females (with 0 #offspring)
  f.frq$M_HOM_ALT <- f0.frq$TWO_ALT_GENO_CTS*0
  rm(f0.frq)
  
  for (i in c(1:14,19)){
    f1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm1.Female",i,".gcount"))
    #Total genotype counts for adult females (with i #offspring), added to adult females with 0 #offspring
    f.frq$CT_TOT <- f.frq$CT_TOT+(f1.frq$HOM_REF_CT+f1.frq$HET_REF_ALT_CTS+f1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HET_F <- f.frq$HET_F+(f1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HOM_ALT_F <- f.frq$HOM_ALT_F+(f1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_REF <- f.frq$M_HOM_REF+(f1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HET <- f.frq$M_HET+(f1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_ALT <- f.frq$M_HOM_ALT+(f1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  ## Frequencies, males
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm1.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c(1:4)]
  #Total genotype counts for adult males (with 0 #offspring)
  m.frq$CT_TOT <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for males (with 0 #offspring)
  m.frq$HET_M <- m0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for males (with 0 #offspring)
  m.frq$HOM_ALT_M <- m0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for males (with 0 #offspring)
  m.frq$M_HOM_REF <- m0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for males (with 0 #offspring)
  m.frq$M_HET <- m0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for males (with 0 #offspring)
  m.frq$M_HOM_ALT <- m0.frq$TWO_ALT_GENO_CTS*0
  rm(m0.frq)
  
  for (i in c(1:15,17)){
    m1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm1.Male",i,".gcount"))
    #Total genotype counts for adult males (with i #offspring), added to adult males with 0 #offspring
    m.frq$CT_TOT <- m.frq$CT_TOT+(m1.frq$HOM_REF_CT+m1.frq$HET_REF_ALT_CTS+m1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HET_M <- m.frq$HET_M+(m1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HOM_ALT_M <- m.frq$HOM_ALT_M+(m1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_REF <- m.frq$M_HOM_REF+(m1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HET <- m.frq$M_HET+(m1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_ALT <- m.frq$M_HOM_ALT+(m1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  #Combine total allele frequency counts and 
  mf.frq <- cbind(f.frq,m.frq[,c(5:10)])
  names(mf.frq)[c(5:16)] <-c("CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M")
  #Reference allele frequency in adult females (Pf)
  mf.frq$p_F_VIABILITY <- 1-(mf.frq$HOM_ALT_F+(0.5*mf.frq$HET_F))/(mf.frq$CT_TOT_F)
  #Reference allele frequency in adult males (Pm)
  mf.frq$p_M_VIABILITY <- 1-(mf.frq$HOM_ALT_M+(0.5*mf.frq$HET_M))/(mf.frq$CT_TOT_M)
  #Adult Fst
  mf.frq$ADULT_FST <- nei_fst_func(mf.frq$p_F_VIABILITY,mf.frq$p_M_VIABILITY)
  #Projected reference allele frequency in adult females (Pf')
  mf.frq$p_F <- (mf.frq$M_HOM_REF_F+(0.5*mf.frq$M_HET_F))/(mf.frq$M_HOM_REF_F+mf.frq$M_HET_F+mf.frq$M_HOM_ALT_F) 
  #Projected reference allele frequency in adult females (Pm')
  mf.frq$p_M <- (mf.frq$M_HOM_REF_M+(0.5*mf.frq$M_HET_M))/(mf.frq$M_HOM_REF_M+mf.frq$M_HET_M+mf.frq$M_HOM_ALT_M)
  #Gametic Fst
  mf.frq$GAMETIC_FST <- nei_fst_func(mf.frq$p_F,mf.frq$p_M)
  
  write.table(mf.frq,paste0(dir,"ukb_fst_v3/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm1_no_further_maf_filter.Fst"),quote=F,row.names=F,sep="\t",col.names=F)
  
  print(paste(j,"completed"))
}


## Imputed data, autosomes
for (j in 1:22){
  print(paste(j,"initiated"))
  ## Genotype frequencies, females
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm1.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c(1:4)]
  #Total genotype counts for adult females (with 0 #offspring)
  f.frq$CT_TOT <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for females (with 0 #offspring)
  f.frq$HET_F <- f0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for females (with 0 #offspring)
  f.frq$HOM_ALT_F <- f0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for females (with 0 #offspring)
  f.frq$M_HOM_REF <- f0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for females (with 0 #offspring)
  f.frq$M_HET <- f0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for females (with 0 #offspring)
  f.frq$M_HOM_ALT <- f0.frq$TWO_ALT_GENO_CTS*0
  rm(f0.frq)
  
  for (i in c(1:14,19)){
    f1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm1.Female",i,".gcount"))
    #Total genotype counts for adult females (with i #offspring), added to adult females with 0 #offspring
    f.frq$CT_TOT <- f.frq$CT_TOT+(f1.frq$HOM_REF_CT+f1.frq$HET_REF_ALT_CTS+f1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HET_F <- f.frq$HET_F+(f1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HOM_ALT_F <- f.frq$HOM_ALT_F+(f1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_REF <- f.frq$M_HOM_REF+(f1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HET <- f.frq$M_HET+(f1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_ALT <- f.frq$M_HOM_ALT+(f1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  ## Genotype frequencies, males
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm1.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c(1:4)]
  #Total genotype counts for adult males (with 0 #offspring)
  m.frq$CT_TOT <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for males (with 0 #offspring)
  m.frq$HET_M <- m0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for males (with 0 #offspring)
  m.frq$HOM_ALT_M <- m0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for males (with 0 #offspring)
  m.frq$M_HOM_REF <- m0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for males (with 0 #offspring)
  m.frq$M_HET <- m0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for males (with 0 #offspring)
  m.frq$M_HOM_ALT <- m0.frq$TWO_ALT_GENO_CTS*0
  rm(m0.frq)
  
  for (i in c(1:15,17)){
    m1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm1.Male",i,".gcount"))
    #Total genotype counts for adult males (with i #offspring), added to adult males with 0 #offspring
    m.frq$CT_TOT <- m.frq$CT_TOT+(m1.frq$HOM_REF_CT+m1.frq$HET_REF_ALT_CTS+m1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HET_M <- m.frq$HET_M+(m1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HOM_ALT_M <- m.frq$HOM_ALT_M+(m1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_REF <- m.frq$M_HOM_REF+(m1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HET <- m.frq$M_HET+(m1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_ALT <- m.frq$M_HOM_ALT+(m1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  #Combine total allele frequency counts and 
  mf.frq <- cbind(f.frq,m.frq[,c(5:10)])
  names(mf.frq)[c(5:16)] <-c("CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M")
  #Reference allele frequency in adult females (Pf)
  mf.frq$p_F_VIABILITY <- 1-(mf.frq$HOM_ALT_F+(0.5*mf.frq$HET_F))/(mf.frq$CT_TOT_F)
  #Reference allele frequency in adult males (Pm)
  mf.frq$p_M_VIABILITY <- 1-(mf.frq$HOM_ALT_M+(0.5*mf.frq$HET_M))/(mf.frq$CT_TOT_M)
  #Adult Fst
  mf.frq$ADULT_FST <- nei_fst_func(mf.frq$p_F_VIABILITY,mf.frq$p_M_VIABILITY)
  #Projected reference allele frequency in adult females (Pf')
  mf.frq$p_F <- (mf.frq$M_HOM_REF_F+(0.5*mf.frq$M_HET_F))/(mf.frq$M_HOM_REF_F+mf.frq$M_HET_F+mf.frq$M_HOM_ALT_F) 
  #Projected reference allele frequency in adult females (Pm')
  mf.frq$p_M <- (mf.frq$M_HOM_REF_M+(0.5*mf.frq$M_HET_M))/(mf.frq$M_HOM_REF_M+mf.frq$M_HET_M+mf.frq$M_HOM_ALT_M)
  #Gametic Fst
  mf.frq$GAMETIC_FST <- nei_fst_func(mf.frq$p_F,mf.frq$p_M)
  
  
  write.table(mf.frq,paste0(dir,"ukb_fst_v3/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm1_no_info_filter_no_further_maf_filter.Fst"),quote=F,row.names=F,sep="\t",col.names=F)
  
  print(paste(j,"completed"))
}


## Pm, pf, pm', pf', permuted sexes (Perm2) ####


## Genotype data, autosomes
for (j in 1:22){
  print(paste(j,"initiated"))
  ## Genotype frequencies, females
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm2.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c(1:4)]
  #Total genotype counts for adult females (with 0 #offspring)
  f.frq$CT_TOT <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for females (with 0 #offspring)
  f.frq$HET_F <- f0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for females (with 0 #offspring)
  f.frq$HOM_ALT_F <- f0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for females (with 0 #offspring)
  f.frq$M_HOM_REF <- f0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for females (with 0 #offspring)
  f.frq$M_HET <- f0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for females (with 0 #offspring)
  f.frq$M_HOM_ALT <- f0.frq$TWO_ALT_GENO_CTS*0
  rm(f0.frq)
  
  for (i in c(1:14)){
    f1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm2.Female",i,".gcount"))
    #Total genotype counts for adult females (with i #offspring), added to adult females with 0 #offspring
    f.frq$CT_TOT <- f.frq$CT_TOT+(f1.frq$HOM_REF_CT+f1.frq$HET_REF_ALT_CTS+f1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HET_F <- f.frq$HET_F+(f1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HOM_ALT_F <- f.frq$HOM_ALT_F+(f1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_REF <- f.frq$M_HOM_REF+(f1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HET <- f.frq$M_HET+(f1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_ALT <- f.frq$M_HOM_ALT+(f1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  ## Genotype frequencies, males
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm2.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c(1:4)]
  #Total genotype counts for adult males (with 0 #offspring)
  m.frq$CT_TOT <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for males (with 0 #offspring)
  m.frq$HET_M <- m0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for males (with 0 #offspring)
  m.frq$HOM_ALT_M <- m0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for males (with 0 #offspring)
  m.frq$M_HOM_REF <- m0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for males (with 0 #offspring)
  m.frq$M_HET <- m0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for males (with 0 #offspring)
  m.frq$M_HOM_ALT <- m0.frq$TWO_ALT_GENO_CTS*0
  rm(m0.frq)
  
  for (i in c(1:15,17,19)){
    m1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm2.Male",i,".gcount"))
    #Total genotype counts for adult males (with i #offspring), added to adult males with 0 #offspring
    m.frq$CT_TOT <- m.frq$CT_TOT+(m1.frq$HOM_REF_CT+m1.frq$HET_REF_ALT_CTS+m1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HET_M <- m.frq$HET_M+(m1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HOM_ALT_M <- m.frq$HOM_ALT_M+(m1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_REF <- m.frq$M_HOM_REF+(m1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HET <- m.frq$M_HET+(m1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_ALT <- m.frq$M_HOM_ALT+(m1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  #Combine total allele frequency counts and 
  mf.frq <- cbind(f.frq,m.frq[,c(5:10)])
  names(mf.frq)[c(5:16)] <-c("CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M")
  #Reference allele frequency in adult females (Pf)
  mf.frq$p_F_VIABILITY <- 1-(mf.frq$HOM_ALT_F+(0.5*mf.frq$HET_F))/(mf.frq$CT_TOT_F)
  #Reference allele frequency in adult males (Pm)
  mf.frq$p_M_VIABILITY <- 1-(mf.frq$HOM_ALT_M+(0.5*mf.frq$HET_M))/(mf.frq$CT_TOT_M)
  #Adult Fst
  mf.frq$ADULT_FST <- nei_fst_func(mf.frq$p_F_VIABILITY,mf.frq$p_M_VIABILITY)
  #Projected reference allele frequency in adult females (Pf')
  mf.frq$p_F <- (mf.frq$M_HOM_REF_F+(0.5*mf.frq$M_HET_F))/(mf.frq$M_HOM_REF_F+mf.frq$M_HET_F+mf.frq$M_HOM_ALT_F) 
  #Projected reference allele frequency in adult females (Pm')
  mf.frq$p_M <- (mf.frq$M_HOM_REF_M+(0.5*mf.frq$M_HET_M))/(mf.frq$M_HOM_REF_M+mf.frq$M_HET_M+mf.frq$M_HOM_ALT_M)
  #Gametic Fst
  mf.frq$GAMETIC_FST <- nei_fst_func(mf.frq$p_F,mf.frq$p_M)
  write.table(mf.frq,paste0(dir,"ukb_fst_v3/ukb52049_chr",j,"_sampleqc1_snpqc1b_perm2_no_further_maf_filter.Fst"),quote=F,row.names=F,sep="\t",col.names=F)
  
  print(paste(j,"completed"))
}

## Imputed data, autosomes
for (j in 1:22){
  print(paste(j,"initiated"))
  ## Genotype frequencies, females
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm2.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c(1:4)]
  #Total genotype counts for adult females (with 0 #offspring)
  f.frq$CT_TOT <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for females (with 0 #offspring)
  f.frq$HET_F <- f0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for females (with 0 #offspring)
  f.frq$HOM_ALT_F <- f0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for females (with 0 #offspring)
  f.frq$M_HOM_REF <- f0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for females (with 0 #offspring)
  f.frq$M_HET <- f0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for females (with 0 #offspring)
  f.frq$M_HOM_ALT <- f0.frq$TWO_ALT_GENO_CTS*0
  rm(f0.frq)
  
  for (i in c(1:14)){
    f1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm2.Female",i,".gcount"))
    #Total genotype counts for adult females (with i #offspring), added to adult females with 0 #offspring
    f.frq$CT_TOT <- f.frq$CT_TOT+(f1.frq$HOM_REF_CT+f1.frq$HET_REF_ALT_CTS+f1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HET_F <- f.frq$HET_F+(f1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$HOM_ALT_F <- f.frq$HOM_ALT_F+(f1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_REF <- f.frq$M_HOM_REF+(f1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HET <- f.frq$M_HET+(f1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate for females (with i #offspring), added to adult females with 0 #offspring
    f.frq$M_HOM_ALT <- f.frq$M_HOM_ALT+(f1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  ## Genotype frequencies, males
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm2.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c(1:4)]
  #Total genotype counts for adult males (with 0 #offspring)
  m.frq$CT_TOT <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  #Heterozygote counts for males (with 0 #offspring)
  m.frq$HET_M <- m0.frq$HET_REF_ALT_CTS
  #Homozygote alternate counts for males (with 0 #offspring)
  m.frq$HOM_ALT_M <- m0.frq$TWO_ALT_GENO_CTS
  #Projected homozygote reference counts for males (with 0 #offspring)
  m.frq$M_HOM_REF <- m0.frq$HOM_REF_CT*0
  #Projected heterozygote counts for males (with 0 #offspring)
  m.frq$M_HET <- m0.frq$HET_REF_ALT_CTS*0
  #Projected homozygote alternate counts for males (with 0 #offspring)
  m.frq$M_HOM_ALT <- m0.frq$TWO_ALT_GENO_CTS*0
  rm(m0.frq)
  
  for (i in c(1:15,17,19)){
    m1.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm2.Male",i,".gcount"))
    #Total genotype counts for adult males (with i #offspring), added to adult males with 0 #offspring
    m.frq$CT_TOT <- m.frq$CT_TOT+(m1.frq$HOM_REF_CT+m1.frq$HET_REF_ALT_CTS+m1.frq$TWO_ALT_GENO_CTS)
    #Heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HET_M <- m.frq$HET_M+(m1.frq$HET_REF_ALT_CTS)
    #Homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$HOM_ALT_M <- m.frq$HOM_ALT_M+(m1.frq$TWO_ALT_GENO_CTS)
    #Projected homozygote reference counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_REF <- m.frq$M_HOM_REF+(m1.frq$HOM_REF_CT*i)
    #Projected heterozygote counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HET <- m.frq$M_HET+(m1.frq$HET_REF_ALT_CTS*i)
    #Projected homozygote alternate counts for males (with i #offspring), added to adult males with 0 #offspring
    m.frq$M_HOM_ALT <- m.frq$M_HOM_ALT+(m1.frq$TWO_ALT_GENO_CTS*i)
  }
  
  #Combine total allele frequency counts and 
  mf.frq <- cbind(f.frq,m.frq[,c(5:10)])
  names(mf.frq)[c(5:16)] <-c("CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M")
  #Reference allele frequency in adult females (Pf)
  mf.frq$p_F_VIABILITY <- 1-(mf.frq$HOM_ALT_F+(0.5*mf.frq$HET_F))/(mf.frq$CT_TOT_F)
  #Reference allele frequency in adult males (Pm)
  mf.frq$p_M_VIABILITY <- 1-(mf.frq$HOM_ALT_M+(0.5*mf.frq$HET_M))/(mf.frq$CT_TOT_M)
  #Adult Fst
  mf.frq$ADULT_FST <- nei_fst_func(mf.frq$p_F_VIABILITY,mf.frq$p_M_VIABILITY)
  #Projected reference allele frequency in adult females (Pf')
  mf.frq$p_F <- (mf.frq$M_HOM_REF_F+(0.5*mf.frq$M_HET_F))/(mf.frq$M_HOM_REF_F+mf.frq$M_HET_F+mf.frq$M_HOM_ALT_F) 
  #Projected reference allele frequency in adult females (Pm')
  mf.frq$p_M <- (mf.frq$M_HOM_REF_M+(0.5*mf.frq$M_HET_M))/(mf.frq$M_HOM_REF_M+mf.frq$M_HET_M+mf.frq$M_HOM_ALT_M)
  #Gametic Fst
  mf.frq$GAMETIC_FST <- nei_fst_func(mf.frq$p_F,mf.frq$p_M)
  
  write.table(mf.frq,paste0(dir,"ukb_fst_v3/ukb_imp_chr",j,"_v3_sampleqc1_snpqc1b_perm2_no_info_filter_no_further_maf_filter.Fst"),quote=F,row.names=F,sep="\t",col.names=F)
  
  print(paste(j,"completed"))
}


#Mean/var LRS (perm2; permuted sexes)####


##Genotyped 
#Female LRS
f.frq.all.chrom.perm <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c("CHROM","ID","REF","ALT")]
  f.frq$NOT_MISSING_CT_0 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female1.gcount"))
  f.frq$NOT_MISSING_CT_1 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female2.gcount"))
  f.frq$NOT_MISSING_CT_2 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female3.gcount"))
  f.frq$NOT_MISSING_CT_3 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female4.gcount"))
  f.frq$NOT_MISSING_CT_4 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female5.gcount"))
  f.frq$NOT_MISSING_CT_5 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female6.gcount"))
  f.frq$NOT_MISSING_CT_6 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female7.gcount"))
  f.frq$NOT_MISSING_CT_7 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female8.gcount"))
  f.frq$NOT_MISSING_CT_8 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female9.gcount"))
  f.frq$NOT_MISSING_CT_9 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female10.gcount"))
  f.frq$NOT_MISSING_CT_10 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female11.gcount"))
  f.frq$NOT_MISSING_CT_11 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female12.gcount"))
  f.frq$NOT_MISSING_CT_12 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female13.gcount"))
  f.frq$NOT_MISSING_CT_13 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Female14.gcount"))
  f.frq$NOT_MISSING_CT_14 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f.frq.all.chrom.perm[[i]] <- f.frq
  print(i)
}
f.frq.all.chrom.perm <- do.call(rbind,f.frq.all.chrom.perm)

#Males
m.frq.all.chrom.perm <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c("CHROM","ID","REF","ALT")]
  m.frq$NOT_MISSING_CT_0 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male1.gcount"))
  m.frq$NOT_MISSING_CT_1 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male2.gcount"))
  m.frq$NOT_MISSING_CT_2 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male3.gcount"))
  m.frq$NOT_MISSING_CT_3 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male4.gcount"))
  m.frq$NOT_MISSING_CT_4 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male5.gcount"))
  m.frq$NOT_MISSING_CT_5 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male6.gcount"))
  m.frq$NOT_MISSING_CT_6 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male7.gcount"))
  m.frq$NOT_MISSING_CT_7 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male8.gcount"))
  m.frq$NOT_MISSING_CT_8 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male9.gcount"))
  m.frq$NOT_MISSING_CT_9 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male10.gcount"))
  m.frq$NOT_MISSING_CT_10 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male11.gcount"))
  m.frq$NOT_MISSING_CT_11 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male12.gcount"))
  m.frq$NOT_MISSING_CT_12 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male13.gcount"))
  m.frq$NOT_MISSING_CT_13 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male14.gcount"))
  m.frq$NOT_MISSING_CT_14 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male15.gcount"))
  m.frq$NOT_MISSING_CT_15 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male17.gcount"))
  m.frq$NOT_MISSING_CT_17 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr",i,"_sampleqc1_snpqc1b_perm2.Male19.gcount"))
  m.frq$NOT_MISSING_CT_19 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m.frq.all.chrom.perm[[i]] <- m.frq
  print(i)
}
m.frq.all.chrom.perm <- do.call(rbind,m.frq.all.chrom.perm)

#Variance and mean^2 fitness for each locus, females
#Column names = LRS
names(f.frq.all.chrom.perm)[5:19] <- as.character(c(0:14))
f.frq.all.chrom2.perm <- f.frq.all.chrom.perm[,5:19]
#Histogram bins
bins <- as.numeric(names(f.frq.all.chrom2.perm))
#Total number of LRS values across bins (n)
n <- rowSums(f.frq.all.chrom2.perm)
#Mean = sum(bin*LRS-per-bin)/n
meanf.geno.perm <- vector()
meanf.geno.perm <- foreach (i=1:nrow(f.frq.all.chrom2.perm), .combine='c') %dopar% {
  (sum(bins*f.frq.all.chrom2.perm[i])/n[i])
}
#Var = sum((bin-mean)^2*LRS-per-bin)/n
varf.geno.perm <- vector()
varf.geno.perm <- foreach (i=1:nrow(f.frq.all.chrom2.perm), .combine='c') %dopar% {
  sum((bins-meanf.geno.perm[i])^2*f.frq.all.chrom2.perm[i])/n[i]
}
#Mean^2 
meansqf.geno.perm <- meanf.geno.perm^2

#Variance and mean^2 fitness for each locus, males
#Column names = LRS
names(m.frq.all.chrom.perm)[5:22] <- as.character(c(0:15,17,19))
m.frq.all.chrom2.perm <- m.frq.all.chrom.perm[,5:22]
#Histogram bins
bins <- as.numeric(names(m.frq.all.chrom2.perm))
#Total number of LRS values across bins (n)
n <- rowSums(m.frq.all.chrom2.perm)
#Mean = sum(bin*LRS-per-bin)/n
meanm.geno.perm <- vector()
for(i in 1:nrow(m.frq.all.chrom2.perm)){
  meanm.geno.perm[i] <- (sum(bins*m.frq.all.chrom2.perm[i])/n[i])
}
#Var = sum((bin-mean)^2*LRS-per-bin)/n
varm.geno.perm <- vector()
for(i in 1:nrow(m.frq.all.chrom2.perm)){
  varm.geno.perm[i] <- sum((bins-meanm.geno.perm[i])^2*m.frq.all.chrom2.perm[i,])/n[i]
}
#Mean^2 
meansqm.geno.perm <- meanm.geno.perm^2


mf.frq.all.chrom.perm <- cbind(f.frq.all.chrom.perm[,1:4],varf.geno.perm,meansqf.geno.perm,varm.geno.perm,meansqm.geno.perm)
names(mf.frq.all.chrom.perm)[5:8] <- c("VAR_LRS_PERM2_F","MEANSQ_LRS_PERM2_F","VAR_LRS_PERM2_M","MEANSQ_LRS_PERM2_M")


#write.table(mf.frq.all.chrom.perm,paste0(dir,"/ukb_locus_by_locus_fitness_var_and_mean_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_perm2.var_and_mean.txt"),quote=F,row.names=F,sep="\t")



##Imputed
registerDoParallel(cores = 4)
#Female LRS
f.frq.chrom.pos.imp <- list()
meansqf.imp <- list()
varf.imp <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female0.gcount"))
  names(f0.frq)[1] <- "CHROM"
  f.frq <- f0.frq[,c("CHROM","ID","REF","ALT")]
  f.frq$NOT_MISSING_CT_0 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female1.gcount"))
  f.frq$NOT_MISSING_CT_1 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female2.gcount"))
  f.frq$NOT_MISSING_CT_2 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female3.gcount"))
  f.frq$NOT_MISSING_CT_3 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female4.gcount"))
  f.frq$NOT_MISSING_CT_4 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female5.gcount"))
  f.frq$NOT_MISSING_CT_5 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female6.gcount"))
  f.frq$NOT_MISSING_CT_6 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female7.gcount"))
  f.frq$NOT_MISSING_CT_7 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female8.gcount"))
  f.frq$NOT_MISSING_CT_8 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female9.gcount"))
  f.frq$NOT_MISSING_CT_9 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female10.gcount"))
  f.frq$NOT_MISSING_CT_10 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female11.gcount"))
  f.frq$NOT_MISSING_CT_11 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female12.gcount"))
  f.frq$NOT_MISSING_CT_12 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female13.gcount"))
  f.frq$NOT_MISSING_CT_13 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  f0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Female14.gcount"))
  f.frq$NOT_MISSING_CT_14 <- f0.frq$HOM_REF_CT+f0.frq$HET_REF_ALT_CTS+f0.frq$TWO_ALT_GENO_CTS
  print(i)
  
  #Keep only MAF/artefact-filtered sites
  f.frq <- subset(f.frq,ID %in% good.ids$ID)
  
  #Variance and mean^2 fitness for each locus
  #Column names = LRS
  names(f.frq)[5:19] <- as.character(c(0:14))
  f.frq2 <- f.frq[,5:19]
  #Histogram bins
  bins <- as.numeric(names(f.frq2))
  #Total number of LRS values across bins (n)
  n <- rowSums(f.frq2)
  #Mean = sum(bin*LRS-per-bin)/n
  meansqf.imp[[i]] <- foreach (j=1:nrow(f.frq2), .combine='rbind') %dopar% {
    (sum(bins*f.frq2[j])/n[j])^2
  }
  #Var = sum((bin-mean)^2*LRS-per-bin)/n
  varf.imp[[i]] <- foreach (j=1:nrow(f.frq2), .combine='rbind') %dopar% {
    sum((bins-sqrt(meansqf.imp[[i]][j]))^2*f.frq2[j])/n[j]
  }
  #Positional information
  f.frq.chrom.pos.imp[[i]] <- f.frq[,1:4]
}

meansqf.perm2.imp <- unlist(meansqf.imp)
varf.perm2.imp <- unlist(varf.imp)
f.frq.chrom.pos.perm2.imp <- do.call(rbind,f.frq.chrom.pos.imp)

#Males
meansqm.imp <- list()
varm.imp <- list()
for (i in 1:22){
  #Number of individuals for each offspring count class
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male0.gcount"))
  names(m0.frq)[1] <- "CHROM"
  m.frq <- m0.frq[,c("CHROM","ID","REF","ALT")]
  m.frq$NOT_MISSING_CT_0 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male1.gcount"))
  m.frq$NOT_MISSING_CT_1 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male2.gcount"))
  m.frq$NOT_MISSING_CT_2 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male3.gcount"))
  m.frq$NOT_MISSING_CT_3 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male4.gcount"))
  m.frq$NOT_MISSING_CT_4 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male5.gcount"))
  m.frq$NOT_MISSING_CT_5 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male6.gcount"))
  m.frq$NOT_MISSING_CT_6 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male7.gcount"))
  m.frq$NOT_MISSING_CT_7 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male8.gcount"))
  m.frq$NOT_MISSING_CT_8 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male9.gcount"))
  m.frq$NOT_MISSING_CT_9 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male10.gcount"))
  m.frq$NOT_MISSING_CT_10 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male11.gcount"))
  m.frq$NOT_MISSING_CT_11 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male12.gcount"))
  m.frq$NOT_MISSING_CT_12 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male13.gcount"))
  m.frq$NOT_MISSING_CT_13 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male14.gcount"))
  m.frq$NOT_MISSING_CT_14 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male15.gcount"))
  m.frq$NOT_MISSING_CT_15 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male17.gcount"))
  m.frq$NOT_MISSING_CT_17 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  m0.frq <- fread(paste0(dir,"/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr",i,"_v3_sampleqc1_snpqc1b_perm2.Male19.gcount"))
  m.frq$NOT_MISSING_CT_19 <- m0.frq$HOM_REF_CT+m0.frq$HET_REF_ALT_CTS+m0.frq$TWO_ALT_GENO_CTS
  
  #Keep only MAF/artefact-filtered sites
  m.frq <- subset(m.frq,ID %in% good.ids$ID)
  
  #Variance and mean^2 fitness for each locus
  #Column names = LRS
  names(m.frq)[5:22] <- as.character(c(0:15,17,19))
  m.frq2 <- m.frq[,5:22]
  #Histogram bins
  bins <- as.numeric(names(m.frq2))
  #Total number of LRS values across bins (n)
  n <- rowSums(m.frq2)
  #Mean = sum(bin*LRS-per-bin)/n
  meansqm.imp[[i]] <- foreach (j=1:nrow(m.frq2), .combine='rbind') %dopar% {
    (sum(bins*m.frq2[j])/n[j])^2
  }
  #Var = sum((bin-mean)^2*LRS-per-bin)/n
  varm.imp[[i]] <- foreach (j=1:nrow(m.frq2), .combine='rbind') %dopar% {
    sum((bins-sqrt(meansqm.imp[[i]][j]))^2*m.frq2[j])/n[j]
  }
}
meansqm.perm2.imp <- unlist(meansqm.imp)
varm.perm2.imp <- unlist(varm.imp)


mf.frq.all.chrom.perm2 <- cbind(f.frq.chrom.pos.perm2.imp[,1:4],varf.perm2.imp,meansqf.perm2.imp,varm.perm2.imp,meansqm.perm2.imp)
names(mf.frq.all.chrom.perm2)[5:8] <- c("VAR_LRS_PERM2_F","MEANSQ_LRS_PERM2_F","VAR_LRS_PERM2_M","MEANSQ_LRS_PERM2_M")


write.table(mf.frq.all.chrom.perm2,paste0(dir,"/ukb_locus_by_locus_fitness_var_and_mean_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm2_info_above0.8_maf_above0.01_artefacts_removed_v3.var_and_mean.txt"),quote=F,row.names=F,sep="\t")


##Prep annotations####


##Genotyped data

#POS for all autosomal genotyped sites (no MAF or artefact filtering)
all.ids_and_pos <- fread(paste0(dir,"/ukb_annotations_v3/ukb52049_chrAUTO.bim"))
names(all.ids_and_pos) <- c("CHROM","ID","Whatever","POS","ALT","REF")

#CHROM/ID/REF/ALT data for MAF + artefact-filtered sites
good.ids <- fread(paste0(dir,"/ukb_annotations_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.4columns"))

#Add positional information to MAF + artefact-filtered IDs
good.ids_and_pos <- merge(good.ids,all.ids_and_pos[,c("CHROM","ID","POS","REF","ALT")],by=c("CHROM","ID","REF","ALT"),all.x=T)
#write.table(good.ids_and_pos[,c("CHROM","POS","ID","REF","ALT")],paste0(dir,"/ukb_annotations_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.vcf"),quote=F,row.names=F,sep="\t")


##Imputed data

#POS for all autosomal genotyped sites (no MAF or artefact filtering)
all.ids_and_pos <- fread(paste0(dir,"/ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b.bim"))
names(all.ids_and_pos) <- c("CHROM","ID","Whatever","POS","ALT","REF")

#CHROM/ID/REF/ALT data for MAF + artefact-filtered sites
good.ids <- fread(paste0(dir,"/ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.4columns"))

#Add positional information to MAF + artefact-filtered IDs
good.ids_and_pos <- merge(good.ids,all.ids_and_pos[,c("CHROM","ID","POS","REF","ALT")],by=c("CHROM","ID","REF","ALT"),all.x=T)
#write.table(good.ids_and_pos[,c("CHROM","POS","ID","REF","ALT")],paste0(dir,"ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.vcf"),quote=F,row.names=F,sep="\t")


## Kasimatis et al. 2020 hits?####
hits_all <- c("rs9870157","rs145369881","rs77638744","rs9508454","rs1048990","rs75745570","rs114928327","rs11598874","rs11032483","rs75212444","rs7298104","rs116890400","rs73196350")
hits_biobank <- c("rs75745570","rs114928327","rs11598874","rs11032483","rs75212444","rs7298104","rs116890400","rs73196350")

##Genotyped data
#Post artefact-filtering
good.ids <-  fread(paste0(dir,"/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.Fst"))
names(good.ids)[1:2] <- c("CHROM","ID")
head(subset(good.ids,ID %in% hits_biobank))
#All removed in genotyped data

#Post artefact-filtering
good.ids <-  fread(paste0(dir,"/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst"))
nrow(subset(good.ids,ID %in% hits_biobank))
#LD pruned
ldpruned <- fread(paste0(dir,"ukb_ldpruning_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_ldpruned_0.2.prune.in"),head=F)
good.ids$LD_PRUNED <- ifelse(good.ids$ID %in% ldpruned$V1,1,0)
good.ids.ld <- subset(good.ids,LD_PRUNED==1)
nrow(subset(good.ids.ld,ID %in% hits_biobank))
#All have been filtered out in LD-pruned imputed data

