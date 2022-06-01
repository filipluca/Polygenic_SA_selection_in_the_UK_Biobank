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


## Import Fst ####

#Imputed data
mf.frq <- fread(paste0(dir,"/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst"))
mf.frq.perm <- fread(paste0(dir,"/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm1_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst"))
mf.frq.perm2 <- fread(paste0(dir,"/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm2_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst"))
names(mf.frq) <- c("CHROM","ID","REF","ALT","CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M","p_F_VIABILITY","p_M_VIABILITY", "FST_VIABILITY","p_F","p_M","FST_GAMETIC","MAF","FIS_P","FIS_M_P","FIS_F_P","MAHOM_F_P","MAHOM_M_P","MAHOM_P","MISSING_DIFF_P") 
mf.frq <- mf.frq[,c("CHROM","ID","REF","ALT","CT_TOT_F","HET_F","HOM_ALT_F","M_HOM_REF_F","M_HET_F","M_HOM_ALT_F","CT_TOT_M","HET_M","HOM_ALT_M","M_HOM_REF_M","M_HET_M","M_HOM_ALT_M","p_F_VIABILITY","p_M_VIABILITY", "FST_VIABILITY","p_F","p_M","FST_GAMETIC","MAF")]
names(mf.frq.perm) <- c("CHROM","ID","REF","ALT","CT_TOT_F_PERM","HET_F_PERM","HOM_ALT_F_PERM","M_HOM_REF_F_PERM","M_HET_F_PERM","M_HOM_ALT_F_PERM","CT_TOT_M_PERM","HET_M_PERM","HOM_ALT_M_PERM","M_HOM_REF_M_PERM","M_HET_M_PERM","M_HOM_ALT_M_PERM","p_F_VIABILITY_PERM","p_M_VIABILITY_PERM", "FST_VIABILITY_PERM","p_F_PERM","p_M_PERM","FST_GAMETIC_PERM")
mf.frq.perm <- mf.frq.perm[,c("CHROM","ID","REF","ALT","p_F_VIABILITY_PERM","p_M_VIABILITY_PERM", "FST_VIABILITY_PERM","p_F_PERM","p_M_PERM","FST_GAMETIC_PERM")]
names(mf.frq.perm2) <- c("CHROM","ID","REF","ALT","CT_TOT_F_PERM2","HET_F_PERM2","HOM_ALT_F_PERM2","M_HOM_REF_F_PERM2","M_HET_F_PERM2","M_HOM_ALT_F_PERM2","CT_TOT_M_PERM2","HET_M_PERM2","HOM_ALT_M_PERM2","M_HOM_REF_M_PERM2","M_HET_M_PERM2","M_HOM_ALT_M_PERM2","p_F_VIABILITY_PERM2","p_M_VIABILITY_PERM2", "FST_VIABILITY_PERM2","p_F_PERM2","p_M_PERM2","FST_GAMETIC_PERM2")
mf.frq.perm2 <- mf.frq.perm2[,c("CHROM","ID","REF","ALT","CT_TOT_F_PERM2","CT_TOT_M_PERM2","p_F_VIABILITY_PERM2","p_M_VIABILITY_PERM2", "FST_VIABILITY_PERM2","p_F_PERM2","p_M_PERM2","FST_GAMETIC_PERM2")]
##Mean-sq and variance of the number of offspring, per-locus
mf.frq.meansq_and_var <- read.table(paste0(dir,"ukb_locus_by_locus_fitness_var_and_mean_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.var_and_mean.txt"),h=T)
mf.frq.perm2.meansq_and_var <- read.table(paste0(dir,"ukb_locus_by_locus_fitness_var_and_mean_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm2_info_above0.8_maf_above0.01_artefacts_removed_v3.var_and_mean.txt"),h=T)
#Merge Fst data into a single table
mf.frqs <- Reduce(function(...) merge(..., all.x=T,by=c("CHROM","ID","REF","ALT")),list(mf.frq,mf.frq.perm,mf.frq.perm2,mf.frq.meansq_and_var,mf.frq.perm2.meansq_and_var))
rm(mf.frq)
rm(mf.frq.perm)
rm(mf.frq.perm2)
rm(mf.frq.meansq_and_var)
rm(mf.frq.perm2.meansq_and_var)
#Annotations
effects <- fread(paste0(dir,"/ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_refcorrected.ann.effect.gene.txt"),h=T)
names(effects) <- c("CHROM","POS","ID","REF","ALT","EFFECT","GENE")
mf.frqs <- merge(mf.frqs,effects,by=c("ID","CHROM"),all.x = T)
rm(effects)
#Remove duplicate REF/ALT in auto.frqs; duplicates occur because REF/ALT order is not always the same between mf.frqs and effects
mf.frqs <- subset(mf.frqs,!duplicated(mf.frqs[,c("ID","CHROM","POS","REF.x","ALT.x")]))
#7,844,719
#make sure that the alleles are the same b/w mf.frqs and effects
mf.frqs <- subset(mf.frqs,(REF.x==REF.y | REF.x==ALT.y) & (ALT.x==REF.y | ALT.x==ALT.y))
#7,842,458
#For some reason, "Affx-8165091" does not have annotation or positional information. Add it manually
mf.frqs[mf.frqs$ID=="Affx-8165091",c("POS")] <- 53011955
mf.frqs <- mf.frqs[order(mf.frqs$CHROM,mf.frqs$ID),]
#Categorise effects into genic/non-genic using 'intergenic region' annotation
mf.frqs$GENIC_OR_NOT <- ifelse(mf.frqs$EFFECT %in% c("intergenic_region"),0,1)
#LD-pruning data
ldpruned <- fread(paste0(dir,"ukb_ldpruning_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_ldpruned_0.2.prune.in"),head=F)
mf.frqs$LD_PRUNED <- ifelse(mf.frqs$ID %in% ldpruned$V1,1,0)
mf.frqs <- subset(mf.frqs,LD_PRUNED==1)
rm(ldpruned)
#Some additional columns
mf.frqs$FIS_F <- fis_func_one_sex_only(mf.frqs$p_F_VIABILITY,(mf.frqs$HET_F)/(mf.frqs$CT_TOT_F))
mf.frqs$FIS_M <- fis_func_one_sex_only(mf.frqs$p_M_VIABILITY,(mf.frqs$HET_M)/(mf.frqs$CT_TOT_M))
mf.frqs$MAF_PERM2 <- with(mf.frqs,ifelse(((p_M_VIABILITY_PERM2*CT_TOT_M_PERM2)+(p_F_VIABILITY_PERM2*CT_TOT_F_PERM2))/(CT_TOT_M_PERM2+CT_TOT_F_PERM2)>0.5,1-((p_M_VIABILITY_PERM2*CT_TOT_M_PERM2)+(p_F_VIABILITY_PERM2*CT_TOT_F_PERM2))/(CT_TOT_M_PERM2+CT_TOT_F_PERM2),((p_M_VIABILITY_PERM2*CT_TOT_M_PERM2)+(p_F_VIABILITY_PERM2*CT_TOT_F_PERM2))/(CT_TOT_M_PERM2+CT_TOT_F_PERM2))) 


## Import Lst####

# Imputed
mf.gwas <- fread(paste0(dir,"/ukb_gwas_v3/using_imputed_modelSNPs/stats_gwas_mf_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab"))
names(mf.gwas) <- c("ID","CHROM","POS","GENPOS","REF","ALT","REF_FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")
mf.gwas.perm <- fread(paste0(dir,"/ukb_gwas_v3/using_imputed_modelSNPs/stats_gwas_mf_perm2_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab"))
names(mf.gwas.perm) <- c("ID","CHROM","POS","GENPOS","REF","ALT","REF_FREQ_PERM2","INFO","BETA_PERM2","SE_PERM2","P_BOLT_LMM_INF_PERM2")
mf.gwases <- merge(mf.gwas,mf.gwas.perm[,c("ID","CHROM","POS","REF_FREQ_PERM2","BETA_PERM2","SE_PERM2","P_BOLT_LMM_INF_PERM2")],by=c("ID","CHROM","POS"),all.x=T)
rm(mf.gwas.perm)
rm(mf.gwas)
#Annotations
effects <- fread(paste0(dir,"/ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_refcorrected.ann.effect.gene.txt"),h=T)
names(effects) <- c("CHROM","POS","ID","REF","ALT","EFFECT","GENE")
mf.gwases <- merge(mf.gwases,effects,by=c("ID","CHROM","POS"),all.x = T)
rm(effects)
#Remove duplicate REF/ALT in auto.frqs; duplicates occur because REF/ALT order is not always the same between mf.gwases and effects
mf.gwases <- subset(mf.gwases,!duplicated(mf.gwases[,c("ID","CHROM","POS","REF.x","ALT.x")]))
#make sure that the alleles are the same b/w mf.frqs and effects
mf.gwases <- subset(mf.gwases,(REF.x==REF.y | REF.x==ALT.y) & (ALT.x==REF.y | ALT.x==ALT.y))
#7,842,458
#For some reason, "Affx-8165091" does not have annotation or positional information. Add it manually
mf.gwases[mf.gwases$ID=="Affx-8165091",c("POS")] <- 53011955
mf.gwases <- mf.gwases[order(mf.gwases$CHROM,mf.gwases$ID),]
#Categorise effects into genic/non-genic using 'intergenic region' annotation
mf.gwases$GENIC_OR_NOT <- ifelse(mf.gwases$EFFECT %in% c("intergenic_region"),0,1)
#LD-pruning data
ldpruned <- fread(paste0(dir,"ukb_ldpruning_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_ldpruned_0.2.prune.in"),head=F)
mf.gwases$LD_PRUNED <- ifelse(mf.gwases$ID %in% ldpruned$V1,1,0)
mf.gwases <- subset(mf.gwases,LD_PRUNED==1)
rm(ldpruned)
#Some additional columns
mf.gwases$MAF <- ifelse(mf.gwases$REF_FREQ<0.5,mf.gwases$REF_FREQ,1-mf.gwases$REF_FREQ)
mf.gwases$MAF_PERM2 <- ifelse(mf.gwases$REF_FREQ_PERM2<0.5,mf.gwases$REF_FREQ_PERM2,1-mf.gwases$REF_FREQ_PERM2)

## Import t####

## Imputed
#Males, observed
lrs.gwas.m <- fread(paste0(dir,"ukb_gwas_v3/using_imputed_modelSNPs/stats_gwas_nchildren_male_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab"))
names(lrs.gwas.m) <- c("ID","CHROM","POS","GENPOS","REF","ALT","REF_FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")
names(lrs.gwas.m)[c(7,9:11)] <- paste(names(lrs.gwas.m)[c(7,9:11)],"M",sep = "_")
#Females, observed
lrs.gwas.f <- fread(paste0(dir,"ukb_gwas_v3/using_imputed_modelSNPs/stats_gwas_nchildren_female_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab"))
names(lrs.gwas.f) <- c("ID","CHROM","POS","GENPOS","REF","ALT","REF_FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")
names(lrs.gwas.f)[c(7,9:11)] <- paste(names(lrs.gwas.f)[c(7,9:11)],"F",sep = "_")
#Males, permuted
lrs.gwas.m.perm <- fread(paste0(dir,"ukb_gwas_v3/using_imputed_modelSNPs/stats_gwas_nchildren_male_perm_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab"))
names(lrs.gwas.m.perm) <- c("ID","CHROM","POS","GENPOS","REF","ALT","REF_FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")
names(lrs.gwas.m.perm)[c(7,9:11)] <- paste(names(lrs.gwas.m.perm)[c(7,9:11)],"M_PERM",sep = "_")
#Females, permuted
lrs.gwas.f.perm <- fread(paste0(dir,"ukb_gwas_v3/using_imputed_modelSNPs/stats_gwas_nchildren_female_perm_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab"))
names(lrs.gwas.f.perm) <- c("ID","CHROM","POS","GENPOS","REF","ALT","REF_FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")
names(lrs.gwas.f.perm)[c(7,9:11)] <- paste(names(lrs.gwas.f.perm)[c(7,9:11)],"F_PERM",sep = "_")
#Combined
mf.gwases.lrs <- Reduce(function(...) merge(...,all.x=T,by=c("ID","CHROM","POS")),list(lrs.gwas.m[,c("ID","CHROM","POS","REF","ALT","REF_FREQ_M","BETA_M","SE_M","P_BOLT_LMM_INF_M")],lrs.gwas.f[,c("ID","CHROM","POS","REF_FREQ_F","BETA_F","SE_F","P_BOLT_LMM_INF_F")],lrs.gwas.m.perm[,c("ID","CHROM","POS","REF_FREQ_M_PERM","BETA_M_PERM","SE_M_PERM","P_BOLT_LMM_INF_M_PERM")],lrs.gwas.f.perm[,c("ID","CHROM","POS","REF_FREQ_F_PERM","BETA_F_PERM","SE_F_PERM","P_BOLT_LMM_INF_F_PERM")]))
#Annotations
effects <- fread(paste0(dir,"/ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_refcorrected.ann.effect.gene.txt"),h=T)
names(effects) <- c("CHROM","POS","ID","REF","ALT","EFFECT","GENE")
mf.gwases.lrs <- merge(mf.gwases.lrs,effects,by=c("ID","CHROM","POS"),all.x = T)
rm(effects)
#Remove duplicate REF/ALT in auto.frqs; duplicates occur because REF/ALT order is not always the same between mf.gwases.lrs and effects
mf.gwases.lrs <- subset(mf.gwases.lrs,!duplicated(mf.gwases.lrs[,c("ID","CHROM","POS","REF.x","ALT.x")]))
#make sure that the alleles are the same b/w mf.frqs and effects
mf.gwases.lrs <- subset(mf.gwases.lrs,(REF.x==REF.y | REF.x==ALT.y) & (ALT.x==REF.y | ALT.x==ALT.y))
#7,842,458
#For some reason, "Affx-8165091" does not have annotation or positional information. Add it manually
mf.gwases.lrs[mf.gwases.lrs$ID=="Affx-8165091",c("POS")] <- 53011955
mf.gwases.lrs <- mf.gwases.lrs[order(mf.gwases.lrs$CHROM,mf.gwases.lrs$ID),]
#Categorise effects into genic/non-genic using 'intergenic region' annotation
mf.gwases.lrs$GENIC_OR_NOT <- ifelse(mf.gwases.lrs$EFFECT %in% c("intergenic_region"),0,1)
#LD-pruning data
ldpruned <- fread(paste0(dir,"ukb_ldpruning_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_ldpruned_0.2.prune.in"),head=F)
mf.gwases.lrs$LD_PRUNED <- ifelse(mf.gwases.lrs$ID %in% ldpruned$V1,1,0)
mf.gwases.lrs <- subset(mf.gwases.lrs,LD_PRUNED==1)
rm(ldpruned)

#Some additional columns
mf.gwases.lrs$MAF <- mf.frqs$MAF
mf.gwases.lrs$MAF_PERM2 <- mf.frqs$MAF_PERM2
#Correlation between genome-wide SNPs
rho <- cor(subset(lrs.gwas.m,LD_PRUNED=1)$BETA_M,subset(lrs.gwas.f,LD_PRUNED=1)$BETA_F)
rho_perm <- cor(subset(lrs.gwas.m.perm,LD_PRUNED=1)$BETA_M_PERM,subset(lrs.gwas.f.perm,LD_PRUNED=1)$BETA_F_PERM)
rm(lrs.gwas.f)
rm(lrs.gwas.m)
rm(lrs.gwas.f.perm)
rm(lrs.gwas.m.perm)


##Viability Fst, Means####

#Theory, adjusted for locus-specific sample sizes
mf.frqs$EXP_MEAN <- 1/(8*mf.frqs$CT_TOT_M)+1/(8*mf.frqs$CT_TOT_F)
set.seed(123)
mf.frqs$RANDOM_CHISQ <- rchisq(nrow(mf.frqs),df=1)
mf.frqs$FST_VIABILITY_THEORY <- mf.frqs$EXP_MEAN*mf.frqs$RANDOM_CHISQ
#Observed, standardised Fst
mf.frqs$FST_VIABILITY_ST <- mf.frqs$FST_VIABILITY/mf.frqs$EXP_MEAN
#Permuted, standardised Fst
mf.frqs$FST_VIABILITY_PERM2_ST <- mf.frqs$FST_VIABILITY_PERM2/mf.frqs$EXP_MEAN

#Observed
mean(mf.frqs$FST_VIABILITY)
#2.104276e-06 (imputed)

#Theory
mean(mf.frqs$FST_VIABILITY_THEORY)
#2.038787e-06 (imputed)

#Permuted null
mean(mf.frqs$FST_VIABILITY_PERM2)
#2.043062e-06 (imputed)

## Viability Fst, Stats####

# Wilcoxon tests
wilcox.test(mf.frqs$FST_VIABILITY_ST,mf.frqs$RANDOM_CHISQ)
#p<0.001 (imputed, genotyped)
wilcox.test(mf.frqs$FST_VIABILITY_ST,mf.frqs$FST_VIABILITY_PERM2_ST)
#p<0.001 (imputed, genotyped)

# Kolmogorov Smirnov tests
ks.test(mf.frqs$FST_VIABILITY_ST,mf.frqs$RANDOM_CHISQ)
#p<0.001 (imputed, genotyped)
ks.test(mf.frqs$FST_VIABILITY_ST,mf.frqs$FST_VIABILITY_PERM2_ST)
#p<0.001 (imputed, genotyped)

## Chisq.test

#Observed vs theory
theory_quantile99 <- quantile(mf.frqs$RANDOM_CHISQ,0.99)
top_1_obs <- sum(mf.frqs$FST_VIABILITY_ST>=theory_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$RANDOM_CHISQ>=theory_quantile99,na.rm=T)
bottom_99_obs <- sum(mf.frqs$FST_VIABILITY_ST<theory_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$RANDOM_CHISQ<theory_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#14.1% excess in top percentile

#Using permuted as null
perm_quantile99 <- quantile(mf.frqs$FST_VIABILITY_PERM2_ST,0.99,na.rm=T)
top_1_obs <- sum(mf.frqs$FST_VIABILITY_ST>=perm_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$FST_VIABILITY_PERM2_ST>=perm_quantile99,na.rm=T)
bottom_99_obs<- sum(mf.frqs$FST_VIABILITY_ST<perm_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$FST_VIABILITY_PERM2_ST<perm_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#13.7% excess in top percentile


## Viability Fst, P and Q-values####
mf.frqs$FST_VIABILITY_P <- pchisq(mf.frqs$FST_VIABILITY_ST,df=1,lower.tail = F)
mf.frqs$FST_VIABILITY_Q <- p.adjust(mf.frqs$FST_VIABILITY_P,"BH")
min(mf.frqs$FST_VIABILITY_P)
#2.236767e-07
min(mf.frqs$FST_VIABILITY_Q)
#0.1761434 (imputed)


## Viability Fst ~ MAF####

#Observed
cor.test(mf.frqs$FST_VIABILITY_ST,mf.frqs$MAF,method="spearman")
#rho=0.009458865 (p<0.001; imputed)

#Theory
cor.test(mf.frqs$RANDOM_CHISQ,mf.frqs$MAF,method="spearman")
#rho=0.0004919547 (p=0.6139; imputed)

#Permuted
cor.test(mf.frqs$FST_VIABILITY_PERM2_ST,mf.frqs$MAF_PERM2,method="spearman")
#rho=-0.0003922775 (p=0.6874; imputed)


##Gametic Fst, Means####

##Theory, adjusted for locus-specific sample sizes ## TO BE UPDATED
mf.frqs$EXP_MEAN <- (1/8) * ( ((1/mf.frqs$CT_TOT_F)*(1+(mf.frqs$VAR_LRS_F/mf.frqs$MEANSQ_LRS_F)) )+( (1/mf.frqs$CT_TOT_M)*(1+(mf.frqs$VAR_LRS_M/mf.frqs$MEANSQ_LRS_M)) ))
set.seed(123)
mf.frqs$RANDOM_CHISQ <- rchisq(nrow(mf.frqs),df=1)
mf.frqs$FST_GAMETIC_THEORY <- mf.frqs$EXP_MEAN*mf.frqs$RANDOM_CHISQ
#Observed, standardised Fst
mf.frqs$FST_GAMETIC_ST <- mf.frqs$FST_GAMETIC/mf.frqs$EXP_MEAN
#Permuted, standardised Fst
mf.frqs$EXP_MEAN <- (1/8) * ( ((1/mf.frqs$CT_TOT_F)*(1+(mf.frqs$VAR_LRS_PERM2_F/mf.frqs$MEANSQ_LRS_PERM2_F)) )+( (1/mf.frqs$CT_TOT_M)*(1+(mf.frqs$VAR_LRS_PERM2_M/mf.frqs$MEANSQ_LRS_PERM2_M)) ))
mf.frqs$FST_GAMETIC_PERM2_ST <- mf.frqs$FST_GAMETIC_PERM2/mf.frqs$EXP_MEAN

#Observed
mean(mf.frqs$FST_GAMETIC)
#2.974048e-06 (imputed)

#Theory 
mean(mf.frqs$FST_GAMETIC_THEORY)
#2.911841e-06 (imputed)

#Permuted2 (sexes)
mean(mf.frqs$FST_GAMETIC_PERM2)
#2.907955e-06 (imputed)

## Gametic Fst, Stats####

# Wilcoxon tests
wilcox.test(mf.frqs$FST_GAMETIC,mf.frqs$FST_GAMETIC_THEORY)
#p<0.001 (imputed, genotyped)
wilcox.test(mf.frqs$FST_GAMETIC_ST,mf.frqs$RANDOM_CHISQ)
#p<0.001 (imputed, genotyped)
wilcox.test(mf.frqs$FST_GAMETIC,mf.frqs$FST_GAMETIC_PERM2)
#p<0.001 (imputed, genotyped)
wilcox.test(mf.frqs$FST_GAMETIC_ST,mf.frqs$FST_GAMETIC_PERM2_ST)
#p<0.001 (imputed, genotyped)

# Kolmogorov Smirnov tests
ks.test(mf.frqs$FST_GAMETIC,mf.frqs$FST_GAMETIC_THEORY)
#p<0.001 (imputed, genotyped)
ks.test(mf.frqs$FST_GAMETIC_ST,mf.frqs$RANDOM_CHISQ)
#p<0.001 (imputed, genotyped)
ks.test(mf.frqs$FST_GAMETIC,mf.frqs$FST_GAMETIC_PERM2)
#p<0.001 (imputed, genotyped)
ks.test(mf.frqs$FST_GAMETIC_ST,mf.frqs$FST_GAMETIC_PERM2_ST)
#p<0.001 (imputed, genotyped)


## Chisq.test

#Observed vs theory
theory_quantile99 <- quantile(mf.frqs$RANDOM_CHISQ,0.99)
top_1_obs <- sum(mf.frqs$FST_GAMETIC_ST>=theory_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$RANDOM_CHISQ>=theory_quantile99,na.rm=T)
bottom_99_obs <- sum(mf.frqs$FST_GAMETIC_ST<theory_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$RANDOM_CHISQ<theory_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#9.0% excess in top percentile

#Using permuted as null
perm_quantile99 <- quantile(mf.frqs$FST_GAMETIC_PERM2_ST,0.99,na.rm=T)
top_1_obs <- sum(mf.frqs$FST_GAMETIC_ST>=perm_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$FST_GAMETIC_PERM2_ST>=perm_quantile99,na.rm=T)
bottom_99_obs<- sum(mf.frqs$FST_GAMETIC_ST<perm_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$FST_GAMETIC_PERM2_ST<perm_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#7.8% excess in top percentile


## Gametic Fst, P and Q-values####

#Observed
mf.frqs$FST_GAMETIC_P <- pchisq(mf.frqs$FST_GAMETIC_ST,df=1,lower.tail = F)
mf.frqs$FST_GAMETIC_Q <- p.adjust(mf.frqs$FST_GAMETIC_P,"BH")
min(mf.frqs$FST_GAMETIC_P)
#4.152214e-06
min(mf.frqs$FST_GAMETIC_Q)
#0.8209188 (imputed)

## Gametic Fst ~ MAF####

#Observed
cor.test(mf.frqs$FST_GAMETIC_ST,mf.frqs$MAF,method="spearman")
#rho=0.006730531 (p<0.001; imputed)

#Theory
cor.test(mf.frqs$RANDOM_CHISQ,mf.frqs$MAF,method="spearman")
#rho=0.001 (p<0.001; imputed)

#Permuted
cor.test(mf.frqs$FST_GAMETIC_PERM2_ST,mf.frqs$MAF_PERM2,method="spearman")
#rho=-0.001593672 (p<0.001; imputed)



## Reprod. Fst, Means####

##Theory, adjusted for locus-specific sample sizes, but no HWE correction
#mf.frqs$EXP_MEAN <- (1/8) * ( ( (1/mf.frqs$CT_TOT_F)*(mf.frqs$VAR_LRS_F/mf.frqs$MEANSQ_LRS_F) ) + ( (1/mf.frqs$CT_TOT_M)*(mf.frqs$VAR_LRS_M/mf.frqs$MEANSQ_LRS_M) ) )
#set.seed(123)
#mf.frqs$RANDOM_CHISQ <- rchisq(nrow(mf.frqs),df=1)
#mf.frqs$FST_REPRODUCTIVE_THEORY_APPROX <- mf.frqs$EXP_MEAN*mf.frqs$RANDOM_CHISQ

##Theory, adjusted for locus-specific sample sizes, with HWE correction
mf.frqs$EXP_MEAN <-(((mf.frqs$p_F_VIABILITY*(1-mf.frqs$p_F_VIABILITY)/(2*mf.frqs$CT_TOT_F))*(mf.frqs$VAR_LRS_F/mf.frqs$MEANSQ_LRS_F)*(1-mf.frqs$FIS_F))+ ((mf.frqs$p_M_VIABILITY*(1-mf.frqs$p_M_VIABILITY)/(2*mf.frqs$CT_TOT_M))*(mf.frqs$VAR_LRS_M/mf.frqs$MEANSQ_LRS_M)*(1-mf.frqs$FIS_M)))/(4*((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2)*(1-((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2)))
set.seed(123)
mf.frqs$RANDOM_CHISQ <- rchisq(nrow(mf.frqs),df=1)
mf.frqs$FST_REPRODUCTIVE_THEORY <- mf.frqs$EXP_MEAN*mf.frqs$RANDOM_CHISQ

#Observed
mf.frqs$FST_REPRODUCTIVE <- ((mf.frqs$p_F-mf.frqs$p_M)-(mf.frqs$p_F_VIABILITY-mf.frqs$p_M_VIABILITY))^2 / (4* ((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2) * (1- ((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2)))
#Observed, standardised
mf.frqs$FST_REPRODUCTIVE_ST <- mf.frqs$FST_REPRODUCTIVE/mf.frqs$EXP_MEAN

#Permuted (offspring)
mf.frqs$FST_REPRODUCTIVE_PERM <- ((mf.frqs$p_F_PERM-mf.frqs$p_M_PERM)-(mf.frqs$p_F_VIABILITY-mf.frqs$p_M_VIABILITY))^2 / (4* ((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2) * (1- ((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2)))
#Permuted (offspring), standardised
mf.frqs$EXP_MEAN <- (((mf.frqs$p_F_VIABILITY*(1-mf.frqs$p_F_VIABILITY)/(2*mf.frqs$CT_TOT_F))*(mf.frqs$VAR_LRS_F/mf.frqs$MEANSQ_LRS_F)*(1-mf.frqs$FIS_F))+ ((mf.frqs$p_M_VIABILITY*(1-mf.frqs$p_M_VIABILITY)/(2*mf.frqs$CT_TOT_M))*(mf.frqs$VAR_LRS_M/mf.frqs$MEANSQ_LRS_M)*(1-mf.frqs$FIS_M)))/(4*((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2)*(1-((mf.frqs$p_M_VIABILITY+mf.frqs$p_F_VIABILITY)/2)))
mf.frqs$FST_REPRODUCTIVE_PERM_ST <- mf.frqs$FST_REPRODUCTIVE_PERM/mf.frqs$EXP_MEAN

#Observed
mean(mf.frqs$FST_REPRODUCTIVE)
#8.900171e-07 (imputed)

#Theory
mean(mf.frqs$FST_REPRODUCTIVE_THEORY)
#8.730574e-07 (imputed)

#Permuted
mean(mf.frqs$FST_REPRODUCTIVE_PERM)
#8.748456e-07 (imputed)

## Reprod. Fst, Stats####

# Wilcoxon tests
wilcox.test(mf.frqs$FST_REPRODUCTIVE_ST,mf.frqs$RANDOM_CHISQ)
#p<0.001 (imputed, genotyped)
wilcox.test(mf.frqs$FST_REPRODUCTIVE_ST,mf.frqs$FST_REPRODUCTIVE_PERM_ST)
#p<0.001 (imputed, genotyped)

# Kolmogorov Smirnov tests
ks.test(mf.frqs$FST_REPRODUCTIVE_ST,mf.frqs$RANDOM_CHISQ)
#p<0.001 (imputed, genotyped)
ks.test(mf.frqs$FST_REPRODUCTIVE_ST,mf.frqs$FST_REPRODUCTIVE_PERM_ST)
#p<0.001 (imputed, genotyped)

##99% quantile

#Observed vs theory
theory_quantile99 <- quantile(mf.frqs$RANDOM_CHISQ,0.99)
top_1_obs <- sum(mf.frqs$FST_REPRODUCTIVE_ST>=theory_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$RANDOM_CHISQ>=theory_quantile99,na.rm=T)
bottom_99_obs <- sum(mf.frqs$FST_REPRODUCTIVE_ST<theory_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$RANDOM_CHISQ<theory_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed), p<0.001 (genotyped)
top_1_obs/top_1_exp
#7.4% excess in top percentile

#Using permuted as null
perm_quantile99 <- quantile(mf.frqs$FST_REPRODUCTIVE_PERM_ST,0.99,na.rm=T)
top_1_obs <- sum(mf.frqs$FST_REPRODUCTIVE_ST>=perm_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$FST_REPRODUCTIVE_PERM_ST>=perm_quantile99,na.rm=T)
bottom_99_obs<- sum(mf.frqs$FST_REPRODUCTIVE_ST<perm_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$FST_REPRODUCTIVE_PERM_ST<perm_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed), p<0.001 (genotyped)
top_1_obs/top_1_exp
#5.0% excess in top percentile



## Reprod. Fst, P and Q-values####

#Observed
mf.frqs$FST_REPRODUCTIVE_P <- pchisq(mf.frqs$FST_REPRODUCTIVE_ST,df=1,lower.tail = F)
mf.frqs$FST_REPRODUCTIVE_Q <- p.adjust(mf.frqs$FST_REPRODUCTIVE_P,"BH")
min(mf.frqs$FST_REPRODUCTIVE_P)
#3.925391e-07
min(mf.frqs$FST_REPRODUCTIVE_Q)
#0.4129311 (imputed)


#Reprod. Fst ~ MAF####

#Observed
cor.test(mf.frqs$FST_REPRODUCTIVE_ST,mf.frqs$MAF,method="spearman")
#rho=0.005903227 (p<0.001; imputed)

#Theory
cor.test(mf.frqs$RANDOM_CHISQ,mf.frqs$MAF,method="spearman")
#rho=0.001 (p<0.001; imputed)

#Permuted
cor.test(mf.frqs$FST_REPRODUCTIVE_PERM_ST,mf.frqs$MAF,method="spearman")
#rho=0.006 (p<0.001; imputed)



## Unfolded Fst, Means####

#Observed 
#Delta1 = (pm'-pm)
mf.frqs$DELTA1 <- (mf.frqs$p_M-mf.frqs$p_M_VIABILITY)/sqrt( ( ( (mf.frqs$p_M_VIABILITY*(1-mf.frqs$p_M_VIABILITY)) / (2*mf.frqs$CT_TOT_M) ) * (mf.frqs$VAR_LRS_M/mf.frqs$MEANSQ_LRS_M) * (1-mf.frqs$FIS_M) ) )
#Delta2 = (pf'-pf)
mf.frqs$DELTA2 <- (mf.frqs$p_F-mf.frqs$p_F_VIABILITY)/sqrt(  ( ((mf.frqs$p_F_VIABILITY*(1-mf.frqs$p_F_VIABILITY))/(2*mf.frqs$CT_TOT_F)) * (mf.frqs$VAR_LRS_F/mf.frqs$MEANSQ_LRS_F) * (1-mf.frqs$FIS_F) ) )
mf.frqs$UNFOLDED_FST <- mf.frqs$DELTA1*mf.frqs$DELTA2

#Permuted
mf.frqs$DELTA1_PERM <- (mf.frqs$p_M_PERM-mf.frqs$p_M_VIABILITY)/sqrt( ( ( (mf.frqs$p_M_VIABILITY*(1-mf.frqs$p_M_VIABILITY)) / (2*mf.frqs$CT_TOT_M) ) * (mf.frqs$VAR_LRS_M/mf.frqs$MEANSQ_LRS_M) * (1-mf.frqs$FIS_M) ) )
mf.frqs$DELTA2_PERM <- (mf.frqs$p_F_PERM-mf.frqs$p_F_VIABILITY)/sqrt(  ( ((mf.frqs$p_F_VIABILITY*(1-mf.frqs$p_F_VIABILITY))/(2*mf.frqs$CT_TOT_F)) * (mf.frqs$VAR_LRS_F/mf.frqs$MEANSQ_LRS_F) * (1-mf.frqs$FIS_F) ) )
mf.frqs$UNFOLDED_FST_PERM <- mf.frqs$DELTA1_PERM*mf.frqs$DELTA2_PERM

#Theory
set.seed(123)
mf.frqs$UNFOLDED_FST_THEORY <- rnorm(n=nrow(mf.frqs),mean=0, sd=1)*rnorm(n=nrow(mf.frqs),mean=0, sd=1)

#Observed, negative
mean(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0])
#-0.6505415 
#Theory
mean(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY<0])
#-0.6354981
#Permuted
mean(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM<0])
#-0.6382583

#Observed, positive
mean(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0])
#0.6943363 
#Theory
mean(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY>0])
#0.6369879 
#Permuted
mean(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM>0])
#0.6401527

## Unfolded Fst, Stats####

# Wilcoxon and KS tests
#Negative unfolded
wilcox.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0],mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY<0])
#p<0.001 (imputed)
wilcox.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0],mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM<0])
#p<0.001 (imputed)
ks.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0],mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY<0])
#p<0.001 (imputed)
ks.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0],mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM<0])
#p<0.001 (imputed)

# Wilcoxon and KS tests
#Positive unfolded
wilcox.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0],mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY>0])
#p<0.001 (imputed, genotyped)
wilcox.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0],mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM>0])
#p<0.001 (imputed, genotyped)
ks.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0],mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY>0])
#p<0.001 (imputed, genotyped)
ks.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0],mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM>0])
#p<0.001 (imputed, genotyped)


##99% quantile

#Positive
#Observed vs theory
theory_quantile99 <- quantile(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY>0],0.99)
top_1_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0]>=theory_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY>0]>=theory_quantile99,na.rm=T)
bottom_99_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0]<theory_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY>0]<theory_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#48.3% excess in top percentile

#Using permuted as null
perm_quantile99 <- quantile(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM>0],0.99,na.rm=T)
top_1_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0]>=perm_quantile99,na.rm=T)
top_1_exp <- sum(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM>0]>=perm_quantile99,na.rm=T)
bottom_99_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0]<perm_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM>0]<perm_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#45.2% excess in top percentile



#Negative
#Observed vs theory
theory_quantile99 <- quantile(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY<0],0.01)
bottom_1_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0]<=theory_quantile99,na.rm=T)
bottom_1_exp <- sum(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY<0]<=theory_quantile99,na.rm=T)
top_99_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0]>theory_quantile99,na.rm=T)
top_99_exp<- sum(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY<0]>theory_quantile99,na.rm=T)
chisq.test(cbind(c(bottom_1_obs,bottom_1_exp),c(top_99_obs,top_99_exp)))
#p<0.001 (imputed)
bottom_1_obs/bottom_1_exp
#9.6% excess in top percentile

#Using permuted as null
perm_quantile99 <- quantile(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM<0],0.01,na.rm=T)
bottom_1_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0]<=perm_quantile99,na.rm=T)
bottom_1_exp <- sum(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM<0]<=perm_quantile99,na.rm=T)
top_99_obs <- sum(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0]>perm_quantile99,na.rm=T)
top_99_exp<- sum(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM<0]>perm_quantile99,na.rm=T)
chisq.test(cbind(c(bottom_1_obs,bottom_1_exp),c(top_99_obs,top_99_exp)))
#p<0.001 (imputed)
bottom_1_obs/bottom_1_exp
#4.2% excess in top percentile



#Unfolded Fst ~ MAF####

#Observed
cor.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST<0],mf.frqs$MAF[mf.frqs$UNFOLDED_FST<0],method="spearman")
#rho=-0.01162673 (p<0.001; imputed);\ 
#Theory
cor.test(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY<0],mf.frqs$MAF[mf.frqs$UNFOLDED_FST_THEORY<0],method="spearman")
#rho=-0.001593896 (p=0.2475; imputed); 
#Permuted
cor.test(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM<0],mf.frqs$MAF[mf.frqs$UNFOLDED_FST_PERM<0],method="spearman")
#rho=-0.0006727542 (p=0.6255; imputed);

cor.test(mf.frqs$UNFOLDED_FST[mf.frqs$UNFOLDED_FST>0],mf.frqs$MAF[mf.frqs$UNFOLDED_FST>0],method="spearman")
#rho=0.03164659 (p<0.001; imputed);
cor.test(mf.frqs$UNFOLDED_FST_THEORY[mf.frqs$UNFOLDED_FST_THEORY>0],mf.frqs$MAF[mf.frqs$UNFOLDED_FST_THEORY>0],method="spearman")
#rho=0.0005189911 (p=0.7064; imputed)
cor.test(mf.frqs$UNFOLDED_FST_PERM[mf.frqs$UNFOLDED_FST_PERM>0],mf.frqs$MAF[mf.frqs$UNFOLDED_FST_PERM>0],method="spearman")
#rho=0.0001074792 (p=0.9378; imputed)


##Lst, Means####

#Absolute odds ratio, transformed to be MAF-independent under the null
mf.gwases$LST <- (sqrt(mf.gwases$REF_FREQ*(1-mf.gwases$REF_FREQ))*mf.gwases$BETA)^2
mf.gwases$LST_PERM2 <- (sqrt(mf.gwases$REF_FREQ_PERM2*(1-mf.gwases$REF_FREQ_PERM2))*mf.gwases$BETA_PERM2)^2

#Observed
mean(mf.gwases$LST)
#5.323197e-07 (imputed)

#Permuted
mean(mf.gwases$LST_PERM2)
#5.235572e-07 (imputed)

## Lst, Stats####
#Statistical tests
wilcox.test(mf.gwases$LST,mf.gwases$LST_PERM2)
#p<0.001 (imputed)
ks.test(mf.gwases$LST,mf.gwases$LST_PERM2)
#p<0.001 (imputed)

#99% quantile
#Using permuted as theory
perm_quantile99 <- quantile(mf.gwases$LST_PERM2,0.99)
top_1_obs <- sum(mf.gwases$LST>=perm_quantile99)
top_1_exp <- sum(mf.gwases$LST_PERM2>=perm_quantile99)
bottom_99_obs <- sum(mf.gwases$LST<perm_quantile99)
bottom_99_exp<- sum(mf.gwases$LST_PERM2<perm_quantile99)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#8.9% excess

##Lst, P- and Q-values####

#Observed
mf.gwases$Q_BOLT_LMM_INF <- p.adjust(mf.gwases$P_BOLT_LMM_INF,"BH")
min(mf.gwases$P_BOLT_LMM_INF,na.rm=T)
min(mf.gwases$Q_BOLT_LMM_INF,na.rm=T)
#0.097 (imputed); 0.145 (genotyped)



## Lst ~ MAF####

#Observed
cor.test(mf.gwases$LST,mf.gwases$MAF,method="spearman")
#rho=0.001207709 (p=0.215; imputed)

#Permuted
cor.test(mf.gwases$LST_PERM2,mf.gwases$MAF,method="spearman")
#rho=-0.010 (p<0.001; imputed)




##|t|, Means####

#t-statistic
mf.gwases.lrs$ABS_T <- abs(with(mf.gwases.lrs,(BETA_F-BETA_M)/sqrt(SE_M^2+SE_F^2-(2*rho*SE_M*SE_F))))
mf.gwases.lrs$ABS_T_PERM <- abs(with(mf.gwases.lrs,(BETA_F_PERM-BETA_M_PERM)/sqrt(SE_M_PERM^2+SE_F_PERM^2-(2*rho_perm*SE_M_PERM*SE_F_PERM))))

#Observed
mean(mf.gwases.lrs$ABS_T)
# 0.8110546 (imputed)

#Permuted
mean(mf.gwases.lrs$ABS_T_PERM)
#0.7963171 (imputed)

## |t|, Stats####

## Statistical tests

#Compare perm to observed
wilcox.test(mf.gwases.lrs$ABS_T,mf.gwases.lrs$ABS_T_PERM)
#p<0.001 (imputed)
ks.test(mf.gwases.lrs$ABS_T,mf.gwases.lrs$ABS_T_PERM)
#P<0.001 (imputed)


#Using permuted as theory
perm_quantile99 <- quantile(mf.gwases.lrs$ABS_T_PERM,0.99,na.rm=T)
top_1_obs <- sum(mf.gwases.lrs$ABS_T>=perm_quantile99,na.rm=T)
top_1_exp <- sum(mf.gwases.lrs$ABS_T_PERM>=perm_quantile99,na.rm=T)
bottom_99_obs<- sum(mf.gwases.lrs$ABS_T<perm_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.gwases.lrs$ABS_T_PERM<perm_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
#p<0.001 (imputed)
top_1_obs/top_1_exp
#11.9% excess


##|t| ~ MAF####

#Observed
cor.test(mf.gwases.lrs$ABS_T,mf.gwases.lrs$MAF,method="spearman")
#rho=0.005180973 (p<0.001; imputed)

#Permuted
cor.test(mf.gwases.lrs$ABS_T_PERM,mf.gwases.lrs$MAF,method="spearman")
#rho=-0.0004088011 (p<0.001; imputed)



##Unfolded t, Means####

#t-statistic
mf.gwases.lrs$UNFOLDED_T <- with(mf.gwases.lrs,(BETA_F*BETA_M)/sqrt(SE_M^2*SE_F^2))
mf.gwases.lrs$UNFOLDED_T_PERM <- with(mf.gwases.lrs,(BETA_F_PERM*BETA_M_PERM)/sqrt(SE_M_PERM^2*SE_F_PERM^2))

#Observed
mean(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T>0])
#0.6921113 (imputed)
mean(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T<0])
#-0.6486949 (imputed)


#Permuted
mean(mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM>0])
#0.6392079 (imputed)
mean(mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM<0])
#-0.6388935 (imputed)

## Unfolded t, Stats####

## Statistical tests

#Compare perm to observed
wilcox.test(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T>0],mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM>0])
#p<0.001 (imputed)
ks.test(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T>0],mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM>0])
#P<0.001 (imputed)

wilcox.test(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T<0],mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM<0])
#p<0.001 (imputed)
ks.test(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T<0],mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM<0])
#P<0.001 (imputed)

#Using permuted as theory
perm_quantile99 <- quantile(mf.gwases.lrs$ABS_T_PERM,0.99,na.rm=T)
top_1_obs <- sum(mf.gwases.lrs$ABS_T>=perm_quantile99,na.rm=T)
top_1_exp <- sum(mf.gwases.lrs$ABS_T_PERM>=perm_quantile99,na.rm=T)
bottom_99_obs<- sum(mf.gwases.lrs$ABS_T<perm_quantile99,na.rm=T)
bottom_99_exp<- sum(mf.gwases.lrs$ABS_T_PERM<perm_quantile99,na.rm=T)
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))$expected
chisq.test(cbind(c(top_1_obs,top_1_exp),c(bottom_99_obs,bottom_99_exp)))$observed
#p<0.001 (imputed), p<0.001 (genotyped)







##Unfolded t ~ MAF####

#Observed
cor.test(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T>0],mf.gwases.lrs$MAF[mf.gwases.lrs$UNFOLDED_T>0],method="spearman")
cor.test(mf.gwases.lrs$UNFOLDED_T[mf.gwases.lrs$UNFOLDED_T<0],mf.gwases.lrs$MAF[mf.gwases.lrs$UNFOLDED_T<0],method="spearman")

#Permuted
cor.test(mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM>0],mf.gwases.lrs$MAF[mf.gwases.lrs$UNFOLDED_T_PERM>0],method="spearman")
cor.test(mf.gwases.lrs$UNFOLDED_T_PERM[mf.gwases.lrs$UNFOLDED_T_PERM<0],mf.gwases.lrs$MAF[mf.gwases.lrs$UNFOLDED_T_PERM<0],method="spearman")
