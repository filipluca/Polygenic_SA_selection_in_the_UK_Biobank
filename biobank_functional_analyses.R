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

#LD SCORE REGRESSION####

#Compatible alleles, prep####

#Make file containing list of all LD-pruned imputed SNPs
ruzickaetal.snplist <- mf.frqs[,c("ID","REF.x","ALT.x")]
names(ruzickaetal.snplist) <- c("SNP","A1","A2")
#write.table(ruzickaetal.snplist,paste0(dir,"/ukb_LDSC_v3/compatible_alleles/ruzickaetal.snplist"),quote=F,sep=" ",col.names=F,row.names=F)

#Make file containing header of .sumstats file
header.sumstats <- c("SNP","N","A1","A2","Z")
#write.table(header.sumstats,paste0(dir,"/ukb_LDSC_v3/compatible_alleles/header.sumstats"),quote=F,sep=" ",col.names=F,row.names=F)

#Import ####
#Objects to import are: mf.frqs, mf.gwases and mf.gwases.lrs
#Make sure they are not LD-pruned before import

#Additional columns, needed for LDSC
mf.frqs$N <- mf.frqs$CT_TOT_F+mf.frqs$CT_TOT_M
mf.frqs$N_PERM2 <- mf.frqs$CT_TOT_F_PERM2+mf.frqs$CT_TOT_M_PERM2
mf.gwases$N <- mf.frqs$CT_TOT_F+mf.frqs$CT_TOT_M
mf.gwases$N_PERM2 <- mf.frqs$CT_TOT_F_PERM2+mf.frqs$CT_TOT_M_PERM2
mf.gwases.lrs$N <- mf.frqs$CT_TOT_F+mf.frqs$CT_TOT_M




#Import compatible alleles
#Import list of compatible alleles (i.e. positions where reference and alternative match with the HapMap 3 snp list that LDSC uses for partitioned heritability and genetic correlations)
#See UNIX code to see how this was made
compatible.alleles <- fread(paste0(dir,"/ukb_LDSC_v3/compatible_alleles/ruzickaetal.snplist.compatible"),h=T)
names(compatible.alleles) <- c("SNP","A1","A2")

#Subset imputed dataframes to only include compatible alleles
mf.frqs.comp <- subset(mf.frqs,ID %in% compatible.alleles$SNP)
mf.gwases.comp <- subset(mf.gwases,ID %in% compatible.alleles$SNP)
mf.gwases.lrs.comp <- subset(mf.gwases.lrs,ID %in% compatible.alleles$SNP)

#Adult Fst -> Z####

#Zscore = sqrt(Chi-sq stat) = sqrt (standardised Fst)
mf.frqs.comp$Z_VIABILITY <- sqrt(mf.frqs.comp$FST_VIABILITY_ST)
mf.frqs.comp$Z_VIABILITY_PERM2 <- sqrt(mf.frqs.comp$FST_VIABILITY_PERM2_ST)

#Is pf>pm? If yes, then z score is multiplied by 1, else -1.
mf.frqs.comp$Z_VIABILITY <- ifelse(mf.frqs.comp$p_F_VIABILITY>mf.frqs.comp$p_M_VIABILITY,mf.frqs.comp$Z_VIABILITY,mf.frqs.comp$Z_VIABILITY*-1)
mf.frqs.comp$Z_VIABILITY_PERM2 <- ifelse(mf.frqs.comp$p_F_VIABILITY_PERM2>mf.frqs.comp$p_M_VIABILITY_PERM2,mf.frqs.comp$Z_VIABILITY_PERM2,mf.frqs.comp$Z_VIABILITY_PERM2*-1)

#Keep only relevant columns, recode REF/ALT as A1/A2, and export them to a .sumstats file 
#REF is coded the "effect" allele (i.e., higher value of REF than ALT for LRS will lead to a positive Z; higher frequency of REF than ALT in females [coded as "2" while males are coded as "1"] leads to a positive Z). The effect allele is A1 in the LDSC notation, while A2 is the non-effect allele

#Observed
mf.frqs.ldsc.adult.obs <- mf.frqs.comp[,c("ID","N","REF.x","ALT.x","Z_VIABILITY")]
names(mf.frqs.ldsc.adult.obs) <- c("SNP","N","A1","A2","Z")
##write.table(mf.frqs.ldsc.adult.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Adult_Fst.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.adult.obs)

#Permuted
mf.frqs.ldsc.adult.perm2 <- mf.frqs.comp[,c("ID","N_PERM2","REF.x","ALT.x","Z_VIABILITY_PERM2")]
names(mf.frqs.ldsc.adult.perm2) <- c("SNP","N","A1","A2","Z")
##write.table(mf.frqs.ldsc.adult.perm2,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Adult_Fst_perm2.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.adult.perm2)


#Gametic Fst -> Z####

#Zscore = sqrt(Chi-sq stat) = sqrt (standardised Fst)
mf.frqs.comp$Z_GAMETIC <- sqrt(mf.frqs.comp$FST_GAMETIC_ST)
mf.frqs.comp$Z_GAMETIC_PERM2 <- sqrt(mf.frqs.comp$FST_GAMETIC_PERM2_ST)

#Is pf'>pm'? If yes, then z score is multiplied by 1, else -1. 
mf.frqs.comp$Z_GAMETIC <- ifelse(mf.frqs.comp$p_F>mf.frqs.comp$p_M,mf.frqs.comp$Z_GAMETIC,mf.frqs.comp$Z_GAMETIC*-1)
mf.frqs.comp$Z_GAMETIC_PERM2 <- ifelse(mf.frqs.comp$p_F_PERM2>mf.frqs.comp$p_M_PERM2,mf.frqs.comp$Z_GAMETIC_PERM2,mf.frqs.comp$Z_GAMETIC_PERM2*-1)

#Keep only relevant columns, recode REF/ALT as A1/A2, and export them to a .sumstats file 
#REF is coded the "effect" allele (i.e., higher value of REF than ALT for LRS will lead to a positive Z; higher frequency of REF than ALT in females [coded as "2" while males are coded as "1"] leads to a positive Z). The effect allele is A1 in the LDSC notation, while A2 is the non-effect allele

#Observed
mf.frqs.ldsc.gametic.obs <- mf.frqs.comp[,c("ID","N","REF.x","ALT.x","Z_GAMETIC")]
names(mf.frqs.ldsc.gametic.obs) <- c("SNP","N","A1","A2","Z")
#write.table(mf.frqs.ldsc.gametic.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Gametic_Fst.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.gametic.obs)

#Permuted
mf.frqs.ldsc.gametic.perm2 <- mf.frqs.comp[,c("ID","N_PERM2","REF.x","ALT.x","Z_GAMETIC_PERM2")]
names(mf.frqs.ldsc.gametic.perm2) <- c("SNP","N","A1","A2","Z")
#write.table(mf.frqs.ldsc.gametic.perm2,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Gametic_Fst_perm2.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.gametic.perm2)


## Reprod. Fst -> Z####

#Zscore = sqrt(Chi-sq stat) = sqrt (standardised Fst)
mf.frqs.comp$Z_REPRODUCTIVE <- sqrt(mf.frqs.comp$FST_REPRODUCTIVE_ST)
mf.frqs.comp$Z_REPRODUCTIVE_PERM <- sqrt(mf.frqs.comp$FST_REPRODUCTIVE_PERM_ST)

#Is (pf'-pf)-(pm'-pm)>0? If yes, then z score is multiplied by 1, else -1. 
mf.frqs.comp$Z_REPRODUCTIVE <- with(mf.frqs.comp,ifelse( ((p_F-p_F_VIABILITY)-(p_M-p_M_VIABILITY))>0,mf.frqs.comp$Z_REPRODUCTIVE,mf.frqs.comp$Z_REPRODUCTIVE*-1))
mf.frqs.comp$Z_REPRODUCTIVE_PERM <- with(mf.frqs.comp,ifelse( ((p_F_PERM-p_F_VIABILITY)-(p_M_PERM-p_M_VIABILITY))>0,mf.frqs.comp$Z_REPRODUCTIVE_PERM,mf.frqs.comp$Z_REPRODUCTIVE_PERM*-1))

#Keep only relevant columns, recode REF/ALT as A1/A2, and export them to a .sumstats file 
#REF is coded the "effect" allele (i.e., higher value of REF than ALT for LRS will lead to a positive Z; higher frequency of REF than ALT in females [coded as "2" while males are coded as "1"] leads to a positive Z). The effect allele is A1 in the LDSC notation, while A2 is the non-effect allele
#Observed
mf.frqs.ldsc.reproductive.obs <- mf.frqs.comp[,c("ID","N","REF.x","ALT.x","Z_REPRODUCTIVE")]
names(mf.frqs.ldsc.reproductive.obs) <- c("SNP","N","A1","A2","Z")
#write.table(mf.frqs.ldsc.reproductive.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Reproductive_Fst.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.reproductive.obs)

mf.frqs.ldsc.reproductive.obs <- mf.frqs.comp[,c("ID","N","REF.x","ALT.x","FST_REPRODUCTIVE_ST")]
names(mf.frqs.ldsc.reproductive.obs) <- c("SNP","N","A1","A2","Z")
#write.table(mf.frqs.ldsc.reproductive.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Reproductive_Fst.2.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.reproductive.obs)

#Permuted
mf.frqs.ldsc.reproductive.perm <- mf.frqs.comp[,c("ID","N","REF.x","ALT.x","Z_REPRODUCTIVE_PERM")]
names(mf.frqs.ldsc.reproductive.perm) <- c("SNP","N","A1","A2","Z")
#write.table(mf.frqs.ldsc.reproductive.perm,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Reproductive_Fst_perm.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.reproductive.perm)



## Unfolded Fst -> Z####


#Zscore = Fst-mean(Fst)/sqrt(Fst)
mf.frqs.comp$Z_UNFOLDED[mf.frqs.comp$UNFOLDED_FST>0] <- (mf.frqs.comp$UNFOLDED_FST[mf.frqs.comp$UNFOLDED_FST>0]-mean(mf.frqs.comp$UNFOLDED_FST[mf.frqs.comp$UNFOLDED_FST>0]))/sqrt(var(mf.frqs.comp$UNFOLDED_FST[mf.frqs.comp$UNFOLDED_FST>0]))
mf.frqs.comp$Z_UNFOLDED[mf.frqs.comp$UNFOLDED_FST<0] <- (mf.frqs.comp$UNFOLDED_FST[mf.frqs.comp$UNFOLDED_FST<0]-mean(mf.frqs.comp$UNFOLDED_FST[mf.frqs.comp$UNFOLDED_FST<0]))/sqrt(var(mf.frqs.comp$UNFOLDED_FST[mf.frqs.comp$UNFOLDED_FST<0]))
mf.frqs.comp$Z_UNFOLDED_PERM[mf.frqs.comp$UNFOLDED_FST_PERM>0] <- (mf.frqs.comp$UNFOLDED_FST_PERM[mf.frqs.comp$UNFOLDED_FST_PERM>0]-mean(mf.frqs.comp$UNFOLDED_FST_PERM[mf.frqs.comp$UNFOLDED_FST_PERM>0]))/sqrt(var(mf.frqs.comp$UNFOLDED_FST_PERM[mf.frqs.comp$UNFOLDED_FST_PERM>0]))
mf.frqs.comp$Z_UNFOLDED_PERM[mf.frqs.comp$UNFOLDED_FST_PERM<0] <- (mf.frqs.comp$UNFOLDED_FST_PERM[mf.frqs.comp$UNFOLDED_FST_PERM<0]-mean(mf.frqs.comp$UNFOLDED_FST_PERM[mf.frqs.comp$UNFOLDED_FST_PERM<0]))/sqrt(var(mf.frqs.comp$UNFOLDED_FST_PERM[mf.frqs.comp$UNFOLDED_FST_PERM<0]))

#Keep only relevant columns, recode REF/ALT as A1/A2, and export them to a .sumstats file 
#REF is coded the "effect" allele (i.e., higher value of REF than ALT for LRS will lead to a positive Z; higher frequency of REF than ALT in females [coded as "2" while males are coded as "1"] leads to a positive Z). The effect allele is A1 in the LDSC notation, while A2 is the non-effect allele
#Observed
mf.frqs.ldsc.unfolded.pos.obs <- mf.frqs.comp[mf.frqs.comp$UNFOLDED_FST>0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED")]
names(mf.frqs.ldsc.unfolded.pos.obs) <- c("SNP","N","A1","A2","Z")
mf.frqs.ldsc.unfolded.pos.obs <- subset(mf.frqs.ldsc.unfolded.pos.obs,!is.na(mf.frqs.ldsc.unfolded.pos.obs$Z))
#write.table(mf.frqs.ldsc.unfolded.pos.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_Fst_positive.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.unfolded.pos.obs)

mf.frqs.ldsc.unfolded.neg.obs <- mf.frqs.comp[mf.frqs.comp$UNFOLDED_FST<0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED")]
names(mf.frqs.ldsc.unfolded.neg.obs) <- c("SNP","N","A1","A2","Z")
mf.frqs.ldsc.unfolded.neg.obs <- subset(mf.frqs.ldsc.unfolded.neg.obs,!is.na(mf.frqs.ldsc.unfolded.neg.obs$Z))
#write.table(mf.frqs.ldsc.unfolded.neg.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_Fst_negative.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.unfolded.neg.obs)

#Permuted
mf.frqs.ldsc.unfolded.pos.perm <- mf.frqs.comp[mf.frqs.comp$UNFOLDED_FST_PERM>0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED_PERM")]
names(mf.frqs.ldsc.unfolded.pos.perm) <- c("SNP","N","A1","A2","Z")
mf.frqs.ldsc.unfolded.pos.perm <- subset(mf.frqs.ldsc.unfolded.pos.perm,!is.na(mf.frqs.ldsc.unfolded.pos.perm$Z))
#write.table(mf.frqs.ldsc.unfolded.pos.perm,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_Fst_positive_perm.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.unfolded.pos.perm)

mf.frqs.ldsc.unfolded.neg.perm <- mf.frqs.comp[mf.frqs.comp$UNFOLDED_FST_PERM<0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED_PERM")]
names(mf.frqs.ldsc.unfolded.neg.perm) <- c("SNP","N","A1","A2","Z")
mf.frqs.ldsc.unfolded.neg.perm <- subset(mf.frqs.ldsc.unfolded.neg.perm,!is.na(mf.frqs.ldsc.unfolded.neg.perm$Z))
#write.table(mf.frqs.ldsc.unfolded.neg.perm,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_Fst_negative_perm.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.frqs.ldsc.unfolded.neg.perm)



#Lst -> Z####


#Zscore = logodds/SE
#The direction of the effect is such that higher Z values indicate that the reference allele is more frequent in females
mf.gwases.comp$Z_LST <- mf.gwases.comp$BETA/mf.gwases.comp$SE
mf.gwases.comp$Z_LST_PERM2 <- mf.gwases.comp$BETA_PERM2/mf.gwases.comp$SE_PERM2

#Keep only relevant columns, recode REF/ALT as A1/A2, and export them to a .sumstats file 
#REF is coded the "effect" allele (i.e., higher value of REF than ALT for LRS will lead to a positive Z; higher frequency of REF than ALT in females [coded as "2" while males are coded as "1"] leads to a positive Z). The effect allele is A1 in the LDSC notation, while A2 is the non-effect allele

#Observed
mf.gwases.ldsc.lst.obs <- mf.gwases.comp[,c("ID","N","REF.x","ALT.x","Z_LST")]
names(mf.gwases.ldsc.lst.obs) <- c("SNP","N","A1","A2","Z")
#write.table(mf.gwases.ldsc.lst.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Lst.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.ldsc.lst.obs)

#Permuted
mf.gwases.ldsc.lst.perm2 <- mf.gwases.comp[,c("ID","N_PERM2","REF.x","ALT.x","Z_LST_PERM2")]
names(mf.gwases.ldsc.lst.perm2) <- c("SNP","N","A1","A2","Z")
#write.table(mf.gwases.ldsc.lst.perm2,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Lst_perm2.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.ldsc.lst.perm2)




#|t| -> signed t####


#Zscore = logodds/SE
#The direction of the effect is such that higher Z values indicate that the reference allele is more frequent in females
mf.gwases.lrs.comp$T <- with(mf.gwases.lrs.comp,(BETA_F-BETA_M)/sqrt(SE_M^2+SE_F^2-(2*rho*SE_M*SE_F)))
mf.gwases.lrs.comp$T_PERM <- with(mf.gwases.lrs.comp,(BETA_F_PERM-BETA_M_PERM)/sqrt(SE_M_PERM^2+SE_F_PERM^2-(2*rho_perm*SE_M_PERM*SE_F_PERM)))

#Keep only relevant columns, recode REF/ALT as A1/A2, and export them to a .sumstats file 
#REF is coded the "effect" allele (i.e., higher value of REF than ALT for LRS will lead to a positive Z; higher frequency of REF than ALT in females [coded as "2" while males are coded as "1"] leads to a positive Z). The effect allele is A1 in the LDSC notation, while A2 is the non-effect allele

#Observed
mf.gwases.lrs.ldsc.t.obs <- mf.gwases.lrs.comp[,c("ID","N","REF.x","ALT.x","T")]
names(mf.gwases.lrs.ldsc.t.obs) <- c("SNP","N","A1","A2","Z")
#write.table(mf.gwases.lrs.ldsc.t.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_t.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.lrs.ldsc.t.obs)

#Permuted
mf.gwases.lrs.ldsc.t.perm <- mf.gwases.lrs.comp[,c("ID","N","REF.x","ALT.x","T_PERM")]
names(mf.gwases.lrs.ldsc.t.perm) <- c("SNP","N","A1","A2","Z")
#write.table(mf.gwases.lrs.ldsc.t.perm,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_t_perm.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.lrs.ldsc.t.perm)


#Unfolded t -> Z####

#Zscore = Fst-mean(Fst)/sqrt(Fst)
mf.gwases.lrs.comp$Z_UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T>0] <- (mf.gwases.lrs.comp$UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T>0]-mean(mf.gwases.lrs.comp$UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T>0]))/sqrt(var(mf.gwases.lrs.comp$UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T>0]))
mf.gwases.lrs.comp$Z_UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T<0] <- (mf.gwases.lrs.comp$UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T<0]-mean(mf.gwases.lrs.comp$UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T<0]))/sqrt(var(mf.gwases.lrs.comp$UNFOLDED_T[mf.gwases.lrs.comp$UNFOLDED_T<0]))
mf.gwases.lrs.comp$Z_UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM>0] <- (mf.gwases.lrs.comp$UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM>0]-mean(mf.gwases.lrs.comp$UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM>0]))/sqrt(var(mf.gwases.lrs.comp$UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM>0]))
mf.gwases.lrs.comp$Z_UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM<0] <- (mf.gwases.lrs.comp$UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM<0]-mean(mf.gwases.lrs.comp$UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM<0]))/sqrt(var(mf.gwases.lrs.comp$UNFOLDED_T_PERM[mf.gwases.lrs.comp$UNFOLDED_T_PERM<0]))

#Keep only relevant columns, recode REF/ALT as A1/A2, and export them to a .sumstats file 
#REF is coded the "effect" allele (i.e., higher value of REF than ALT for LRS will lead to a positive Z; higher frequency of REF than ALT in females [coded as "2" while males are coded as "1"] leads to a positive Z). The effect allele is A1 in the LDSC notation, while A2 is the non-effect allele
#Observed
mf.gwases.lrs.ldsc.unfolded.t.pos.obs <- mf.gwases.lrs.comp[mf.gwases.lrs.comp$UNFOLDED_T>0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED_T")]
names(mf.gwases.lrs.ldsc.unfolded.t.pos.obs) <- c("SNP","N","A1","A2","Z")
mf.gwases.lrs.ldsc.unfolded.t.pos.obs <- subset(mf.gwases.lrs.ldsc.unfolded.t.pos.obs,!is.na(mf.gwases.lrs.ldsc.unfolded.t.pos.obs$Z))
#write.table(mf.gwases.lrs.ldsc.unfolded.t.pos.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_t_positive.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.lrs.ldsc.unfolded.t.pos.obs)

mf.gwases.lrs.ldsc.unfolded.t.neg.obs <- mf.gwases.lrs.comp[mf.gwases.lrs.comp$UNFOLDED_T<0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED_T")]
names(mf.gwases.lrs.ldsc.unfolded.t.neg.obs) <- c("SNP","N","A1","A2","Z")
mf.gwases.lrs.ldsc.unfolded.t.neg.obs <- subset(mf.gwases.lrs.ldsc.unfolded.t.neg.obs,!is.na(mf.gwases.lrs.ldsc.unfolded.t.neg.obs$Z))
#write.table(mf.gwases.lrs.ldsc.unfolded.t.neg.obs,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_t_negative.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.lrs.ldsc.unfolded.t.neg.obs)

#Permuted
mf.gwases.lrs.ldsc.unfolded.t.pos.perm <- mf.gwases.lrs.comp[mf.gwases.lrs.comp$UNFOLDED_T_PERM>0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED_T_PERM")]
names(mf.gwases.lrs.ldsc.unfolded.t.pos.perm) <- c("SNP","N","A1","A2","Z")
mf.gwases.lrs.ldsc.unfolded.t.pos.perm <- subset(mf.gwases.lrs.ldsc.unfolded.t.pos.perm,!is.na(mf.gwases.lrs.ldsc.unfolded.t.pos.perm$Z))
#write.table(mf.gwases.lrs.ldsc.unfolded.t.pos.perm,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_t_positive_perm.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.lrs.ldsc.unfolded.t.pos.perm)

mf.gwases.lrs.ldsc.unfolded.t.neg.perm <- mf.gwases.lrs.comp[mf.gwases.lrs.comp$UNFOLDED_T_PERM<0,c("ID","N","REF.x","ALT.x","Z_UNFOLDED_T_PERM")]
names(mf.gwases.lrs.ldsc.unfolded.t.neg.perm) <- c("SNP","N","A1","A2","Z")
mf.gwases.lrs.ldsc.unfolded.t.neg.perm <- subset(mf.gwases.lrs.ldsc.unfolded.t.neg.perm,!is.na(mf.gwases.lrs.ldsc.unfolded.t.neg.perm$Z))
#write.table(mf.gwases.lrs.ldsc.unfolded.t.neg.perm,paste0(dir,"/ukb_LDSC_v3/sumstats/PASS_Unfolded_t_negative_perm.sumstats"),quote=F,sep=" ",row.names=F,col.names=T)
rm(mf.gwases.lrs.ldsc.unfolded.t.neg.perm)
