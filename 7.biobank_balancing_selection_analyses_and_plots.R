
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



#ALLELE AGE####

#Import####

#Genotyped
#allele.ages <- fread(paste0(dir,"allele_age_data/atlas.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3_with_header.csv"),skip=3)
#allele.ages <- allele.ages[,c("VariantID","Chromosome","Position","AlleleRef","AlleleAlt","DataSource","AgeMedian_Jnt")]
#names(allele.ages)[1] <- "ID"

#Imputed
allele.ages <- fread(paste0(dir,"allele_age_data/atlas.ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_with_header.csv"),skip=3)
allele.ages <- allele.ages[,c("VariantID","Chromosome","Position","AlleleRef","AlleleAlt","DataSource","AgeMedian_Jnt")]
names(allele.ages)[1] <- "ID"

#Make mf.frqs table smaller, for merging
#mf.frqs <- mf.frqs[,c("CHROM","ID","REF.x","ALT.x","MAF","FST_VIABILITY","FST_VIABILITY_PERM","FST_VIABILITY_THEORY","FST_GAMETIC","FST_GAMETIC_PERM2","FST_GAMETIC_THEORY","FST_REPRODUCTIVE","FST_REPRODUCTIVE_THEORY_APPROX","FST_REPRODUCTIVE_PERM","POS","REF.y","ALT.y","EFFECT","GENE","ALT_FREQ")]

#Merge estimates of allele age with estimates of Fst
#allele.ages <- merge(allele.ages,mf.frqs,by.x="VariantID",by.y="ID",all.y=T)
allele.ages <- Reduce(function(...) merge(...,by="ID",all.x=T),list(mf.frqs,mf.gwases[,c("ID","LST","LST_PERM2")],mf.gwases.lrs[,c("ID","ABS_T","ABS_T_PERM","UNFOLDED_T","UNFOLDED_T_PERM")],allele.ages))

#Remove positions where alternative and reference alleles do not match, or Chromosomes do not match
allele.ages <- subset(allele.ages,(AlleleRef==REF.x | AlleleRef==ALT.x) & (AlleleAlt==REF.x | AlleleAlt==ALT.x) & Chromosome==CHROM)

#Keep only one set of allele ages for each SNP
allele.ages.tgp <- subset(allele.ages,DataSource=="TGP")
allele.ages.sgdp <- subset(allele.ages,DataSource=="SGDP")
allele.ages.comb <- subset(allele.ages,DataSource=="Combined")
rm(allele.ages)
allele.ages2 <- Reduce(function(...) merge(...,by="ID",all=T),list(allele.ages.tgp,allele.ages.sgdp[,c("ID","DataSource","AgeMedian_Jnt")],allele.ages.comb[,c("ID","DataSource","AgeMedian_Jnt")]))
allele.ages2$AgeMedian_Jnt_definitive <- with(allele.ages2,ifelse(!is.na(AgeMedian_Jnt),AgeMedian_Jnt,ifelse(!is.na(AgeMedian_Jnt.x),AgeMedian_Jnt.x,ifelse(!is.na(AgeMedian_Jnt.y),AgeMedian_Jnt.y,NA))))
rm(allele.ages.tgp)
rm(allele.ages.sgdp)
rm(allele.ages.comb)

#Remove sites which are not present in Biobank data
allele.ages2 <- subset(allele.ages2,!is.na(Chromosome) & !is.na(Position))

#What is the frequency of alternative allele (according to Atlas of Variant Age) in the UK Biobank?
allele.ages2$AlleleAlt_Freq <- with(allele.ages2,ifelse(AlleleAlt==ALT.x,1-(((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F)),(((p_M_VIABILITY*CT_TOT_M)+(p_F_VIABILITY*CT_TOT_F))/(CT_TOT_M+CT_TOT_F))))


#Statistical tests####

#Rank correlations
#Viability Fst
cor.test(allele.ages2$AgeMedian_Jnt_definitive/allele.ages2$AlleleAlt_Freq,allele.ages2$FST_VIABILITY_ST,method = "spearman")
#rho=-0.000826759, p=0.4016  

#Reproductive Fst
cor.test(allele.ages2$AgeMedian_Jnt_definitive/allele.ages2$AlleleAlt_Freq,allele.ages2$FST_REPRODUCTIVE_ST,method = "spearman") 
#rho=0.0003841976, p=0.6967

#Gametic Fst
cor.test(allele.ages2$AgeMedian_Jnt_definitive/allele.ages2$AlleleAlt_Freq,allele.ages2$FST_GAMETIC_ST,method = "spearman")
#rho=0.001290932, p=0.1903

#Unfolded Fst
cor.test(allele.ages2$AgeMedian_Jnt_definitive[allele.ages2$UNFOLDED_FST<0]/allele.ages2$AlleleAlt_Freq[allele.ages2$UNFOLDED_FST<0],allele.ages2$UNFOLDED_FST[allele.ages2$UNFOLDED_FST<0],method = "spearman")
#rho=-0.003820654, p=0.006609
cor.test(allele.ages2$AgeMedian_Jnt_definitive[allele.ages2$UNFOLDED_FST>0]/allele.ages2$AlleleAlt_Freq[allele.ages2$UNFOLDED_FST>0],allele.ages2$UNFOLDED_FST[allele.ages2$UNFOLDED_FST>0],method = "spearman")
#rho=-0.005438573,p=8.294e-05

#Lst
cor.test(allele.ages2$AgeMedian_Jnt_definitive/allele.ages2$AlleleAlt_Freq,allele.ages2$LST,method = "spearman") 
#rho=-0.002, p=0.0206

#|t|
cor.test(allele.ages2$AgeMedian_Jnt_definitive/allele.ages2$AlleleAlt_Freq,allele.ages2$ABS_T,method = "spearman") 
#rho=-7.410957e-05, p=0.9401

#Unfolded t
cor.test(allele.ages2$AgeMedian_Jnt_definitive[allele.ages2$UNFOLDED_T<0]/allele.ages2$AlleleAlt_Freq[allele.ages2$UNFOLDED_T<0],allele.ages2$UNFOLDED_T[allele.ages2$UNFOLDED_T<0],method = "spearman")
#rho=-0.002852482, p=0.0426
cor.test(allele.ages2$AgeMedian_Jnt_definitive[allele.ages2$UNFOLDED_T>0]/allele.ages2$AlleleAlt_Freq[allele.ages2$UNFOLDED_T>0],allele.ages2$UNFOLDED_T[allele.ages2$UNFOLDED_T>0],method = "spearman")
#rho=-0.003932987,p=0.00442


#Quantile-proportion prep####

#No alternative allele correction

#Null quantiles
quantiles.fst.theory <- quantile(allele.ages2$RANDOM_CHISQ, prob = seq(0, 1, length = 101), type = 1)
quantiles.fst.theory[1] <- 0
quantiles.fst.theory[101] <- 10^6

quantiles.unfolded.theory <- quantile(allele.ages2$UNFOLDED_FST_THEORY, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.theory[1] <- -10^6
quantiles.unfolded.theory[101] <- 10^6

quantiles.lst.perm <- quantile(allele.ages2$LST_PERM2, prob = seq(0, 1, length = 101), type = 1)
quantiles.lst.perm[1] <- 0
quantiles.lst.perm[101] <- 10^6

quantiles.t.perm <- quantile(allele.ages2$ABS_T_PERM, prob = seq(0, 1, length = 101), type = 1)
quantiles.t.perm[1] <- 0
quantiles.t.perm[101] <- 10^6

quantiles.unfolded.t.perm <- quantile(allele.ages2$UNFOLDED_T, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.t.perm[1] <- -10^6
quantiles.unfolded.t.perm[101] <- 10^6

#Ages
ages.viability <- vector()
ages.viability.perm <- vector()
#ages.viability.theory <- vector()

ages.reproductive <- vector()
ages.reproductive.perm <-vector()
#ages.reproductive.theory <-vector()

ages.total <- vector()
ages.total.perm <- vector()
#ages.total.theory <- vector()

ages.unfolded <- vector()
ages.unfolded.perm <- vector()
#ages.unfolded.theory <- vector()

ages.lst <- vector()
ages.lst.perm <-vector()

ages.t <- vector()
ages.t.perm <- vector()

ages.unfolded.t <- vector()
ages.unfolded.t.perm <- vector()

  for (j in 1:100){
    
    ages.viability[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$FST_VIABILITY_ST>quantiles.fst.theory[j] & allele.ages2$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]])
    ages.viability.perm[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & allele.ages2$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]])
    #ages.viability.theory[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$RANDOM_CHISQ>quantiles.fst.theory[j] & allele.ages2$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]])
    
    ages.reproductive[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & allele.ages2$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]])
    ages.reproductive.perm[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & allele.ages2$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]])
    #ages.reproductive.theory[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$RANDOM_CHISQ>quantiles.fst.theory[j] & allele.ages2$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]])
    
    ages.total[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$FST_GAMETIC_ST>quantiles.fst.theory[j] & allele.ages2$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]])
    ages.total.perm[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & allele.ages2$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]])
    #ages.total.theory[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$RANDOM_CHISQ>quantiles.fst.theory[j] & allele.ages2$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]])
    
    ages.unfolded[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$UNFOLDED_FST>quantiles.unfolded.theory[j] & allele.ages2$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]])
    ages.unfolded.perm[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & allele.ages2$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]])
    #ages.unfolded.theory[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & allele.ages2$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]])
    
    ages.lst[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$LST>quantiles.lst.perm[j] & allele.ages2$LST<quantiles.lst.perm[(j+1)]])
    ages.lst.perm[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$LST_PERM2>quantiles.lst.perm[j] & allele.ages2$LST_PERM2<quantiles.lst.perm[(j+1)]])
    
    ages.t[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$ABS_T>quantiles.t.perm[j] & allele.ages2$ABS_T<quantiles.t.perm[(j+1)]])
    ages.t.perm[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$ABS_T_PERM>quantiles.t.perm[j] & allele.ages2$ABS_T_PERM<quantiles.t.perm[(j+1)]])
    
    ages.unfolded.t[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$UNFOLDED_T>quantiles.unfolded.t.perm[j] & allele.ages2$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]])
    ages.unfolded.t.perm[j] <- mean(allele.ages2$AgeMedian_Jnt_definitive[ allele.ages2$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & allele.ages2$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]])
    
    print(j)

  }
 

## Alternative allele correction
quantiles.maf <- quantile(allele.ages2$AlleleAlt_Freq, prob = seq(0, 1, length = 21), type = 1)
quantiles.maf[1] <- 0.01
quantiles.maf[21] <- 0.99

ages.viability.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
ages.viability.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
#ages.viability.theory.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

ages.reproductive.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
ages.reproductive.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
#ages.reproductive.theory.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

ages.total.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
ages.total.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
#ages.total.theory.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

ages.unfolded.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
ages.unfolded.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
#ages.unfolded.theory.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

ages.lst.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
ages.lst.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

ages.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
ages.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

ages.unfolded.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
ages.unfolded.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

for (i in 1:20){
  for (j in 1:100){
    allele.ages3 <- subset(allele.ages2,AlleleAlt_Freq>quantiles.maf[i] & AlleleAlt_Freq<quantiles.maf[(i+1)])
    
    ages.viability.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$FST_VIABILITY_ST>quantiles.fst.theory[j] & allele.ages3$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]])
    ages.viability.perm.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & allele.ages3$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]])
    #ages.viability.theory.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$RANDOM_CHISQ>quantiles.fst.theory[j] & allele.ages3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]])
    
    ages.reproductive.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & allele.ages3$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]])
    ages.reproductive.perm.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & allele.ages3$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]])
    #ages.reproductive.theory.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$RANDOM_CHISQ>quantiles.fst.theory[j] & allele.ages3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]])
    
    ages.total.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$FST_GAMETIC_ST>quantiles.fst.theory[j] & allele.ages3$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]])
    ages.total.perm.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & allele.ages3$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]])
    #ages.total.theory.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$RANDOM_CHISQ>quantiles.fst.theory[j] & allele.ages3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]])
    
    ages.unfolded.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$UNFOLDED_FST>quantiles.unfolded.theory[j] & allele.ages3$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]])
    ages.unfolded.perm.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & allele.ages3$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]])
    #ages.unfolded.theory.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & allele.ages3$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]])
    
    ages.lst.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$LST>quantiles.lst.perm[j] & allele.ages3$LST<quantiles.lst.perm[(j+1)]])
    ages.lst.perm.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$LST_PERM2>quantiles.lst.perm[j] & allele.ages3$LST_PERM2<quantiles.lst.perm[(j+1)]])
    
    ages.t.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$ABS_T>quantiles.t.perm[j] & allele.ages3$ABS_T<quantiles.t.perm[(j+1)]])
    ages.t.perm.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$ABS_T_PERM>quantiles.t.perm[j] & allele.ages3$ABS_T_PERM<quantiles.t.perm[(j+1)]])
    
    ages.unfolded.t.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$UNFOLDED_T>quantiles.unfolded.t.perm[j] & allele.ages3$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]])
    ages.unfolded.t.perm.mafcorrected[i,j] <- mean(allele.ages3$AgeMedian_Jnt_definitive[ allele.ages3$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & allele.ages3$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]])
    
  }
  print(i)
}  


#Plots####
ages.viability.mafcorrected <- colMeans(ages.viability.mafcorrected[,1:100])
ages.viability.perm.mafcorrected <- colMeans(ages.viability.perm.mafcorrected[,1:100])
#ages.viability.theory.mafcorrected <- colMeans(ages.viability.theory.mafcorrected[,1:100])

ages.reproductive.mafcorrected <- colMeans(ages.reproductive.mafcorrected[,1:100])
ages.reproductive.perm.mafcorrected <- colMeans(ages.reproductive.perm.mafcorrected[,1:100])
#ages.reproductive.theory.mafcorrected <- colMeans(ages.reproductive.theory.mafcorrected[,1:100])

ages.total.mafcorrected <- colMeans(ages.total.mafcorrected[,1:100])
ages.total.perm.mafcorrected <- colMeans(ages.total.perm.mafcorrected[,1:100])
#ages.total.theory.mafcorrected <- colMeans(ages.total.theory.mafcorrected[,1:100])

ages.unfolded.mafcorrected <- colMeans(ages.unfolded.mafcorrected[,1:100])
ages.unfolded.perm.mafcorrected <- colMeans(ages.unfolded.perm.mafcorrected[,1:100])
#ages.unfolded.theory.mafcorrected <- colMeans(ages.unfolded.theory.mafcorrected[,1:100])


ages.lst.mafcorrected <- colMeans(ages.lst.mafcorrected[,1:100])
ages.lst.perm.mafcorrected <- colMeans(ages.lst.perm.mafcorrected[,1:100])

ages.t.mafcorrected <- colMeans(ages.t.mafcorrected[,1:100])
ages.t.perm.mafcorrected <- colMeans(ages.t.perm.mafcorrected[,1:100])

ages.unfolded.t.mafcorrected <- colMeans(ages.unfolded.t.mafcorrected[,1:100])
ages.unfolded.t.perm.mafcorrected <- colMeans(ages.unfolded.t.perm.mafcorrected[,1:100])


dd.quantiles <- data.frame(c(ages.viability,ages.viability.mafcorrected,
                             ages.viability.perm,ages.viability.perm.mafcorrected,
                             ages.reproductive,ages.reproductive.mafcorrected,
                             ages.reproductive.perm,ages.reproductive.perm.mafcorrected,
                             ages.total,ages.total.mafcorrected,
                             ages.total.perm,ages.total.perm.mafcorrected,
                             ages.unfolded,ages.unfolded.mafcorrected,
                             ages.unfolded.perm,ages.unfolded.perm.mafcorrected,
                             ages.lst,ages.lst.mafcorrected,
                             ages.lst.perm,ages.lst.perm.mafcorrected,
                             ages.t,ages.t.mafcorrected,
                             ages.t.perm,ages.t.perm.mafcorrected,
                             ages.unfolded.t,ages.unfolded.t.mafcorrected,
                             ages.unfolded.t.perm,ages.unfolded.t.perm.mafcorrected))
dd.quantiles$Type <- factor(rep(c(rep("Observed",200),rep("Permuted",200)),7))
dd.quantiles$Quantiles <- rep(c(1:100),28)
dd.quantiles$AltAF_corrected <- factor(rep(c(rep("No AltAF correction",100),rep("AltAF corrected",100)),14))
dd.quantiles$Metric <- factor(c(rep("Adult Fst",400),rep("Reproductive Fst",400),rep("Gametic Fst",400),rep("Unfolded Fst",400),rep("Lst",400),rep("|t|",400),rep("Unfolded t",400)))
names(dd.quantiles) <- c("Age","Type","Quantiles","AltAF_corrected","Metric")
dd.quantiles$AltAF_corrected <- factor(dd.quantiles$AltAF_corrected, levels = c("No AltAF correction", "AltAF corrected"))
dd.quantiles$Type <- factor(dd.quantiles$Type, levels = c("Permuted", "Observed"))
dd.quantiles <- dd.quantiles[order(dd.quantiles$Type,dd.quantiles$AltAF_corrected),]

ggplot(subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "darkorange","Adult Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "darkorange","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Adult Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="darkorange")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Adult Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Lst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="darkorange")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & AltAF_corrected=="AltAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("|t|" = "forestgreen","Reproductive Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "forestgreen","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Reproductive Fst" = 21,"|t|" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="forestgreen")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Reproductive Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="|t|") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="forestgreen")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & AltAF_corrected=="AltAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Unfolded t" = "#40B0A6","Unfolded Fst"="black"))+
  scale_fill_manual(values=c("Permuted" = "#40B0A6","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Unfolded Fst" = 21,"Unfolded t" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & AltAF_corrected=="AltAF corrected" & Type=="Permuted"), method='loess',formula=y~x,se=T,size=1.5,col="#40B0A6")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Permuted"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded t") & AltAF_corrected=="AltAF corrected" & Type=="Permuted"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="#40B0A6")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Gametic Fst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "purple","Gametic Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "purple","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Gametic Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="purple")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))


#Sup plots####
ggplot(subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "darkorange","Adult Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "darkorange","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Adult Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="darkorange")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Adult Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Lst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="darkorange")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & AltAF_corrected=="AltAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("|t|" = "forestgreen","Reproductive Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "forestgreen","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Reproductive Fst" = 21,"|t|" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="forestgreen")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Reproductive Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="|t|") & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="forestgreen")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & AltAF_corrected=="AltAF corrected" & Type=="Permuted"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Unfolded t" = "#40B0A6","Unfolded Fst"="black"))+
  scale_fill_manual(values=c("Permuted" = "#40B0A6","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Unfolded Fst" = 21,"Unfolded t" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & AltAF_corrected=="AltAF corrected" & Type=="Permuted"), method='loess',formula=y~x,se=T,size=1.5,col="#40B0A6")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Permuted"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded t") & AltAF_corrected=="AltAF corrected" & Type=="Permuted"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="#40B0A6")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Gametic Fst") & AltAF_corrected=="AltAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Age,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "purple","Gametic Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "purple","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Gametic Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="purple")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & AltAF_corrected=="AltAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  scale_y_continuous(limits=c(11900, 12900),)+
  ylab(expression(Mean~allele~age))+
  xlab(expression(Null~quantile))



#1,000 bootstraps#### 
diff.boots <- matrix(data=NA,ncol=10,nrow=10)
allele.ages2.s <- allele.ages2[,c("AgeMedian_Jnt_definitive","FST_VIABILITY_ST","FST_VIABILITY_PERM2_ST","FST_REPRODUCTIVE_ST","FST_REPRODUCTIVE_PERM_ST","FST_GAMETIC_ST","FST_GAMETIC_PERM2_ST","UNFOLDED_FST","UNFOLDED_FST_PERM","UNFOLDED_FST_THEORY","AlleleAlt_Freq")]

set.seed(123)
for (i in 1:10){
  boot <- sample(1:nrow(allele.ages2.s),replace=T)
  allele.ages2.s2 <- allele.ages2.s[boot,]
  
  diff.boots[i,1:10] <- c(cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_VIABILITY_ST,method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_VIABILITY_PERM2_ST,method = "spearman"),cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_VIABILITY_ST,method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$RANDOM_CHISQ,method = "spearman"),
                    cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_REPRODUCTIVE_ST,method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_REPRODUCTIVE_PERM_ST,method = "spearman"),cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_REPRODUCTIVE_ST,method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$RANDOM_CHISQ,method = "spearman"),
                         cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_GAMETIC_ST,method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_GAMETIC_PERM2_ST,method = "spearman"),cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$FST_GAMETIC_ST,method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive/allele.ages2.s2$AlleleAlt_Freq,allele.ages2.s2$RANDOM_CHISQ,method = "spearman"),
                    cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST>0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST>0],allele.ages2.s2$UNFOLDED_FST[allele.ages2.s2$UNFOLDED_FST>0],method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST_PERM>0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST_PERM>0],allele.ages2.s2$UNFOLDED_FST_PERM[allele.ages2.s2$UNFOLDED_FST_PERM>0],method = "spearman"),cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST>0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST>0],allele.ages2.s2$UNFOLDED_FST[allele.ages2.s2$UNFOLDED_FST>0],method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST_THEORY>0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST_THEORY>0],allele.ages2.s2$UNFOLDED_FST_THEORY[allele.ages2.s2$UNFOLDED_FST_THEORY>0],method = "spearman"),
                    cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST<0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST<0],allele.ages2.s2$UNFOLDED_FST[allele.ages2.s2$UNFOLDED_FST<0],method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST_PERM<0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST_PERM<0],allele.ages2.s2$UNFOLDED_FST_PERM[allele.ages2.s2$UNFOLDED_FST_PERM<0],method = "spearman"),cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST<0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST<0],allele.ages2.s2$UNFOLDED_FST[allele.ages2.s2$UNFOLDED_FST<0],method = "spearman")-cor(allele.ages2.s2$AgeMedian_Jnt_definitive[allele.ages2.s2$UNFOLDED_FST_THEORY<0]/allele.ages2.s2$AlleleAlt_Freq[allele.ages2.s2$UNFOLDED_FST_THEORY<0],allele.ages2.s2$UNFOLDED_FST_THEORY[allele.ages2.s2$UNFOLDED_FST_THEORY<0],method = "spearman"))
  print(i)
}

dd.boots.ages <- data.frame(melt(diff.boots))
names(dd.boots.ages) <- c("Replicate","Type","Difference")
dd.boots.ages$Type <- as.factor(dd.boots.ages$Type)
levels(dd.boots.ages$Type) <- c("Permuted_adult","Theory_adult","Permuted_reproductive","Theory_reproductive","Permuted_gametic","Theory_gametic","Permuted_unfoldedPos","Theory_unfoldedPos","Permuted_unfoldedNeg","Theory_unfoldedNeg")
dd.boots.ages$Type2 <- as.factor(do.call(rbind,strsplit(as.character(dd.boots.ages$Type),split="_"))[,2])
levels(dd.boots.ages$Type2) <- c("Adult","Gametic","Reproductive","UnfoldedPos","UnfoldedNeg")
levels(dd.boots.ages$Type) <- rep(c("Permuted","Theory"),5)

ggplot(subset(dd.boots.ages,Type2=="Adult"),aes(x=Difference,col=Type,fill=Type))+
  geom_histogram(alpha=0.5, position="identity")+
  theme_classic()+
  scale_fill_manual(values=c("darkorange","darkorange"))+
  scale_color_manual(values=c("darkgrey","black"))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20),strip.text=element_blank())+
  ylab(expression(Frequency))+
  xlab(expression(Observed~rho(hat(F)[ST],AltAge/AltAF)-Null~rho(hat(F)[ST],AltAge/AltAF)~(1000~boots.)))+
  ggtitle(expression(Adult~hat(F)[ST]))+
  xlim(c(min(subset(dd.boots.ages,Type2=="Adult")$Difference),max(subset(dd.boots.ages,Type2=="Adult")$Difference)))+
  geom_vline(xintercept=0,linetype="dashed")+
  facet_grid(Type~.,scale="free")
#p-value
sum(subset(dd.boots.ages,Type2=="Adult" &Type=="Permuted")$Difference<0)
#0.262
sum(subset(dd.boots.ages,Type2=="Adult" &Type=="Theory")$Difference<0)
#0.021

ggplot(subset(dd.boots.ages,Type2=="Reproductive"),aes(x=Difference,col=Type,fill=Type))+
  geom_histogram(alpha=0.5, position="identity")+
  theme_classic()+
  scale_fill_manual(values=c("forestgreen","forestgreen"))+
  scale_color_manual(values=c("darkgrey","black"))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20),strip.text=element_blank())+
  ylab(expression(Frequency))+
  xlab(expression(Observed~rho(hat(F)[ST],AltAge/AltAF)-Null~rho(hat(F)[ST],AltAge/AltAF)~(1000~boots.)))+
  ggtitle(expression(Reproductive~hat(F)[ST]))+
  xlim(c(min(subset(dd.boots.ages,Type2=="Reproductive")$Difference),max(subset(dd.boots.ages,Type2=="Reproductive")$Difference)))+
  geom_vline(xintercept=0,linetype="dashed")+
  facet_grid(Type~.,scale="free")
#p-value
sum(subset(dd.boots.ages,Type2=="Reproductive" &Type=="Permuted")$Difference<0)
#0.167
sum(subset(dd.boots.ages,Type2=="Reproductive" &Type=="Theory")$Difference<0)
#0.082

ggplot(subset(dd.boots.ages,Type2=="Gametic"),aes(x=Difference,col=Type,fill=Type))+
  geom_histogram(alpha=0.5, position="identity")+
  theme_classic()+
  scale_fill_manual(values=c("purple","purple"))+
  scale_color_manual(values=c("darkgrey","black"))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20),strip.text=element_blank())+
  ylab(expression(Frequency))+
  xlab(expression(Observed~rho(hat(F)[ST],AltAge/AltAF)-Null~rho(hat(F)[ST],AltAge/AltAF)~(1000~boots.)))+
  ggtitle(expression(Gametic~hat(F)[ST]))+
  xlim(c(min(subset(dd.boots.ages,Type2=="Gametic")$Difference),max(subset(dd.boots.ages,Type2=="Gametic")$Difference)))+
  geom_vline(xintercept=0,linetype="dashed")+
  facet_grid(Type~.,scale="free")
#p-value
sum(subset(dd.boots.ages,Type2=="Gametic" &Type=="Permuted")$Difference<0)
#0.097
sum(subset(dd.boots.ages,Type2=="Gametic" &Type=="Theory")$Difference<0)
#0.009


ggplot(subset(dd.boots.ages,Type2=="UnfoldedPos"),aes(x=Difference,col=Type,fill=Type))+
  geom_histogram(alpha=0.5, position="identity")+
  theme_classic()+
  scale_fill_manual(values=c("#40B0A6","#40B0A6"))+
  scale_color_manual(values=c("darkgrey","black"))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20),strip.text=element_blank())+
  ylab(expression(Frequency))+
  xlab(expression(Observed~rho(hat(F)[ST],AltAge/AltAF)-Null~rho(hat(F)[ST],AltAge/AltAF)~(1000~boots.)))+
  ggtitle(expression(Unfolded~reproductive~hat(F)[ST]))+
  xlim(c(min(subset(dd.boots.ages,Type2=="UnfoldedPos")$Difference),max(subset(dd.boots.ages,Type2=="UnfoldedPos")$Difference)))+
  geom_vline(xintercept=0,linetype="dashed")+
  facet_grid(Type~.,scale="free")
#p-value
sum(subset(dd.boots.ages,Type2=="Gametic" &Type=="Permuted")$Difference<0)
#0.097
sum(subset(dd.boots.ages,Type2=="Gametic" &Type=="Theory")$Difference<0)
#0.009


ggplot(subset(dd.boots.ages,Type2=="UnfoldedNeg"),aes(x=Difference,col=Type,fill=Type))+
  geom_histogram(alpha=0.5, position="identity")+
  theme_classic()+
  scale_fill_manual(values=c("#40B0A6","#40B0A6"))+
  scale_color_manual(values=c("darkgrey","black"))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20),strip.text=element_blank())+
  ylab(expression(Frequency))+
  xlab(expression(Observed~rho(hat(F)[ST],AltAge/AltAF)-Null~rho(hat(F)[ST],AltAge/AltAF)~(1000~boots.)))+
  ggtitle(expression(Unfolded~reproductive~hat(F)[ST]))+
  xlim(c(min(subset(dd.boots.ages,Type2=="UnfoldedNeg")$Difference),max(subset(dd.boots.ages,Type2=="UnfoldedNeg")$Difference)))+
  geom_vline(xintercept=0,linetype="dashed")+
  facet_grid(Type~.,scale="free")
#p-value
sum(subset(dd.boots.ages,Type2=="Gametic" &Type=="Permuted")$Difference<0)
#0.097
sum(subset(dd.boots.ages,Type2=="Gametic" &Type=="Theory")$Difference<0)
#0.009

#...####
#...####
#...####

#1000G, Fst####

#Prep####
dir <- "~/Dropbox/dropbox_work/data/biobank/1000_genomes_data/"

## Select relevant samples from each sub-population

vcf_sampleids <- read.table(paste0(dir,"sample_idsA.txt"))

sample_info <- read.table(paste0(dir,"igsr_samples_edited.tsv"),h=T,sep="\t")
sample_info <- subset(sample_info,Sample.name %in% vcf_sampleids$V1)

for (i in levels(sample_info$Population.code)[2:30]){
write.table(subset(sample_info,Population.code==i)$Sample.name,paste0(dir,"sample_ids_",i,".txt"),quote=F,col.names=F,row.names=F)
}

##Allele frequencies among subpopulations
pop_ids <- levels(sample_info$Population.code)[c(2:7,9:12,14:26,28)]
write.table(pop_ids,paste0(dir,"pop_ids"),quote=F,col.names=F,row.names=F)


#Import####

dir <- "~/Dropbox/dropbox_work/data/biobank/"

#Genotyped
tmp2 <- list()
for (i in c("CHS","PUR","TSI","GIH","YRI")){
  tmp <- read.table(paste0(dir,"1000_genomes_data/ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v2.",i))
  names(tmp) <- c("CHROM","POS","ID","REF","ALT","AC","AN")
  tmp$AC <- as.numeric(as.character(tmp$AC))
  tmp$MAF <- ifelse(tmp$AC/tmp$AN<0.5,tmp$AC/tmp$AN,1-tmp$AC/tmp$AN)
  tmp2[[i]] <- tmp
}
all.frq <- do.call(cbind,tmp2)
all.frq <- all.frq[,c(1:5,8,16,24,32,40)]
names(all.frq)[1:5] <- c("CHROM","POS","ID","REF2","ALT2")
rm(tmp2)
rm(tmp)

#Imputed
tmp2 <- list()
for (i in c("CHS","PUR","TSI","GIH","YRI")){
  tmp <- read.table(paste0(dir,"1000_genomes_data/ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.",i))
  names(tmp) <- c("CHROM","POS","ID","REF","ALT","AC","AN")
  tmp$AC <- as.numeric(as.character(tmp$AC))
  tmp$MAF <- ifelse(tmp$AC/tmp$AN<0.5,tmp$AC/tmp$AN,1-tmp$AC/tmp$AN)
  tmp2[[i]] <- tmp
}
all.frq <- do.call(cbind,tmp2)
all.frq <- all.frq[,c(1:5,8,16,24,32,40)]
names(all.frq)[1:5] <- c("CHROM","POS","ID","REF2","ALT2")
rm(tmp2)
rm(tmp)


#Merge estimates of between-pop Fst with estimates of between-sex Fst
mf.frqs2 <- merge(mf.frqs,all.frq,by="ID",all.x=T)
mf.gwases2 <- merge(mf.gwases,all.frq,by="ID",all.x=T)
mf.gwases.lrs2 <- merge(mf.gwases.lrs,all.frq,by="ID",all.x=T)
rm(all.frq)

#Remove sites which do not match ALT and REF, and CHROM
mf.frqs2 <- subset(mf.frqs2,(REF2==REF.y | REF2==ALT.y) & (ALT2==REF.y | ALT2==ALT.y) & CHROM.x==CHROM.y)
mf.gwases2 <- subset(mf.gwases2,(REF2==REF.y | REF2==ALT.y) & (ALT2==REF.y | ALT2==ALT.y) & CHROM.x==CHROM.y)
mf.gwases.lrs2 <- subset(mf.gwases.lrs2,(REF2==REF.y | REF2==ALT.y) & (ALT2==REF.y | ALT2==ALT.y) & CHROM.x==CHROM.y)

#all.frq2$FST_CHS_PUR <- with(all.frq2,ifelse(!is.na(CHS.MAF) & !is.na(PUR.MAF),nei_fst_func(CHS.MAF,PUR.MAF),NA))
#all.frq2$FST_CHS_GIH <- with(all.frq2,ifelse(!is.na(CHS.MAF) & !is.na(GIH.MAF),nei_fst_func(CHS.MAF,GIH.MAF),NA))
#all.frq2$FST_CHS_YRI <- with(all.frq2,ifelse(!is.na(CHS.MAF) & !is.na(YRI.MAF),nei_fst_func(CHS.MAF,YRI.MAF),NA))
#all.frq2$FST_PUR_GIH <- with(all.frq2,ifelse(!is.na(PUR.MAF) & !is.na(GIH.MAF),nei_fst_func(PUR.MAF,GIH.MAF),NA))
#all.frq2$FST_PUR_YRI <- with(all.frq2,ifelse(!is.na(PUR.MAF) & !is.na(YRI.MAF),nei_fst_func(PUR.MAF,YRI.MAF),NA))
mf.frqs2$FST_GIH_YRI <- with(mf.frqs2,ifelse(!is.na(GIH.MAF) & !is.na(YRI.MAF),nei_fst_func(GIH.MAF,YRI.MAF),NA))
mf.gwases2$FST_GIH_YRI <- with(mf.gwases2,ifelse(!is.na(GIH.MAF) & !is.na(YRI.MAF),nei_fst_func(GIH.MAF,YRI.MAF),NA))
mf.gwases.lrs2$FST_GIH_YRI <- with(mf.gwases.lrs2,ifelse(!is.na(GIH.MAF) & !is.na(YRI.MAF),nei_fst_func(GIH.MAF,YRI.MAF),NA))



#Statistical tests####

#Viability Fst
summary(glm(data=subset(mf.frqs2,is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~FST_VIABILITY_ST+MAF))
#-0.0002049, p=0.879

#Reproductive Fst
summary(glm(data=subset(mf.frqs2,is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~FST_REPRODUCTIVE_ST+MAF))
#0.0004504, p=0.743

#Gametic Fst
summary(glm(data=subset(mf.frqs2,is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~FST_GAMETIC_ST+MAF))
#0.001650, p=0.227

#Unfolded Fst
summary(glm(data=subset(mf.frqs2,UNFOLDED_FST<0 & is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~UNFOLDED_FST+MAF))
#-0.007285, p=0.0419
summary(glm(data=subset(mf.frqs2,UNFOLDED_FST>0 & is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~UNFOLDED_FST+MAF))
#0.006853, p=0.030

#Lst
summary(glm(data=subset(mf.gwases2,is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~LST+MAF))
#-2.373e+02, p=0.000273

#|t|
summary(glm(data=subset(mf.gwases.lrs2,is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~ABS_T+MAF))
#3.919e-05, p=0.629

#Unfolded t
summary(glm(data=subset(mf.gwases.lrs2,UNFOLDED_T<0 & is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~UNFOLDED_T+MAF))
#-0.009006, p=0.012
summary(glm(data=subset(mf.gwases.lrs2,UNFOLDED_T>0 & is.finite(log(FST_GIH_YRI))),log(FST_GIH_YRI)~UNFOLDED_T+MAF))
#0.006571, p=0.0363


#Quantile-proportion plot####

#No MAF correction
quantiles.fst.theory <- quantile(mf.frqs2$RANDOM_CHISQ, prob = seq(0, 1, length = 101), type = 1)
quantiles.fst.theory[1] <- 0
quantiles.fst.theory[101] <- 10^6

quantiles.unfolded.theory <- quantile(mf.frqs2$UNFOLDED_FST_THEORY, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.theory[1] <- -10^6
quantiles.unfolded.theory[101] <- 10^6

quantiles.lst.perm <- quantile(mf.gwases2$LST_PERM2, prob = seq(0, 1, length = 101), type = 1)
quantiles.lst.perm[1] <- 0
quantiles.lst.perm[101] <- 10^6

quantiles.t.perm <- quantile(mf.gwases.lrs2$ABS_T_PERM, prob = seq(0, 1, length = 101), type = 1)
quantiles.t.perm[1] <- 0
quantiles.t.perm[101] <- 10^6

quantiles.unfolded.t.perm <- quantile(mf.gwases.lrs2$UNFOLDED_T, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.t.perm[1] <- -10^6
quantiles.unfolded.t.perm[101] <- 10^6


btwn.viability <- vector("numeric")
btwn.viability.perm <- vector("numeric")
btwn.reproductive <- vector("numeric")
btwn.reproductive.perm <- vector("numeric")
btwn.total <- vector("numeric")
btwn.total.perm <- vector("numeric")
btwn.unfolded <- vector("numeric")
btwn.unfolded.perm <- vector("numeric")
btwn.lst <- vector()
btwn.lst.perm <- vector()
btwn.t <- vector()
btwn.t.perm <- vector()
btwn.unfolded.t <- vector("numeric")
btwn.unfolded.t.perm <- vector("numeric")

for (j in 1:100){
  
  btwn.viability[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs2$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  btwn.viability.perm[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs2$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #btwn.viability.theory[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs2$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  btwn.reproductive[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs2$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  btwn.reproductive.perm[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs2$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #btwn.reproductive.theory[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs2$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  btwn.total[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs2$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  btwn.total.perm[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs2$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #btwn.total.theory[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs2$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  btwn.unfolded[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs2$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  btwn.unfolded.perm[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs2$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  #btwn.unfolded.theory[j] <- mean(mf.frqs2$FST_GIH_YRI[ mf.frqs2$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs2$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  
  btwn.lst[j] <- mean(mf.gwases2$FST_GIH_YRI[ mf.gwases2$LST>quantiles.lst.perm[j] & mf.gwases2$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
  btwn.lst.perm[j] <- mean(mf.gwases2$FST_GIH_YRI[ mf.gwases2$LST_PERM2>quantiles.lst.perm[j] & mf.gwases2$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
  
  btwn.t[j] <- mean(mf.gwases.lrs2$FST_GIH_YRI[ mf.gwases.lrs2$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs2$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
  btwn.t.perm[j] <- mean(mf.gwases.lrs2$FST_GIH_YRI[ mf.gwases.lrs2$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs2$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
  
  btwn.unfolded.t[j] <- mean(mf.gwases.lrs2$FST_GIH_YRI[ mf.gwases.lrs2$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs2$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  btwn.unfolded.t.perm[j] <- mean(mf.gwases.lrs2$FST_GIH_YRI[ mf.gwases.lrs2$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs2$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  
  
  print(j)
  
}


## MAF correction
quantiles.maf <- quantile(mf.frqs2$MAF, prob = seq(0, 1, length = 21), type = 1)
quantiles.maf[1] <- 0.01
quantiles.maf[21] <- 0.5

btwn.viability.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
btwn.viability.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

btwn.reproductive.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
btwn.reproductive.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

btwn.total.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
btwn.total.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

btwn.unfolded.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
btwn.unfolded.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

btwn.lst.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
btwn.lst.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

btwn.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
btwn.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

btwn.unfolded.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
btwn.unfolded.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)


for (i in 1:20){
  for (j in 1:100){
    mf.frqs3 <- subset(mf.frqs2,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases3 <- subset(mf.gwases2,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases.lrs3 <- subset(mf.gwases.lrs2,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    
    btwn.viability.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    btwn.viability.perm.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #btwn.viability.theory[j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    btwn.reproductive.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    btwn.reproductive.perm.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #btwn.reproductive.theory[j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    btwn.total.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    btwn.total.perm.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #btwn.total.theory[j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    btwn.unfolded.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    btwn.unfolded.perm.mafcorrected[i,j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    #btwn.unfolded.theory[j] <- mean(mf.frqs3$FST_GIH_YRI[ mf.frqs3$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    
    btwn.lst.mafcorrected[i,j] <- mean(mf.gwases3$FST_GIH_YRI[ mf.gwases3$LST>quantiles.lst.perm[j] & mf.gwases3$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
    btwn.lst.perm.mafcorrected[i,j] <- mean(mf.gwases3$FST_GIH_YRI[ mf.gwases3$LST_PERM2>quantiles.lst.perm[j] & mf.gwases3$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
    
    btwn.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$FST_GIH_YRI[ mf.gwases.lrs3$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
    btwn.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$FST_GIH_YRI[ mf.gwases.lrs3$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
    
    btwn.unfolded.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$FST_GIH_YRI[ mf.gwases.lrs3$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
    btwn.unfolded.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$FST_GIH_YRI[ mf.gwases.lrs3$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  }
  print(i)
}  


#Plots####
btwn.viability.mafcorrected <- colMeans(btwn.viability.mafcorrected[,1:100])
btwn.viability.perm.mafcorrected <- colMeans(btwn.viability.perm.mafcorrected[,1:100])
#btwn.viability.theory.mafcorrected <- colMeans(btwn.viability.theory.mafcorrected[,1:100])

btwn.reproductive.mafcorrected <- colMeans(btwn.reproductive.mafcorrected[,1:100])
btwn.reproductive.perm.mafcorrected <- colMeans(btwn.reproductive.perm.mafcorrected[,1:100])
#btwn.reproductive.theory.mafcorrected <- colMeans(btwn.reproductive.theory.mafcorrected[,1:100])

btwn.total.mafcorrected <- colMeans(btwn.total.mafcorrected[,1:100])
btwn.total.perm.mafcorrected <- colMeans(btwn.total.perm.mafcorrected[,1:100])
#btwn.total.theory.mafcorrected <- colMeans(btwn.total.theory.mafcorrected[,1:100])

btwn.unfolded.mafcorrected <- colMeans(btwn.unfolded.mafcorrected[,1:100])
btwn.unfolded.perm.mafcorrected <- colMeans(btwn.unfolded.perm.mafcorrected[,1:100])
#btwn.unfolded.theory.mafcorrected <- colMeans(btwn.unfolded.theory.mafcorrected[,1:100])

btwn.lst.mafcorrected <- colMeans(btwn.lst.mafcorrected[,1:100])
btwn.lst.perm.mafcorrected <- colMeans(btwn.lst.perm.mafcorrected[,1:100])

btwn.t.mafcorrected <- colMeans(btwn.t.mafcorrected[,1:100])
btwn.t.perm.mafcorrected <- colMeans(btwn.t.perm.mafcorrected[,1:100])

btwn.unfolded.t.mafcorrected <- colMeans(btwn.unfolded.t.mafcorrected[,1:100])
btwn.unfolded.t.perm.mafcorrected <- colMeans(btwn.unfolded.t.perm.mafcorrected[,1:100])


dd.quantiles <- data.frame(c(btwn.viability,btwn.viability.mafcorrected,
                             btwn.viability.perm,btwn.viability.perm.mafcorrected,
                             btwn.reproductive,btwn.reproductive.mafcorrected,
                             btwn.reproductive.perm,btwn.reproductive.perm.mafcorrected,
                             btwn.total,btwn.total.mafcorrected,
                             btwn.total.perm,btwn.total.perm.mafcorrected,
                             btwn.unfolded,btwn.unfolded.mafcorrected,
                             btwn.unfolded.perm,btwn.unfolded.perm.mafcorrected,
                             btwn.lst,btwn.lst.mafcorrected,
                             btwn.lst.perm,btwn.lst.perm.mafcorrected,
                             btwn.t,btwn.t.mafcorrected,
                             btwn.t.perm,btwn.t.perm.mafcorrected,
                             btwn.unfolded.t,btwn.unfolded.t.mafcorrected,
                             btwn.unfolded.t.perm,btwn.unfolded.t.perm.mafcorrected))
dd.quantiles$Type <- factor(rep(c(rep("Observed",200),rep("Permuted",200)),7))
dd.quantiles$Quantiles <- rep(c(1:100),28)
dd.quantiles$MAF_corrected <- factor(rep(c(rep("No MAF correction",100),rep("MAF corrected",100)),14))
dd.quantiles$Metric <- factor(c(rep("Adult Fst",400),rep("Reproductive Fst",400),rep("Gametic Fst",400),rep("Unfolded Fst",400),rep("Lst",400),rep("|t|",400),rep("Unfolded t",400)))
names(dd.quantiles) <- c("Proportion","Type","Quantiles","MAF_corrected","Metric")
dd.quantiles$MAF_corrected <- factor(dd.quantiles$MAF_corrected, levels = c("No MAF correction", "MAF corrected"))
dd.quantiles$Type <- factor(dd.quantiles$Type, levels = c("Permuted", "Observed"))
dd.quantiles <- dd.quantiles[order(dd.quantiles$Type,dd.quantiles$MAF_corrected),]

ggplot(subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "darkorange","Adult Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "darkorange","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Adult Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="darkorange")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Adult Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="darkorange")+
  scale_y_continuous(limits=c(0.024, 0.029),)+
  ylab(expression(hat(F)[ST[GIH-YRI]]))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("|t|" = "forestgreen","Reproductive Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "forestgreen","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Reproductive Fst" = 21,"|t|" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="forestgreen")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Reproductive Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="forestgreen")+
  scale_y_continuous(limits=c(0.024, 0.029),)+
  ylab(expression(hat(F)[ST[GIH-YRI]]))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Unfolded t" = "#40B0A6","Unfolded Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "#40B0A6","Observed"="darkgrey"))+
  scale_shape_manual(values=c("Unfolded Fst" = 21,"Unfolded t" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Unfolded Fst" | Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="#40B0A6")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="#40B0A6")+
  scale_y_continuous(limits=c(0.024, 0.029),)+
  ylab(expression(hat(F)[ST[GIH-YRI]]))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Gametic Fst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "purple","Gametic Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "purple","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Gametic Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="purple")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  scale_y_continuous(limits=c(0.024, 0.029),)+
  ylab(expression(hat(F)[ST[GIH-YRI]]))+
  xlab(expression(Null~quantile))



#...####
#...####
#...####

#BALANCING SELECTION CANDIDATES####
dir <- "~/Dropbox/dropbox_work/data/biobank/"

#Import####

#Make mf.frqs table smaller, if needed
#mf.frqs <- mf.frqs[,c("CHROM","ID","REF.x","ALT.x","FST_VIABILITY","MAF","FST","FST_PERM2","FST_VIABILITY_PERM","FST_THEORY","FST_ADULT_THEORY_APPROX","FST_ADULT_PERM","FST_ADULT","POS","REF.y","ALT.y","EFFECT","GENE")]


bitarello.genes <- fread(paste0(dir,"balancing_selection_data/Bitarello_et_al_2018_Candidate_windows/bitarello_genes_cleaned.csv"))
andres.genes <- fread(paste0(dir,"balancing_selection_data/Andres_et_al_2009_Candidate_genes/andres_genes_cleaned.csv"))
degiorgio.genes <- fread(paste0(dir,"balancing_selection_data/DeGiorgio_et_al_2014_Candidate_genes/degiorgio_genes_cleaned.csv"))

#Columns that identify balanced genes from each dataset
#Bitarello et al. 2018
mf.frqs$Bitarello_candidate <- ifelse(mf.frqs$GENE %in% bitarello.genes$Acronym,1,0)
mf.gwases$Bitarello_candidate <- ifelse(mf.gwases$GENE %in% bitarello.genes$Acronym,1,0)
mf.gwases.lrs$Bitarello_candidate <- ifelse(mf.gwases.lrs$GENE %in% bitarello.genes$Acronym,1,0)
#Andres et al 2009
mf.frqs$Andres_candidate <- ifelse(mf.frqs$GENE %in% andres.genes$Acronym,1,0)
mf.gwases$Andres_candidate <- ifelse(mf.gwases$GENE %in% andres.genes$Acronym,1,0)
mf.gwases.lrs$Andres_candidate <- ifelse(mf.gwases.lrs$GENE %in% andres.genes$Acronym,1,0)
#DeGiorgio et al 2014
mf.frqs$DeGiorgio_candidate <- ifelse(mf.frqs$GENE %in% degiorgio.genes$Acronym,1,0)
mf.gwases$DeGiorgio_candidate <- ifelse(mf.gwases$GENE %in% degiorgio.genes$Acronym,1,0)
mf.gwases.lrs$DeGiorgio_candidate <- ifelse(mf.gwases.lrs$GENE %in% degiorgio.genes$Acronym,1,0)


#...####

#Bitarello, Statistical tests####

#Viability Fst
summary(glm(data=mf.frqs,Bitarello_candidate~FST_VIABILITY_ST+MAF,family="binomial"))
#0.000853, p=0.662

#Reproductive Fst
summary(glm(data=mf.frqs,Bitarello_candidate~FST_REPRODUCTIVE_ST+MAF,family="binomial"))
#0.004464, p=0.0236

#Gametic Fst
summary(glm(data=mf.frqs,Bitarello_candidate~FST_GAMETIC_ST+MAF,family="binomial"))
#0.0007403, p=0.708

#Unfolded Fst
summary(glm(data=subset(mf.frqs,UNFOLDED_FST<0),Bitarello_candidate~UNFOLDED_FST+MAF,family="binomial"))
#0.001542, p=0.766
summary(glm(data=subset(mf.frqs,UNFOLDED_FST>0),Bitarello_candidate~UNFOLDED_FST+MAF,family="binomial"))
#0.011127, p=0.0146

#Lst
summary(glm(data=mf.gwases,Bitarello_candidate~LST+MAF,family="binomial"))
#-5.995e+02, p=0.873

#|t|
summary(glm(data=mf.gwases.lrs,Bitarello_candidate~ABS_T+MAF,family="binomial"))
#0.008170, p=0.0801

#Unfolded t
summary(glm(data=subset(mf.gwases.lrs,UNFOLDED_T<0),Bitarello_candidate~UNFOLDED_T+MAF,family="binomial"))
#0.0005491, p=0.916
summary(glm(data=subset(mf.gwases.lrs,UNFOLDED_T>0),Bitarello_candidate~UNFOLDED_T+MAF,family="binomial"))
#0.011537, p=0.0104



#Quantile-proportion plot####

#No MAF correction
quantiles.fst.theory <- quantile(mf.frqs$RANDOM_CHISQ, prob = seq(0, 1, length = 101), type = 1)
quantiles.fst.theory[1] <- 0
quantiles.fst.theory[101] <- 10^6

quantiles.unfolded.theory <- quantile(mf.frqs$UNFOLDED_FST_THEORY, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.theory[1] <- -10^6
quantiles.unfolded.theory[101] <- 10^6

quantiles.lst.perm <- quantile(mf.gwases$LST_PERM2, prob = seq(0, 1, length = 101), type = 1)
quantiles.lst.perm[1] <- 0
quantiles.lst.perm[101] <- 10^6

quantiles.t.perm <- quantile(mf.gwases.lrs$ABS_T_PERM, prob = seq(0, 1, length = 101), type = 1)
quantiles.t.perm[1] <- 0
quantiles.t.perm[101] <- 10^6

quantiles.unfolded.t.perm <- quantile(mf.gwases.lrs$UNFOLDED_T, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.t.perm[1] <- -10^6
quantiles.unfolded.t.perm[101] <- 10^6


bal.viability <- vector("numeric")
bal.viability.perm <- vector("numeric")
bal.reproductive <- vector("numeric")
bal.reproductive.perm <- vector("numeric")
bal.total <- vector("numeric")
bal.total.perm <- vector("numeric")
bal.unfolded <- vector("numeric")
bal.unfolded.perm <- vector("numeric")
bal.lst <- vector()
bal.lst.perm <- vector()
bal.t <- vector()
bal.t.perm <- vector()
bal.unfolded.t <- vector("numeric")
bal.unfolded.t.perm <- vector("numeric")

for (j in 1:100){
  
  bal.viability[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.viability.perm[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.viability.theory[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.reproductive[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.reproductive.perm[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.reproductive.theory[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.total[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.total.perm[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.total.theory[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.unfolded[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  bal.unfolded.perm[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  #bal.unfolded.theory[j] <- mean(mf.frqs$Bitarello_candidate[ mf.frqs$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  
  bal.lst[j] <- mean(mf.gwases$Bitarello_candidate[ mf.gwases$LST>quantiles.lst.perm[j] & mf.gwases$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
  bal.lst.perm[j] <- mean(mf.gwases$Bitarello_candidate[ mf.gwases$LST_PERM2>quantiles.lst.perm[j] & mf.gwases$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
  
  bal.t[j] <- mean(mf.gwases.lrs$Bitarello_candidate[ mf.gwases.lrs$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
  bal.t.perm[j] <- mean(mf.gwases.lrs$Bitarello_candidate[ mf.gwases.lrs$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
  
  bal.unfolded.t[j] <- mean(mf.gwases.lrs$Bitarello_candidate[ mf.gwases.lrs$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  bal.unfolded.t.perm[j] <- mean(mf.gwases.lrs$Bitarello_candidate[ mf.gwases.lrs$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  
  
  print(j)
  
}


## MAF correction
quantiles.maf <- quantile(mf.frqs$MAF, prob = seq(0, 1, length = 21), type = 1)
quantiles.maf[1] <- 0.01
quantiles.maf[21] <- 0.5

bal.viability.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.viability.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.reproductive.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.reproductive.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.total.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.total.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.unfolded.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.unfolded.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.lst.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.lst.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.unfolded.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.unfolded.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)


for (i in 1:20){
  for (j in 1:100){
    mf.frqs3 <- subset(mf.frqs,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases3 <- subset(mf.gwases,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases.lrs3 <- subset(mf.gwases.lrs,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    
    bal.viability.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.viability.perm.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.viability.theory[j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.reproductive.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.reproductive.perm.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.reproductive.theory[j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.total.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.total.perm.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.total.theory[j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.unfolded.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    bal.unfolded.perm.mafcorrected[i,j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    #bal.unfolded.theory[j] <- mean(mf.frqs3$Bitarello_candidate[ mf.frqs3$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    
    bal.lst.mafcorrected[i,j] <- mean(mf.gwases3$Bitarello_candidate[ mf.gwases3$LST>quantiles.lst.perm[j] & mf.gwases3$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
    bal.lst.perm.mafcorrected[i,j] <- mean(mf.gwases3$Bitarello_candidate[ mf.gwases3$LST_PERM2>quantiles.lst.perm[j] & mf.gwases3$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
    
    bal.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Bitarello_candidate[ mf.gwases.lrs3$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
    bal.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Bitarello_candidate[ mf.gwases.lrs3$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
    
    bal.unfolded.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Bitarello_candidate[ mf.gwases.lrs3$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
    bal.unfolded.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Bitarello_candidate[ mf.gwases.lrs3$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  }
  print(i)
}  


#Plots####
bal.viability.mafcorrected <- colMeans(bal.viability.mafcorrected[,1:100])
bal.viability.perm.mafcorrected <- colMeans(bal.viability.perm.mafcorrected[,1:100])
#bal.viability.theory.mafcorrected <- colMeans(bal.viability.theory.mafcorrected[,1:100])

bal.reproductive.mafcorrected <- colMeans(bal.reproductive.mafcorrected[,1:100])
bal.reproductive.perm.mafcorrected <- colMeans(bal.reproductive.perm.mafcorrected[,1:100])
#bal.reproductive.theory.mafcorrected <- colMeans(bal.reproductive.theory.mafcorrected[,1:100])

bal.total.mafcorrected <- colMeans(bal.total.mafcorrected[,1:100])
bal.total.perm.mafcorrected <- colMeans(bal.total.perm.mafcorrected[,1:100])
#bal.total.theory.mafcorrected <- colMeans(bal.total.theory.mafcorrected[,1:100])

bal.unfolded.mafcorrected <- colMeans(bal.unfolded.mafcorrected[,1:100])
bal.unfolded.perm.mafcorrected <- colMeans(bal.unfolded.perm.mafcorrected[,1:100])
#bal.unfolded.theory.mafcorrected <- colMeans(bal.unfolded.theory.mafcorrected[,1:100])

bal.lst.mafcorrected <- colMeans(bal.lst.mafcorrected[,1:100])
bal.lst.perm.mafcorrected <- colMeans(bal.lst.perm.mafcorrected[,1:100])

bal.t.mafcorrected <- colMeans(bal.t.mafcorrected[,1:100])
bal.t.perm.mafcorrected <- colMeans(bal.t.perm.mafcorrected[,1:100])

bal.unfolded.t.mafcorrected <- colMeans(bal.unfolded.t.mafcorrected[,1:100])
bal.unfolded.t.perm.mafcorrected <- colMeans(bal.unfolded.t.perm.mafcorrected[,1:100])


dd.quantiles <- data.frame(c(bal.viability,bal.viability.mafcorrected,
                             bal.viability.perm,bal.viability.perm.mafcorrected,
                             bal.reproductive,bal.reproductive.mafcorrected,
                             bal.reproductive.perm,bal.reproductive.perm.mafcorrected,
                             bal.total,bal.total.mafcorrected,
                             bal.total.perm,bal.total.perm.mafcorrected,
                             bal.unfolded,bal.unfolded.mafcorrected,
                             bal.unfolded.perm,bal.unfolded.perm.mafcorrected,
                             bal.lst,bal.lst.mafcorrected,
                             bal.lst.perm,bal.lst.perm.mafcorrected,
                             bal.t,bal.t.mafcorrected,
                             bal.t.perm,bal.t.perm.mafcorrected,
                             bal.unfolded.t,bal.unfolded.t.mafcorrected,
                             bal.unfolded.t.perm,bal.unfolded.t.perm.mafcorrected))
dd.quantiles$Type <- factor(rep(c(rep("Observed",200),rep("Permuted",200)),7))
dd.quantiles$Quantiles <- rep(c(1:100),28)
dd.quantiles$MAF_corrected <- factor(rep(c(rep("No MAF correction",100),rep("MAF corrected",100)),14))
dd.quantiles$Metric <- factor(c(rep("Adult Fst",400),rep("Reproductive Fst",400),rep("Gametic Fst",400),rep("Unfolded Fst",400),rep("Lst",400),rep("|t|",400),rep("Unfolded t",400)))
names(dd.quantiles) <- c("Proportion","Type","Quantiles","MAF_corrected","Metric")
dd.quantiles$MAF_corrected <- factor(dd.quantiles$MAF_corrected, levels = c("No MAF correction", "MAF corrected"))
dd.quantiles$Type <- factor(dd.quantiles$Type, levels = c("Permuted", "Observed"))
dd.quantiles <- dd.quantiles[order(dd.quantiles$Type,dd.quantiles$MAF_corrected),]

ggplot(subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "darkorange","Adult Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "darkorange","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Adult Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="darkorange")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Adult Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="darkorange")+
  scale_y_continuous(limits=c(0.12, 0.145),)+
  ylab(expression(Proportion~candidates~(Bitarello~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("|t|" = "forestgreen","Reproductive Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "forestgreen","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Reproductive Fst" = 21,"|t|" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="forestgreen")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Reproductive Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="forestgreen")+
  scale_y_continuous(limits=c(0.12, 0.145),)+
  ylab(expression(Proportion~candidates~(Bitarello~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Unfolded t" = "#40B0A6","Unfolded Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "#40B0A6","Observed"="darkgrey"))+
  scale_shape_manual(values=c("Unfolded Fst" = 21,"Unfolded t" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="#40B0A6")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="#40B0A6")+
  scale_y_continuous(limits=c(0.12, 0.145),)+
  ylab(expression(Proportion~candidates~(Bitarello~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Gametic Fst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "purple","Gametic Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "purple","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Gametic Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="purple")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  scale_y_continuous(limits=c(0.12, 0.145),)+
  ylab(expression(Proportion~candidates~(Bitarello~et~al.)))+
  xlab(expression(Null~quantile))




#...####

##Andres, Statistical tests####

#Viability Fst
summary(glm(data=mf.frqs,Andres_candidate~FST_VIABILITY_ST+MAF,family="binomial"))
#-0.01577, p=0.391

#Reproductive Fst
summary(glm(data=mf.frqs,Andres_candidate~FST_REPRODUCTIVE_ST+MAF,family="binomial"))
#-0.01927, p=0.306

#Gametic Fst
summary(glm(data=mf.frqs,Andres_candidate~FST_GAMETIC_ST+MAF,family="binomial"))
#-0.008314, p=0.651

#Unfolded Fst
summary(glm(data=subset(mf.frqs,UNFOLDED_FST<0),Andres_candidate~UNFOLDED_FST+MAF,family="binomial"))
#-0.003074, p=0.9496
summary(glm(data=subset(mf.frqs,UNFOLDED_FST>0),Andres_candidate~UNFOLDED_FST+MAF,family="binomial"))
#0.01481, p=0.706

#Lst
summary(glm(data=mf.gwases,Andres_candidate~LST+MAF,family="binomial"))
#-3.508e+04, p=0.326

#|t|
summary(glm(data=mf.gwases.lrs,Andres_candidate~ABS_T+MAF,family="binomial"))
#-0.02263, p=0.599

#Unfolded t
summary(glm(data=subset(mf.gwases.lrs,UNFOLDED_T<0),Andres_candidate~UNFOLDED_T+MAF,family="binomial"))
#0.01457, p=0.7656
summary(glm(data=subset(mf.gwases.lrs,UNFOLDED_T>0),Andres_candidate~UNFOLDED_T+MAF,family="binomial"))
#0.04055, p=0.252

#Quantile-proportion plot####

#No MAF correction
quantiles.fst.theory <- quantile(mf.frqs$RANDOM_CHISQ, prob = seq(0, 1, length = 101), type = 1)
quantiles.fst.theory[1] <- 0
quantiles.fst.theory[101] <- 10^6

quantiles.unfolded.theory <- quantile(mf.frqs$UNFOLDED_FST_THEORY, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.theory[1] <- -10^6
quantiles.unfolded.theory[101] <- 10^6

quantiles.lst.perm <- quantile(allele.ages2$LST_PERM2, prob = seq(0, 1, length = 101), type = 1)
quantiles.lst.perm[1] <- 0
quantiles.lst.perm[101] <- 10^6

quantiles.t.perm <- quantile(allele.ages2$ABS_T_PERM, prob = seq(0, 1, length = 101), type = 1)
quantiles.t.perm[1] <- 0
quantiles.t.perm[101] <- 10^6

quantiles.unfolded.t.perm <- quantile(allele.ages2$UNFOLDED_T, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.t.perm[1] <- -10^6
quantiles.unfolded.t.perm[101] <- 10^6


bal.viability <- vector("numeric")
bal.viability.perm <- vector("numeric")
bal.reproductive <- vector("numeric")
bal.reproductive.perm <- vector("numeric")
bal.total <- vector("numeric")
bal.total.perm <- vector("numeric")
bal.unfolded <- vector("numeric")
bal.unfolded.perm <- vector("numeric")
bal.lst <- vector()
bal.lst.perm <- vector()
bal.t <- vector()
bal.t.perm <- vector()
bal.unfolded.t <- vector("numeric")
bal.unfolded.t.perm <- vector("numeric")

for (j in 1:100){
  
  bal.viability[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.viability.perm[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.viability.theory[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.reproductive[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.reproductive.perm[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.reproductive.theory[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.total[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.total.perm[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.total.theory[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.unfolded[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  bal.unfolded.perm[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  #bal.unfolded.theory[j] <- mean(mf.frqs$Andres_candidate[ mf.frqs$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  
  bal.lst[j] <- mean(mf.gwases$Andres_candidate[ mf.gwases$LST>quantiles.lst.perm[j] & mf.gwases$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
  bal.lst.perm[j] <- mean(mf.gwases$Andres_candidate[ mf.gwases$LST_PERM2>quantiles.lst.perm[j] & mf.gwases$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
  
  bal.t[j] <- mean(mf.gwases.lrs$Andres_candidate[ mf.gwases.lrs$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
  bal.t.perm[j] <- mean(mf.gwases.lrs$Andres_candidate[ mf.gwases.lrs$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
  
  bal.unfolded.t[j] <- mean(mf.gwases.lrs$Andres_candidate[ mf.gwases.lrs$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  bal.unfolded.t.perm[j] <- mean(mf.gwases.lrs$Andres_candidate[ mf.gwases.lrs$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  
  
  print(j)
  
}


## MAF correction
quantiles.maf <- quantile(mf.frqs$MAF, prob = seq(0, 1, length = 21), type = 1)
quantiles.maf[1] <- 0.01
quantiles.maf[21] <- 0.5

bal.viability.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.viability.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.reproductive.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.reproductive.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.total.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.total.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.unfolded.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.unfolded.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.lst.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.lst.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.unfolded.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.unfolded.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)


for (i in 1:20){
  for (j in 1:100){
    mf.frqs3 <- subset(mf.frqs,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases3 <- subset(mf.gwases,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases.lrs3 <- subset(mf.gwases.lrs,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    
    bal.viability.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.viability.perm.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.viability.theory[j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.reproductive.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.reproductive.perm.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.reproductive.theory[j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.total.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.total.perm.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.total.theory[j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.unfolded.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    bal.unfolded.perm.mafcorrected[i,j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    #bal.unfolded.theory[j] <- mean(mf.frqs3$Andres_candidate[ mf.frqs3$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    
    bal.lst.mafcorrected[i,j] <- mean(mf.gwases3$Andres_candidate[ mf.gwases3$LST>quantiles.lst.perm[j] & mf.gwases3$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
    bal.lst.perm.mafcorrected[i,j] <- mean(mf.gwases3$Andres_candidate[ mf.gwases3$LST_PERM2>quantiles.lst.perm[j] & mf.gwases3$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
    
    bal.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Andres_candidate[ mf.gwases.lrs3$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
    bal.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Andres_candidate[ mf.gwases.lrs3$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
    
    bal.unfolded.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Andres_candidate[ mf.gwases.lrs3$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
    bal.unfolded.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$Andres_candidate[ mf.gwases.lrs3$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  }
  print(i)
}  


#Plots####
bal.viability.mafcorrected <- colMeans(bal.viability.mafcorrected[,1:100])
bal.viability.perm.mafcorrected <- colMeans(bal.viability.perm.mafcorrected[,1:100])
#bal.viability.theory.mafcorrected <- colMeans(bal.viability.theory.mafcorrected[,1:100])

bal.reproductive.mafcorrected <- colMeans(bal.reproductive.mafcorrected[,1:100])
bal.reproductive.perm.mafcorrected <- colMeans(bal.reproductive.perm.mafcorrected[,1:100])
#bal.reproductive.theory.mafcorrected <- colMeans(bal.reproductive.theory.mafcorrected[,1:100])

bal.total.mafcorrected <- colMeans(bal.total.mafcorrected[,1:100])
bal.total.perm.mafcorrected <- colMeans(bal.total.perm.mafcorrected[,1:100])
#bal.total.theory.mafcorrected <- colMeans(bal.total.theory.mafcorrected[,1:100])

bal.unfolded.mafcorrected <- colMeans(bal.unfolded.mafcorrected[,1:100])
bal.unfolded.perm.mafcorrected <- colMeans(bal.unfolded.perm.mafcorrected[,1:100])
#bal.unfolded.theory.mafcorrected <- colMeans(bal.unfolded.theory.mafcorrected[,1:100])

bal.lst.mafcorrected <- colMeans(bal.lst.mafcorrected[,1:100])
bal.lst.perm.mafcorrected <- colMeans(bal.lst.perm.mafcorrected[,1:100])

bal.t.mafcorrected <- colMeans(bal.t.mafcorrected[,1:100])
bal.t.perm.mafcorrected <- colMeans(bal.t.perm.mafcorrected[,1:100])

bal.unfolded.t.mafcorrected <- colMeans(bal.unfolded.t.mafcorrected[,1:100])
bal.unfolded.t.perm.mafcorrected <- colMeans(bal.unfolded.t.perm.mafcorrected[,1:100])


dd.quantiles <- data.frame(c(bal.viability,bal.viability.mafcorrected,
                             bal.viability.perm,bal.viability.perm.mafcorrected,
                             bal.reproductive,bal.reproductive.mafcorrected,
                             bal.reproductive.perm,bal.reproductive.perm.mafcorrected,
                             bal.total,bal.total.mafcorrected,
                             bal.total.perm,bal.total.perm.mafcorrected,
                             bal.unfolded,bal.unfolded.mafcorrected,
                             bal.unfolded.perm,bal.unfolded.perm.mafcorrected,
                             bal.lst,bal.lst.mafcorrected,
                             bal.lst.perm,bal.lst.perm.mafcorrected,
                             bal.t,bal.t.mafcorrected,
                             bal.t.perm,bal.t.perm.mafcorrected,
                             bal.unfolded.t,bal.unfolded.t.mafcorrected,
                             bal.unfolded.t.perm,bal.unfolded.t.perm.mafcorrected))
dd.quantiles$Type <- factor(rep(c(rep("Observed",200),rep("Permuted",200)),7))
dd.quantiles$Quantiles <- rep(c(1:100),28)
dd.quantiles$MAF_corrected <- factor(rep(c(rep("No MAF correction",100),rep("MAF corrected",100)),14))
dd.quantiles$Metric <- factor(c(rep("Adult Fst",400),rep("Reproductive Fst",400),rep("Gametic Fst",400),rep("Unfolded Fst",400),rep("Lst",400),rep("|t|",400),rep("Unfolded t",400)))
names(dd.quantiles) <- c("Proportion","Type","Quantiles","MAF_corrected","Metric")
dd.quantiles$MAF_corrected <- factor(dd.quantiles$MAF_corrected, levels = c("No MAF correction", "MAF corrected"))
dd.quantiles$Type <- factor(dd.quantiles$Type, levels = c("Permuted", "Observed"))
dd.quantiles <- dd.quantiles[order(dd.quantiles$Type,dd.quantiles$MAF_corrected),]

ggplot(subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "darkorange","Adult Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "darkorange","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Adult Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="darkorange")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Adult Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="darkorange")+
  scale_y_continuous(limits=c(0, 0.003),)+
  ylab(expression(Proportion~candidates~(Andres~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("|t|" = "forestgreen","Reproductive Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "forestgreen","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Reproductive Fst" = 21,"|t|" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="forestgreen")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Reproductive Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="forestgreen")+
  scale_y_continuous(limits=c(0, 0.003),)+
  ylab(expression(Proportion~candidates~(Andres~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Unfolded t" = "#40B0A6","Unfolded Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "#40B0A6","Observed"="darkgrey"))+
  scale_shape_manual(values=c("Unfolded Fst" = 21,"Unfolded t" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="#40B0A6")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="#40B0A6")+
  scale_y_continuous(limits=c(0, 0.003),)+
  ylab(expression(Proportion~candidates~(Andres~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Gametic Fst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "purple","Gametic Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "purple","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Gametic Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="purple")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  scale_y_continuous(limits=c(0, 0.003),)+
  ylab(expression(Proportion~candidates~(Andres~et~al.)))+
  xlab(expression(Null~quantile))



#...####
 

#DeGiorgio, Statistical tests####

#Viability Fst
summary(glm(data=mf.frqs,DeGiorgio_candidate~FST_VIABILITY_ST+MAF,family="binomial"))
#-0.005832, p=0.207

#Reproductive Fst
summary(glm(data=mf.frqs,DeGiorgio_candidate~FST_REPRODUCTIVE_ST+MAF,family="binomial"))
#0.005246, p=0.254

#Gametic Fst
summary(glm(data=mf.frqs,DeGiorgio_candidate~FST_GAMETIC_ST+MAF,family="binomial"))
#0.002775, p=0.547

#Unfolded Fst
summary(glm(data=subset(mf.frqs,UNFOLDED_FST<0),DeGiorgio_candidate~UNFOLDED_FST+MAF,family="binomial"))
#-0.00594, p=0.623
summary(glm(data=subset(mf.frqs,UNFOLDED_FST>0),DeGiorgio_candidate~UNFOLDED_FST+MAF,family="binomial"))
#0.079487, p<0.001

#Lst
summary(glm(data=mf.gwases,DeGiorgio_candidate~LST+MAF,family="binomial"))
#-3.573e+03, p=0.686

#|t|
summary(glm(data=mf.gwases.lrs,DeGiorgio_candidate~ABS_T+MAF,family="binomial"))
#-0.003157, p=0.773

#Unfolded t
summary(glm(data=subset(mf.gwases.lrs,UNFOLDED_T<0),DeGiorgio_candidate~UNFOLDED_T+MAF,family="binomial"))
#-0.01028, p=0.395885
summary(glm(data=subset(mf.gwases.lrs,UNFOLDED_T>0),DeGiorgio_candidate~UNFOLDED_T+MAF,family="binomial"))
#0.079545, p<0.001


#Quantile-proportion plot####

#No MAF correction
quantiles.fst.theory <- quantile(mf.frqs$RANDOM_CHISQ, prob = seq(0, 1, length = 101), type = 1)
quantiles.fst.theory[1] <- 0
quantiles.fst.theory[101] <- 10^6

quantiles.unfolded.theory <- quantile(mf.frqs$UNFOLDED_FST_THEORY, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.theory[1] <- -10^6
quantiles.unfolded.theory[101] <- 10^6

quantiles.lst.perm <- quantile(allele.ages2$LST_PERM2, prob = seq(0, 1, length = 101), type = 1)
quantiles.lst.perm[1] <- 0
quantiles.lst.perm[101] <- 10^6

quantiles.t.perm <- quantile(allele.ages2$ABS_T_PERM, prob = seq(0, 1, length = 101), type = 1)
quantiles.t.perm[1] <- 0
quantiles.t.perm[101] <- 10^6

quantiles.unfolded.t.perm <- quantile(allele.ages2$UNFOLDED_T, prob = seq(0, 1, length = 101), type = 1)
quantiles.unfolded.t.perm[1] <- -10^6
quantiles.unfolded.t.perm[101] <- 10^6


bal.viability <- vector("numeric")
bal.viability.perm <- vector("numeric")
bal.reproductive <- vector("numeric")
bal.reproductive.perm <- vector("numeric")
bal.total <- vector("numeric")
bal.total.perm <- vector("numeric")
bal.unfolded <- vector("numeric")
bal.unfolded.perm <- vector("numeric")
bal.lst <- vector()
bal.lst.perm <- vector()
bal.t <- vector()
bal.t.perm <- vector()
bal.unfolded.t <- vector("numeric")
bal.unfolded.t.perm <- vector("numeric")

for (j in 1:100){
  
  bal.viability[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.viability.perm[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.viability.theory[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.reproductive[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.reproductive.perm[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.reproductive.theory[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.total[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  bal.total.perm[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
  #bal.total.theory[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
  
  bal.unfolded[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  bal.unfolded.perm[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  #bal.unfolded.theory[j] <- mean(mf.frqs$DeGiorgio_candidate[ mf.frqs$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
  
  bal.lst[j] <- mean(mf.gwases$DeGiorgio_candidate[ mf.gwases$LST>quantiles.lst.perm[j] & mf.gwases$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
  bal.lst.perm[j] <- mean(mf.gwases$DeGiorgio_candidate[ mf.gwases$LST_PERM2>quantiles.lst.perm[j] & mf.gwases$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
  
  bal.t[j] <- mean(mf.gwases.lrs$DeGiorgio_candidate[ mf.gwases.lrs$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
  bal.t.perm[j] <- mean(mf.gwases.lrs$DeGiorgio_candidate[ mf.gwases.lrs$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
  
  bal.unfolded.t[j] <- mean(mf.gwases.lrs$DeGiorgio_candidate[ mf.gwases.lrs$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  bal.unfolded.t.perm[j] <- mean(mf.gwases.lrs$DeGiorgio_candidate[ mf.gwases.lrs$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  
  
  print(j)
  
}


## MAF correction
quantiles.maf <- quantile(mf.frqs$MAF, prob = seq(0, 1, length = 21), type = 1)
quantiles.maf[1] <- 0.01
quantiles.maf[21] <- 0.5

bal.viability.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.viability.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.reproductive.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.reproductive.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.total.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.total.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.unfolded.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.unfolded.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.lst.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.lst.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)

bal.unfolded.t.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)
bal.unfolded.t.perm.mafcorrected <- matrix(nrow=20,ncol=100,data=NA)


for (i in 1:20){
  for (j in 1:100){
    mf.frqs3 <- subset(mf.frqs,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases3 <- subset(mf.gwases,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    mf.gwases.lrs3 <- subset(mf.gwases.lrs,MAF>quantiles.maf[i] & MAF<quantiles.maf[(i+1)])
    
    bal.viability.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$FST_VIABILITY_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.viability.perm.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$FST_VIABILITY_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_VIABILITY_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.viability.theory[j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.reproductive.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$FST_REPRODUCTIVE_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.reproductive.perm.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$FST_REPRODUCTIVE_PERM_ST>quantiles.fst.theory[j] & mf.frqs3$FST_REPRODUCTIVE_PERM_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.reproductive.theory[j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.total.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$FST_GAMETIC_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    bal.total.perm.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$FST_GAMETIC_PERM2_ST>quantiles.fst.theory[j] & mf.frqs3$FST_GAMETIC_PERM2_ST<quantiles.fst.theory[(j+1)]],na.rm=T)
    #bal.total.theory[j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$RANDOM_CHISQ>quantiles.fst.theory[j] & mf.frqs3$RANDOM_CHISQ<quantiles.fst.theory[(j+1)]],na.rm=T)
    
    bal.unfolded.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$UNFOLDED_FST>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    bal.unfolded.perm.mafcorrected[i,j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$UNFOLDED_FST_PERM>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_PERM<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    #bal.unfolded.theory[j] <- mean(mf.frqs3$DeGiorgio_candidate[ mf.frqs3$UNFOLDED_FST_THEORY>quantiles.unfolded.theory[j] & mf.frqs3$UNFOLDED_FST_THEORY<quantiles.unfolded.theory[(j+1)]],na.rm=T)
    
    bal.lst.mafcorrected[i,j] <- mean(mf.gwases3$DeGiorgio_candidate[ mf.gwases3$LST>quantiles.lst.perm[j] & mf.gwases3$LST<quantiles.lst.perm[(j+1)]],na.rm=T)
    bal.lst.perm.mafcorrected[i,j] <- mean(mf.gwases3$DeGiorgio_candidate[ mf.gwases3$LST_PERM2>quantiles.lst.perm[j] & mf.gwases3$LST_PERM2<quantiles.lst.perm[(j+1)]],na.rm=T)
    
    bal.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$DeGiorgio_candidate[ mf.gwases.lrs3$ABS_T>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T<quantiles.t.perm[(j+1)]],na.rm=T)
    bal.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$DeGiorgio_candidate[ mf.gwases.lrs3$ABS_T_PERM>quantiles.t.perm[j] & mf.gwases.lrs3$ABS_T_PERM<quantiles.t.perm[(j+1)]],na.rm=T)
    
    bal.unfolded.t.mafcorrected[i,j] <- mean(mf.gwases.lrs3$DeGiorgio_candidate[ mf.gwases.lrs3$UNFOLDED_T>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
    bal.unfolded.t.perm.mafcorrected[i,j] <- mean(mf.gwases.lrs3$DeGiorgio_candidate[ mf.gwases.lrs3$UNFOLDED_T_PERM>quantiles.unfolded.t.perm[j] & mf.gwases.lrs3$UNFOLDED_T_PERM<quantiles.unfolded.t.perm[(j+1)]],na.rm=T)
  }
  print(i)
}  


#Plots####
bal.viability.mafcorrected <- colMeans(bal.viability.mafcorrected[,1:100])
bal.viability.perm.mafcorrected <- colMeans(bal.viability.perm.mafcorrected[,1:100])
#bal.viability.theory.mafcorrected <- colMeans(bal.viability.theory.mafcorrected[,1:100])

bal.reproductive.mafcorrected <- colMeans(bal.reproductive.mafcorrected[,1:100])
bal.reproductive.perm.mafcorrected <- colMeans(bal.reproductive.perm.mafcorrected[,1:100])
#bal.reproductive.theory.mafcorrected <- colMeans(bal.reproductive.theory.mafcorrected[,1:100])

bal.total.mafcorrected <- colMeans(bal.total.mafcorrected[,1:100])
bal.total.perm.mafcorrected <- colMeans(bal.total.perm.mafcorrected[,1:100])
#bal.total.theory.mafcorrected <- colMeans(bal.total.theory.mafcorrected[,1:100])

bal.unfolded.mafcorrected <- colMeans(bal.unfolded.mafcorrected[,1:100])
bal.unfolded.perm.mafcorrected <- colMeans(bal.unfolded.perm.mafcorrected[,1:100])
#bal.unfolded.theory.mafcorrected <- colMeans(bal.unfolded.theory.mafcorrected[,1:100])

bal.lst.mafcorrected <- colMeans(bal.lst.mafcorrected[,1:100])
bal.lst.perm.mafcorrected <- colMeans(bal.lst.perm.mafcorrected[,1:100])

bal.t.mafcorrected <- colMeans(bal.t.mafcorrected[,1:100])
bal.t.perm.mafcorrected <- colMeans(bal.t.perm.mafcorrected[,1:100])

bal.unfolded.t.mafcorrected <- colMeans(bal.unfolded.t.mafcorrected[,1:100])
bal.unfolded.t.perm.mafcorrected <- colMeans(bal.unfolded.t.perm.mafcorrected[,1:100])


dd.quantiles <- data.frame(c(bal.viability,bal.viability.mafcorrected,
                             bal.viability.perm,bal.viability.perm.mafcorrected,
                             bal.reproductive,bal.reproductive.mafcorrected,
                             bal.reproductive.perm,bal.reproductive.perm.mafcorrected,
                             bal.total,bal.total.mafcorrected,
                             bal.total.perm,bal.total.perm.mafcorrected,
                             bal.unfolded,bal.unfolded.mafcorrected,
                             bal.unfolded.perm,bal.unfolded.perm.mafcorrected,
                             bal.lst,bal.lst.mafcorrected,
                             bal.lst.perm,bal.lst.perm.mafcorrected,
                             bal.t,bal.t.mafcorrected,
                             bal.t.perm,bal.t.perm.mafcorrected,
                             bal.unfolded.t,bal.unfolded.t.mafcorrected,
                             bal.unfolded.t.perm,bal.unfolded.t.perm.mafcorrected))
dd.quantiles$Type <- factor(rep(c(rep("Observed",200),rep("Permuted",200)),7))
dd.quantiles$Quantiles <- rep(c(1:100),28)
dd.quantiles$MAF_corrected <- factor(rep(c(rep("No MAF correction",100),rep("MAF corrected",100)),14))
dd.quantiles$Metric <- factor(c(rep("Adult Fst",400),rep("Reproductive Fst",400),rep("Gametic Fst",400),rep("Unfolded Fst",400),rep("Lst",400),rep("|t|",400),rep("Unfolded t",400)))
names(dd.quantiles) <- c("Proportion","Type","Quantiles","MAF_corrected","Metric")
dd.quantiles$MAF_corrected <- factor(dd.quantiles$MAF_corrected, levels = c("No MAF correction", "MAF corrected"))
dd.quantiles$Type <- factor(dd.quantiles$Type, levels = c("Permuted", "Observed"))
dd.quantiles <- dd.quantiles[order(dd.quantiles$Type,dd.quantiles$MAF_corrected),]

ggplot(subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "darkorange","Adult Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "darkorange","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Adult Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Adult Fst" |Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="darkorange")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Adult Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Lst") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="darkorange")+
  scale_y_continuous(limits=c(0.015, 0.03),)+
  ylab(expression(Proportion~candidates~(DeGiorgio~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("|t|" = "forestgreen","Reproductive Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "forestgreen","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Reproductive Fst" = 21,"|t|" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Reproductive Fst" |Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="forestgreen")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Reproductive Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="|t|") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="forestgreen")+
  scale_y_continuous(limits=c(0.015, 0.03),)+
  ylab(expression(Proportion~candidates~(DeGiorgio~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Unfolded t" = "#40B0A6","Unfolded Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "#40B0A6","Observed"="darkgrey"))+
  scale_shape_manual(values=c("Unfolded Fst" = 21,"Unfolded t" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Unfolded Fst" |Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="#40B0A6")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Unfolded t") & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.4,col="#40B0A6")+
  scale_y_continuous(limits=c(0.015, 0.03),)+
  ylab(expression(Proportion~candidates~(DeGiorgio~et~al.)))+
  xlab(expression(Null~quantile))

ggplot(subset(dd.quantiles,(Metric=="Gametic Fst") & MAF_corrected=="MAF corrected" & Type=="Observed"),aes(x=Quantiles,y=Proportion,col=Metric,fill=Type,shape=Metric))+
  geom_point(size=3)+
  #geom_line(aes(color=Type),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values=c("Lst" = "purple","Gametic Fst"="black"))+
  scale_fill_manual(values=c("Observed" = "purple","Permuted"="darkgrey"))+
  scale_shape_manual(values=c("Gametic Fst" = 21,"Lst" = 23))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),legend.position="none",title = element_text(size=20))+
  geom_smooth(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), method='loess',formula=y~x,se=T,size=1.5,col="purple")+
  geom_ribbon(data=subset(dd.quantiles,(Metric=="Gametic Fst" ) & MAF_corrected=="MAF corrected" & Type=="Observed"), stat="smooth",method='loess',formula=y~x,se=T,alpha=0,linetype="dotted",size=0.6,col="black")+
  scale_y_continuous(limits=c(0.015, 0.03),)+
  ylab(expression(Proportion~candidates~(DeGiorgio~et~al.)))+
  xlab(expression(Null~quantile))



#...####
#...####
#...####


#1000G, Tajima's D####

#Import####

#Imputed
Tajima.D.YRI <- read.table(paste0(dir,"/1000_genomes_data/chrAUTO.YRI.Tajima.D"),h=T)
Tajima.D.YRI <- data.frame(sapply( Tajima.D.YRI, function(x)as.numeric(as.character(x)) ))
names(Tajima.D.YRI)[3:4] <- c("N_SNPS_YRI","TajimaD_YRI")
Tajima.D.GIH <- read.table(paste0(dir,"/1000_genomes_data/chrAUTO.GIH.Tajima.D"),h=T)
Tajima.D.GIH <- data.frame(sapply( Tajima.D.GIH, function(x)as.numeric(as.character(x)) ))
names(Tajima.D.GIH)[3:4] <- c("N_SNPS_GIH","TajimaD_GIH")


#10kb windows
mf.frqs$BIN_START <- (cut(mf.frqs$POS,breaks=seq(0,249216384,by=10000),labels = F)-1)*10000
mf.gwases$BIN_START <- (cut(mf.gwases$POS,breaks=seq(0,249216384,by=10000),labels = F)-1)*10000
mf.gwases.lrs$BIN_START <- (cut(mf.gwases.lrs$POS,breaks=seq(0,249216384,by=10000),labels = F)-1)*10000

#Take mean of each metric of sex-differential selection and MAF in each 10kb bin
mf.frqs.agg <- aggregate(cbind(FST_VIABILITY_ST,FST_VIABILITY_PERM2_ST,FST_REPRODUCTIVE_ST,FST_REPRODUCTIVE_PERM_ST,FST_GAMETIC_ST,FST_GAMETIC_PERM2_ST,MAF,RANDOM_CHISQ,UNFOLDED_FST_THEORY,UNFOLDED_FST,UNFOLDED_FST_PERM)~CHROM+BIN_START,data=mf.frqs,mean,na.rm=T)
mf.frqs.neg <- subset(mf.frqs,UNFOLDED_FST<0)
mf.frqs.agg.neg <- aggregate(cbind(UNFOLDED_FST)~CHROM+BIN_START,data=mf.frqs.neg,mean,na.rm=T)
names(mf.frqs.agg.neg)[3] <- "UNFOLDED_FST_NEGATIVE"
mf.frqs.pos <- subset(mf.frqs,UNFOLDED_FST>0)
mf.frqs.agg.pos <- aggregate(cbind(UNFOLDED_FST)~CHROM+BIN_START,data=mf.frqs.pos,mean,na.rm=T)
names(mf.frqs.agg.pos)[3] <- "UNFOLDED_FST_POSITIVE"

mf.gwases.agg <- aggregate(cbind(LST,LST_PERM2,MAF)~CHROM+BIN_START,data=mf.gwases,mean,na.rm=T)

mf.gwases.lrs.agg <- aggregate(cbind(ABS_T,ABS_T_PERM,MAF,UNFOLDED_T_PERM,UNFOLDED_T)~CHROM+BIN_START,data=mf.gwases.lrs,mean,na.rm=T)
mf.gwases.lrs.neg <- subset(mf.gwases.lrs,UNFOLDED_T<0)
mf.gwases.lrs.agg.neg <- aggregate(cbind(UNFOLDED_T)~CHROM+BIN_START,data=mf.gwases.lrs.neg,mean,na.rm=T)
names(mf.gwases.lrs.agg.neg)[3] <- "UNFOLDED_T_NEGATIVE"
mf.gwases.lrs.pos <- subset(mf.gwases.lrs,UNFOLDED_T>0)
mf.gwases.lrs.agg.pos <- aggregate(cbind(UNFOLDED_T)~CHROM+BIN_START,data=mf.gwases.lrs.pos,mean,na.rm=T)
names(mf.gwases.lrs.agg.pos)[3] <- "UNFOLDED_T_POSITIVE"

#Merge sex-differential and MAF info with Tajima's D values
mf.frqs.agg2 <- Reduce(function(...) merge(..., all.x=T,by=c("CHROM","BIN_START")),list(mf.frqs.agg,Tajima.D.YRI,Tajima.D.GIH,mf.frqs.agg.pos,mf.frqs.agg.neg))
mf.gwases.agg2 <- Reduce(function(...) merge(..., all.x=T,by=c("CHROM","BIN_START")),list(mf.gwases.agg,Tajima.D.YRI,Tajima.D.GIH))
mf.gwases.lrs.agg2 <- Reduce(function(...) merge(..., all.x=T,by=c("CHROM","BIN_START")),list(mf.gwases.lrs.agg,Tajima.D.YRI,Tajima.D.GIH,mf.gwases.lrs.agg.pos,mf.gwases.lrs.agg.neg))

#Statistical tests####

#Viability Fst
summary(glm(data=mf.frqs.agg2,TajimaD_YRI~FST_VIABILITY_ST+MAF))
#-0.003782, p=0.0151
summary(glm(data=mf.frqs.agg2,TajimaD_GIH~FST_VIABILITY_ST+MAF))
#-0.0007409, p=0.737

#Reproductive Fst
summary(glm(data=mf.frqs.agg2,TajimaD_YRI~FST_REPRODUCTIVE_ST+MAF))
#0.000665, p=0.677
summary(glm(data=mf.frqs.agg2,TajimaD_GIH~FST_REPRODUCTIVE_ST+MAF))
#0.004997, p=0.0274

#Gametic Fst
summary(glm(data=mf.frqs.agg2,TajimaD_YRI~FST_GAMETIC_ST+MAF))
#-0.0007156, p=0.652
summary(glm(data=mf.frqs.agg2,TajimaD_GIH~FST_GAMETIC_ST+MAF))
#0.003843, p=0.0879

#Unfolded Fst
summary(glm(data=mf.frqs.agg2,TajimaD_YRI~UNFOLDED_FST_NEGATIVE+MAF))
#0.000687, 0.779
summary(glm(data=mf.frqs.agg2,TajimaD_GIH~UNFOLDED_FST_NEGATIVE+MAF))
#-0.0003981, p=0.908
summary(glm(data=mf.frqs.agg2,TajimaD_YRI~UNFOLDED_FST_POSITIVE+MAF))
#-0.006919, p=0.00214
summary(glm(data=mf.frqs.agg2,TajimaD_GIH~UNFOLDED_FST_POSITIVE+MAF))
#-0.008728, p=0.006

#Lst
summary(glm(data=mf.gwases.agg2,TajimaD_YRI~LST+MAF))
#-4.436e+03, p=0.139
summary(glm(data=mf.gwases.agg2,TajimaD_GIH~LST+MAF))
#1.131e+03, p=0.791

#|t|
summary(glm(data=mf.gwases.lrs.agg2,TajimaD_YRI~ABS_T+MAF))
#0.004345, p=0.248
summary(glm(data=mf.gwases.lrs.agg2,TajimaD_GIH~ABS_T+MAF))
#0.007511, p=0.16

#Unfolded t
summary(glm(data=subset(mf.gwases.lrs.agg2),TajimaD_YRI~UNFOLDED_T_NEGATIVE+MAF))
#-3.133e-05, p=0.99
summary(glm(data=subset(mf.gwases.lrs.agg2),TajimaD_GIH~UNFOLDED_T_NEGATIVE+MAF))
#0.0000735, p=0.983
summary(glm(data=subset(mf.gwases.lrs.agg2),TajimaD_YRI~UNFOLDED_T_POSITIVE+MAF))
#-0.005514, p=0.0138
summary(glm(data=subset(mf.gwases.lrs.agg2),TajimaD_GIH~UNFOLDED_T_POSITIVE+MAF))
#-0.006751, p=0.0325
