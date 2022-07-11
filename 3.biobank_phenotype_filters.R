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

## Phenotype filters####
bd <- read.table(paste0(dir,"ukb37327.tab"), head=TRUE, sep="\t")

lvl.0009 <- c(0,1)
lbl.0009 <- c("Female","Male")
bd$f.31.0.0 <- ordered(bd$f.31.0.0, levels=lvl.0009, labels=lbl.0009)
lvl.0008 <- c(1,2,3,4,5,6,7,8,9,10,11,12)
lbl.0008 <- c("January","February","March","April","May","June","July","August","September","October","November","December")
bd$f.52.0.0 <- ordered(bd$f.52.0.0, levels=lvl.0008, labels=lbl.0008)
bd$f.53.0.0 <- as.Date(bd$f.53.0.0)
bd$f.53.1.0 <- as.Date(bd$f.53.1.0)
bd$f.53.2.0 <- as.Date(bd$f.53.2.0)
bd$f.53.3.0 <- as.Date(bd$f.53.3.0)
bd$f.191.0.0 <- as.Date(bd$f.191.0.0)
lvl.100306 <- c(-3,-2,-1)
lbl.100306 <- c("Prefer not to answer","Never went to school","Do not know")
lvl.100508 <- c(-3,-1,1,2,3,4)
lbl.100508 <- c("Prefer not to answer","Do not know","Excellent","Good","Fair","Poor")
bd$f.2178.0.0 <- ordered(bd$f.2178.0.0, levels=lvl.100508, labels=lbl.100508)
bd$f.2178.1.0 <- ordered(bd$f.2178.1.0, levels=lvl.100508, labels=lbl.100508)
bd$f.2178.2.0 <- ordered(bd$f.2178.2.0, levels=lvl.100508, labels=lbl.100508)
bd$f.2178.3.0 <- ordered(bd$f.2178.3.0, levels=lvl.100508, labels=lbl.100508)
lvl.100570 <- c(-3,-1,1,2,3)
lbl.100570 <- c("Prefer not to answer","Do not know","Younger than average","About average age","Older than average")
bd$f.2385.0.0 <- ordered(bd$f.2385.0.0, levels=lvl.100570, labels=lbl.100570)
bd$f.2385.1.0 <- ordered(bd$f.2385.1.0, levels=lvl.100570, labels=lbl.100570)
bd$f.2385.2.0 <- ordered(bd$f.2385.2.0, levels=lvl.100570, labels=lbl.100570)
lvl.100291 <- c(-3,-1)
lbl.100291 <- c("Prefer not to answer","Do not know")
lvl.100349 <- c(-3,-1,0,1)
lbl.100349 <- c("Prefer not to answer","Do not know","No","Yes")
bd$f.2443.0.0 <- ordered(bd$f.2443.0.0, levels=lvl.100349, labels=lbl.100349)
bd$f.2443.1.0 <- ordered(bd$f.2443.1.0, levels=lvl.100349, labels=lbl.100349)
bd$f.2443.2.0 <- ordered(bd$f.2443.2.0, levels=lvl.100349, labels=lbl.100349)
bd$f.2443.3.0 <- ordered(bd$f.2443.3.0, levels=lvl.100349, labels=lbl.100349)
lvl.100584 <- c(-3)
lbl.100584 <- c("Prefer not to answer")
lvl.100582 <- c(-6,-3,-1)
lbl.100582 <- c("Irregular cycle","Prefer not to answer","Do not know")
bd$f.3799.0.0 <- ordered(bd$f.3799.0.0, levels=lvl.100349, labels=lbl.100349)
bd$f.3799.1.0 <- ordered(bd$f.3799.1.0, levels=lvl.100349, labels=lbl.100349)
bd$f.3799.2.0 <- ordered(bd$f.3799.2.0, levels=lvl.100349, labels=lbl.100349)
bd$f.3799.3.0 <- ordered(bd$f.3799.3.0, levels=lvl.100349, labels=lbl.100349)
lvl.100696 <- c(-1)
lbl.100696 <- c("Abandoned")
lvl.100305 <- c(-7,-3,1,2,3,4,5,6)
lbl.100305 <- c("None of the above","Prefer not to answer","College or University degree","A levels/AS levels or equivalent","O levels/GCSEs or equivalent","CSEs or equivalent","NVQ or HND or HNC or equivalent","Other professional qualifications eg: nursing, teaching")
bd$f.6138.0.0 <- ordered(bd$f.6138.0.0, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.0.1 <- ordered(bd$f.6138.0.1, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.0.2 <- ordered(bd$f.6138.0.2, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.0.3 <- ordered(bd$f.6138.0.3, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.0.4 <- ordered(bd$f.6138.0.4, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.0.5 <- ordered(bd$f.6138.0.5, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.1.0 <- ordered(bd$f.6138.1.0, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.1.1 <- ordered(bd$f.6138.1.1, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.1.2 <- ordered(bd$f.6138.1.2, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.1.3 <- ordered(bd$f.6138.1.3, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.1.4 <- ordered(bd$f.6138.1.4, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.1.5 <- ordered(bd$f.6138.1.5, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.2.0 <- ordered(bd$f.6138.2.0, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.2.1 <- ordered(bd$f.6138.2.1, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.2.2 <- ordered(bd$f.6138.2.2, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.2.3 <- ordered(bd$f.6138.2.3, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.2.4 <- ordered(bd$f.6138.2.4, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.2.5 <- ordered(bd$f.6138.2.5, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.3.0 <- ordered(bd$f.6138.3.0, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.3.1 <- ordered(bd$f.6138.3.1, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.3.2 <- ordered(bd$f.6138.3.2, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.3.3 <- ordered(bd$f.6138.3.3, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.3.4 <- ordered(bd$f.6138.3.4, levels=lvl.100305, labels=lbl.100305)
bd$f.6138.3.5 <- ordered(bd$f.6138.3.5, levels=lvl.100305, labels=lbl.100305)
lvl.0013 <- c(-3,-1)
lbl.0013 <- c("Preferred not to answer","Date uncertain or unknown")
lvl.1543 <- c(0,1,2)
lbl.1543 <- c("non-myopic","moderate/low myopia","highly myopic")
bd$f.20262.0.0 <- ordered(bd$f.20262.0.0, levels=lvl.1543, labels=lbl.1543)
lvl.1001 <- c(-3,-1,1,2,3,4,5,6,1001,1002,1003,2001,2002,2003,2004,3001,3002,3003,3004,4001,4002,4003)
lbl.1001 <- c("Prefer not to answer","Do not know","White","Mixed","Asian or Asian British","Black or Black British","Chinese","Other ethnic group","British","Irish","Any other white background","White and Black Caribbean","White and Black African","White and Asian","Any other mixed background","Indian","Pakistani","Bangladeshi","Any other Asian background","Caribbean","African","Any other Black background")
bd$f.21000.0.0 <- ordered(bd$f.21000.0.0, levels=lvl.1001, labels=lbl.1001)
bd$f.21000.1.0 <- ordered(bd$f.21000.1.0, levels=lvl.1001, labels=lbl.1001)
bd$f.21000.2.0 <- ordered(bd$f.21000.2.0, levels=lvl.1001, labels=lbl.1001)
bd$f.22001.0.0 <- ordered(bd$f.22001.0.0, levels=lvl.0009, labels=lbl.0009)
lvl.1002 <- c(1)
lbl.1002 <- c("Caucasian")
bd$f.22006.0.0 <- ordered(bd$f.22006.0.0, levels=lvl.1002, labels=lbl.1002)
lvl.0001 <- c(1)
lbl.0001 <- c("Yes")
bd$f.22019.0.0 <- ordered(bd$f.22019.0.0, levels=lvl.0001, labels=lbl.0001)
bd$f.22020.0.0 <- ordered(bd$f.22020.0.0, levels=lvl.0001, labels=lbl.0001)
lvl.0682 <- c(-1,0,1,10)
lbl.0682 <- c("Participant excluded from kinship inference process","No kinship found","At least one relative identified","Ten or more third-degree relatives identified")
bd$f.22021.0.0 <- ordered(bd$f.22021.0.0, levels=lvl.0682, labels=lbl.0682)
bd$f.22027.0.0 <- ordered(bd$f.22027.0.0, levels=lvl.0001, labels=lbl.0001)
lvl.100264 <- c(0,1)
lbl.100264 <- c("No","Yes")
bd$f.22028.0.0 <- ordered(bd$f.22028.0.0, levels=lvl.100264, labels=lbl.100264)
bd$f.22029.0.0 <- ordered(bd$f.22029.0.0, levels=lvl.100264, labels=lbl.100264)
bd$f.22030.0.0 <- ordered(bd$f.22030.0.0, levels=lvl.100264, labels=lbl.100264)
bd$f.40000.0.0 <- as.Date(bd$f.40000.0.0)
bd$f.40000.1.0 <- as.Date(bd$f.40000.1.0)
bd$f.40005.0.0 <- as.Date(bd$f.40005.0.0)
bd$f.40005.1.0 <- as.Date(bd$f.40005.1.0)
bd$f.40005.2.0 <- as.Date(bd$f.40005.2.0)
bd$f.40005.3.0 <- as.Date(bd$f.40005.3.0)
bd$f.40005.4.0 <- as.Date(bd$f.40005.4.0)
bd$f.40005.5.0 <- as.Date(bd$f.40005.5.0)
bd$f.40005.6.0 <- as.Date(bd$f.40005.6.0)
bd$f.40005.7.0 <- as.Date(bd$f.40005.7.0)
bd$f.40005.8.0 <- as.Date(bd$f.40005.8.0)
bd$f.40005.9.0 <- as.Date(bd$f.40005.9.0)
bd$f.40005.10.0 <- as.Date(bd$f.40005.10.0)
bd$f.40005.11.0 <- as.Date(bd$f.40005.11.0)
bd$f.40005.12.0 <- as.Date(bd$f.40005.12.0)
bd$f.40005.13.0 <- as.Date(bd$f.40005.13.0)
bd$f.40005.14.0 <- as.Date(bd$f.40005.14.0)
bd$f.40005.15.0 <- as.Date(bd$f.40005.15.0)
bd$f.40005.16.0 <- as.Date(bd$f.40005.16.0)
lvl.0261 <- c(1,2,7,19,55)
lbl.0261 <- c("IC Death Format (2011 and earlier)","IC Death Format (2012 onwards)","Scottish Morbidity Record (SMR)","Scottish Morbidity Record (SMR) 99B - 2015","IC Scottish deaths (2017 onwards)")
bd$f.40018.0.0 <- ordered(bd$f.40018.0.0, levels=lvl.0261, labels=lbl.0261)
bd$f.40018.1.0 <- ordered(bd$f.40018.1.0, levels=lvl.0261, labels=lbl.0261)
lvl.0272 <- c(1)
lbl.0272 <- c("Date is unknown")
bd$f.42020.0.0 <- as.Date(bd$f.42020.0.0)
bd$f.42032.0.0 <- as.Date(bd$f.42032.0.0)

## Write out non-sample-filtered phenotype file
#setwd("/Users/fruz0001/Documents/biobank/")
#write.table(bd,"pheno.txt",quote=F,row.names=F,sep="\t")

## Basic filters####

#Keep only unrelated individuals
bd1 <- subset(bd,f.22021.0.0=="No kinship found")
#Remove non-'White British ancestry' individuals
bd2 <- subset(bd1,f.22006.0.0=="Caucasian")
rm(bd1)
#Remove samples with high heterozygosity or missing rates
bd3 <- subset(bd2,is.na(f.22027.0.0))
rm(bd2)
#Remove samples where inferred sex does not match reported sex
bd4 <- subset(bd3,f.31.0.0==f.22001.0.0)
rm(bd3)
#Remove aneuploid samples
bd5 <- subset(bd4,is.na(f.22019.0.0))
rm(bd4)
#write.table(bd5,"pheno_filtered_v1.txt",quote=F,row.names=F,sep="\t")
#write.table(bd5[,c(1,1)],"indiv_filtered_v1.txt",quote=F,row.names=F,sep="\t",col.names=F)

#Only keep individuals aged 45 and above 
bd6 <- subset(bd5,f.21003.0.0>44)
names(bd6)[1] <- "FID"
rm(bd5)
rm(bd)
#write.table(bd6,"pheno_filtered_v2.txt",quote=F,row.names=F,sep="\t")
#write.table(bd6[,c(1,1)],"ukb_phenotypes_v3/indiv_filtered_v2.txt",quote=F,row.names=F,sep="\t",col.names=F)

#Fitness filters####
indiv <- read.table(paste0(dir,"ukb_phenotypes_v3/indiv_filtered_v2.txt"))
names(indiv) <- c("FID","IID")
indiv.extra <- merge(indiv,bd6[c("FID","f.31.0.0","f.22001.0.0","f.2405.0.0","f.2405.1.0","f.2405.2.0","f.2405.3.0","f.2734.0.0","f.2734.1.0","f.2734.2.0","f.2734.3.0")],by="FID",all.x=T)
## Males, quality-filtering fitness data
indiv.extra.m <- subset(indiv.extra,f.31.0.0=="Male")
#Maximum value among four sampling points
indiv.extra.m$Maximum_value <- apply(indiv.extra.m[,c("f.2405.0.0","f.2405.1.0","f.2405.2.0","f.2405.3.0")],1,function(x)max(x,na.rm=T))
#Is there a negative value among any sampling point (-1=do not know, -3=prefer not to answer)
indiv.extra.m$Negative_value <- with(indiv.extra.m,ifelse(f.2405.0.0<0 | f.2405.1.0<0 | f.2405.2.0<0 | f.2405.3.0<0,"Yes","No"))
#Is the condition: number of children (time point 4) >= number of children (time point 3) >= number of children (time point 2) >= number of children (time point 1) met?  
indiv.extra.m$Lost_children <- with(indiv.extra.m,ifelse((f.2405.3.0<f.2405.2.0)|(f.2405.3.0<f.2405.1.0) |(f.2405.2.0<f.2405.1.0) |(f.2405.2.0<f.2405.0.0)|(f.2405.3.0<f.2405.0.0)|(f.2405.1.0<f.2405.0.0),"Yes","No"))
#Remove unreliable individuals (negative sampling value, 'lost' children, non-finite because no value inputted), and pick the maximum value among those that remain
indiv.extra.m$f.2405.max <- as.numeric(with(indiv.extra.m,ifelse((is.na(Lost_children) | Lost_children=="No") & (is.na(Negative_value) | Negative_value=="No") & is.finite(Maximum_value),Maximum_value,"Bad")))
#Implausible number of children
indiv.extra.m$f.2405.definitive <- with(indiv.extra.m,ifelse(f.2405.max>=20,NA,f.2405.max))
## Females, quality-filtering fitness data
indiv.extra.f <- subset(indiv.extra,f.31.0.0=="Female")
#Maximum value among four sampling points
indiv.extra.f$Maximum_value <- apply(indiv.extra.f[,c("f.2734.0.0","f.2734.1.0","f.2734.2.0","f.2734.3.0")],1,function(x)max(x,na.rm=T))
#Is there a negative value among any sampling point (-1=do not know, -3=prefer not to answer)
indiv.extra.f$Negative_value <- with(indiv.extra.f,ifelse(f.2734.0.0<0 | f.2734.1.0<0 | f.2734.2.0<0 | f.2734.3.0<0,"Yes","No"))
#Is the condition: number of children (time point 4) >= number of children (time point 3) >= number of children (time point 2) >= number of children (time point 1) met?  
indiv.extra.f$Lost_children <- with(indiv.extra.f,ifelse((f.2734.3.0<f.2734.2.0)|(f.2734.3.0<f.2734.1.0) |(f.2734.2.0<f.2734.1.0) |(f.2734.2.0<f.2734.0.0)|(f.2734.3.0<f.2734.0.0)|(f.2734.1.0<f.2734.0.0),"Yes","No"))
#Remove unreliable individuals (negative sampling value, lost children, non-finite because no value inputted) and pick the maximum value among those that remain
indiv.extra.f$f.2734.max <- as.numeric(with(indiv.extra.f,ifelse((is.na(Lost_children) | Lost_children=="No")&(is.na(Negative_value) | Negative_value=="No")&is.finite(Maximum_value),Maximum_value,"Bad")))
#Implausible number of children
indiv.extra.f$f.2734.definitive <- with(indiv.extra.f,ifelse(f.2734.max>=20,NA,f.2734.max))

#Observed data####

#Phenotype file, excluding any individual with NA fitness
indiv.extra.m2 <- indiv.extra.m[,c("FID","IID","f.31.0.0","f.22001.0.0","f.2405.definitive")]
names(indiv.extra.m2)[5] <- "Number_of_children"
indiv.extra.f2 <- indiv.extra.f[,c("FID","IID","f.31.0.0","f.22001.0.0","f.2734.definitive")]
names(indiv.extra.f2)[5] <- "Number_of_children"
indiv.extra2 <- rbind(indiv.extra.m2,indiv.extra.f2)
indiv.extra2 <- indiv.extra2[order(indiv.extra2$FID),]
bd9 <- merge(bd6,indiv.extra2[,c("FID","IID","Number_of_children")],by=c("FID"),all.x=T)
bd9 <- bd9[!is.na(bd9$Number_of_children),]
#write.table(bd9,paste0(dir,"/ukb_phenotypes_v3/pheno_filtered_v3.txt"),quote=F,row.names=F,sep="\t")
#write.table(bd9[,c(1,1)],paste0(dir,"/ukb_phenotypes_v3/indiv_filtered_v3.txt"),quote=F,row.names=F,sep="\t",col.names=F)


#Phenotype file with FID/IIDs for each sex and number of children
bd9$SexChildren <- paste0(bd9$f.31.0.0,bd9$Number_of_children)
#write.table(bd9[,c("FID","IID","SexChildren")],paste0(dir,"ukb_phenotypes_v3/indiv_filtered_v3_SexChildren_cluster.txt"),quote=F,row.names=F)


#BOLT-LMM

#Phenotype file, excluding any individual with NA fitness, for input in BOLT-LMM 
# (white-space delimited, with only  relevant phenotypes for logistic regression analysis, i.e. Sex=f.31.0.0.numeric, Assessment centre=f.54.0.0, Age=f.21003.0.0, + 20 PCs)
bd9$f.31.0.0.numeric <- ifelse(bd9$f.31.0.0=="Female",2,1)
#Import principal components file
pcs <- fread(paste0(dir,"/ukb_kinship_and_pca_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.eigenvec"))
names(pcs)[1] <- "FID"
bd12 <- merge(bd9,pcs,all.x=T,by=c("FID","IID"))
#write.table(bd12[,c("FID","IID","f.31.0.0.numeric","f.54.0.0","f.21003.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")],paste0(dir,"/ukb_phenotypes_v3/pheno_filtered_v3_for_boltlmm.txt"),quote=F,row.names=F,sep=" ")
#write.table(bd12[,c("FID","IID","f.31.0.0.numeric","f.54.0.0","f.21003.0.0","Number_of_children","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")],paste0(dir,"/ukb_phenotypes_v3/pheno_filtered_v3_for_boltlmm_nchildren.txt"),quote=F,row.names=F,sep=" ")
rm(bd12)

#List of male and female individuals (with fitness data)
bd10 <- bd9[,c(1,1:ncol(bd9))]
names(bd10)[1:2] <- c("FID","IID")
rm(bd9)
#write.table(bd10[bd10$f.31.0.0=="Male",c(1,2)],paste0(dir,"ukb_phenotypes/list_of_male_ids_v3.txt"),quote=F,row.names=F,col.names=F,sep=" ")
#write.table(bd10[bd10$f.31.0.0=="Female",c(1,2)],paste0(dir,"ukb_phenotypes/list_of_female_ids_v3.txt"),quote=F,row.names=F,col.names=F,sep=" ")

#List of individuals without fitness data
bd11 <- merge(bd6,indiv.extra2[,c("FID","IID","Number_of_children")],by=c("FID"),all.x=T)
#write.table(bd11[is.na(bd11$Number_of_children),c("FID","IID")],paste0(dir,"/ukb_phenotypes_v3/indiv_without_LRS_v3.txt"),quote=F,row.names=F,col.names=F,sep=" ")
rm(bd6)


## Permuted LRS####
#(sample number of offspring *within* each sex, once, without replacement)
indiv3 <- read.table(paste0(dir,"/ukb_phenotypes_v3/indiv_filtered_v3.txt"))
names(indiv3) <- c("FID","IID")
indiv.extra3 <- merge(indiv3,bd10[c("FID","f.31.0.0","f.31.0.0.numeric","f.22001.0.0","f.54.0.0","f.21003.0.0","Number_of_children")],by="FID",all.x=T)
indiv.extra3.m <- subset(indiv.extra3,f.31.0.0=="Male")
indiv.extra3.f <- subset(indiv.extra3,f.31.0.0=="Female")
set.seed(123)
indiv.extra3.m$Number_of_children_permuted <- sample(indiv.extra3.m$Number_of_children,replace=F,size = nrow(indiv.extra3.m))
set.seed(123)
indiv.extra3.f$Number_of_children_permuted <- sample(indiv.extra3.f$Number_of_children,replace=F,size = nrow(indiv.extra3.f))

## Phenotype file with FID/IIDs for each sex and number of children, permuted
indiv.extra3 <- rbind(indiv.extra3.m,indiv.extra3.f)
indiv.extra3 <- indiv.extra3[order(indiv.extra3$FID),]
indiv.extra3$SexChildrenPerm <- paste0(indiv.extra3$f.31.0.0,indiv.extra3$Number_of_children_permuted)
#write.table(indiv.extra3[,c("FID","IID","SexChildrenPerm")],paste0(dir,"/ukb_phenotypes_v3/indiv_filtered_v3_SexChildrenPerm_cluster.txt"),quote=F,row.names=F)

##BOLT-LMM
#Phenotype file, excluding any individual with NA fitness, for input in BOLT-LMM 
# (white-space delimited, with only relevant phenotypes for GWAS analysis, i.e. Sex=f.31.0.0.numeric, Assessment centre=f.54.0.0, Age=f.21003.0.0, Offspring=Number_of_children, PermutedOffspring=Number_of_children_permuted, + 20 PCs)
indiv.extra3.bolt <- merge(indiv.extra3,pcs,all.x=T,by=c("FID","IID"))
#write.table(indiv.extra3.bolt[,c("FID","IID","f.31.0.0.numeric","f.54.0.0","f.21003.0.0","Number_of_children","Number_of_children_permuted","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")],paste0(dir,"ukb_phenotypes_v3/pheno_filtered_v3_perm1_for_boltlmm_nchildren.txt"),quote=F,row.names=F,sep=" ")





## Permuted sex####
#(sample sex, once, without replacement)
set.seed(123)
indiv.extra3$SexPerm2 <- sample(indiv.extra3$f.31.0.0,replace=F,size = nrow(indiv.extra3))

#Phenotype file with FID/IIDs for each permuted sex and number of children
indiv.extra3$SexChildrenPerm2 <- paste0(indiv.extra3$SexPerm2,indiv.extra3$Number_of_children)
#write.table(indiv.extra3[,c("FID","IID","SexChildrenPerm2")],paste0(dir,"ukb_phenotypes_v3/indiv_filtered_v3_SexChildrenPerm2_cluster.txt"),quote=F,row.names=F)

##BOLT-LMM
#Phenotype file, excluding any individual with NA fitness, for input in BOLT-LMM 
# (white-space delimited, with only relevant phenotypes for logistic regression analysis, i.e. Sex=f.31.0.0.numeric, Assessment centre=f.54.0.0, Age=f.21003.0.0, PermutedSex=SexPerm2.numeric, + 20 PCs)
indiv.extra3$SexPerm2.numeric <- ifelse(indiv.extra3$SexPerm2=="Female",2,1)
indiv.extra3.bolt2 <- merge(indiv.extra3,pcs,all.x=T,by=c("FID","IID"))
#write.table(indiv.extra3.bolt2[,c("FID","IID","f.31.0.0.numeric","f.54.0.0","f.21003.0.0","SexPerm2.numeric","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")],paste0(dir,"ukb_phenotypes_v3/pheno_filtered_v3_perm2_for_boltlmm.txt"),quote=F,row.names=F,sep=" ")





