
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



#Plots####

#Import partitions####
partitions.adult.fst <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.adultfst.results"),h=T)
partitions.adult.fst$Metric <- "Adult Fst"
partitions.adult.fst$Type <- "Observed"
partitions.adult.fst$Type2 <- "Viability"
partitions.adult.fst$Type3 <- "Fst"
partitions.adult.fst.perm2 <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.adultfst.perm2.results"),h=T)
partitions.adult.fst.perm2$Metric <- "Adult Fst"
partitions.adult.fst.perm2$Type <- "Permuted"
partitions.adult.fst.perm2$Type2 <- "Viability"
partitions.adult.fst.perm2$Type3 <- "Fst"
partitions.reproductive.fst <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.reproductivefst.results"),h=T)
partitions.reproductive.fst$Metric <- "Reproductive Fst"
partitions.reproductive.fst$Type <- "Observed"
partitions.reproductive.fst$Type2 <- "Reproduction"
partitions.reproductive.fst$Type3 <- "Fst"
partitions.reproductive.fst.perm <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.reproductivefst.perm.results"),h=T)
partitions.reproductive.fst.perm$Metric <- "Reproductive Fst"
partitions.reproductive.fst.perm$Type <- "Permuted"
partitions.reproductive.fst.perm$Type2 <- "Reproduction"
partitions.reproductive.fst.perm$Type3 <- "Fst"
partitions.gametic.fst <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.gameticfst.results"),h=T)
partitions.gametic.fst$Metric <- "Gametic Fst"
partitions.gametic.fst$Type <- "Observed"
partitions.gametic.fst$Type2 <- "Total"
partitions.gametic.fst$Type3 <- "Fst"
partitions.gametic.fst.perm2 <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.gameticfst.perm2.results"),h=T)
partitions.gametic.fst.perm2$Metric <- "Gametic Fst"
partitions.gametic.fst.perm2$Type <- "Permuted"
partitions.gametic.fst.perm2$Type2 <- "Total"
partitions.gametic.fst.perm2$Type3 <- "Fst"
partitions.unfolded.fst.pos <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstpos.results"),h=T)
partitions.unfolded.fst.pos$Metric <- "Unfolded Fst [Fst>0]"
partitions.unfolded.fst.pos$Type <- "Observed"
partitions.unfolded.fst.pos$Type2 <- "Reproduction 2"
partitions.unfolded.fst.pos$Type3 <- "Fst"
partitions.unfolded.fst.pos.perm <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstpos.perm.results"),h=T)
partitions.unfolded.fst.pos.perm$Metric <- "Unfolded Fst [Fst>0]"
partitions.unfolded.fst.pos.perm$Type <- "Permuted"
partitions.unfolded.fst.pos.perm$Type2 <- "Reproduction 2"
partitions.unfolded.fst.pos.perm$Type3 <- "Fst"
partitions.unfolded.fst.neg <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstneg.results"),h=T)
partitions.unfolded.fst.neg$Metric <- "Unfolded Fst [Fst<0]"
partitions.unfolded.fst.neg$Type <- "Observed"
partitions.unfolded.fst.neg$Type2 <- "Reproduction 2"
partitions.unfolded.fst.neg$Type3 <- "Fst"
partitions.unfolded.fst.neg.perm <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstneg.perm.results"),h=T)
partitions.unfolded.fst.neg.perm$Metric <- "Unfolded Fst [Fst<0]"
partitions.unfolded.fst.neg.perm$Type <- "Permuted"
partitions.unfolded.fst.neg.perm$Type2 <- "Reproduction 2"
partitions.unfolded.fst.neg.perm$Type3 <- "Fst"
partitions.lst <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.lst.results"),h=T)
partitions.lst$Metric <- "Lst"
partitions.lst$Type <- "Observed"
partitions.lst$Type2 <- "Viability"
partitions.lst$Type3 <- "non-Fst"
partitions.lst.perm2 <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.lst.perm2.results"),h=T)
partitions.lst.perm2$Metric <- "Lst"
partitions.lst.perm2$Type <- "Permuted"
partitions.lst.perm2$Type2 <- "Viability"
partitions.lst.perm2$Type3 <- "non-Fst"
partitions.t <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.t.results"),h=T)
partitions.t$Metric <- "|t|"
partitions.t$Type <- "Observed"
partitions.t$Type2 <- "Reproduction"
partitions.t$Type3 <- "non-Fst"
partitions.t.perm <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.t.perm.results"),h=T)
partitions.t.perm$Metric <- "|t|"
partitions.t.perm$Type <- "Permuted"
partitions.t.perm$Type2 <- "Reproduction"
partitions.t.perm$Type3 <- "non-Fst"
partitions.unfolded.t.pos <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtpos.results"),h=T)
partitions.unfolded.t.pos$Metric <- "Unfolded t [t>0]"
partitions.unfolded.t.pos$Type <- "Observed"
partitions.unfolded.t.pos$Type2 <- "Reproduction 2"
partitions.unfolded.t.pos$Type3 <- "non-Fst"
partitions.unfolded.t.pos.perm <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtpos.perm.results"),h=T)
partitions.unfolded.t.pos.perm$Metric <- "Unfolded t [t>0]"
partitions.unfolded.t.pos.perm$Type <- "Permuted"
partitions.unfolded.t.pos.perm$Type2 <- "Reproduction 2"
partitions.unfolded.t.pos.perm$Type3 <- "non-Fst"
partitions.unfolded.t.neg <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtneg.results"),h=T)
partitions.unfolded.t.neg$Metric <- "Unfolded t [t<0]"
partitions.unfolded.t.neg$Type <- "Observed"
partitions.unfolded.t.neg$Type2 <- "Reproduction 2"
partitions.unfolded.t.neg$Type3 <- "non-Fst"
partitions.unfolded.t.neg.perm <- read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtneg.perm.results"),h=T)
partitions.unfolded.t.neg.perm$Metric <- "Unfolded t [t<0]"
partitions.unfolded.t.neg.perm$Type <- "Permuted"
partitions.unfolded.t.neg.perm$Type2 <- "Reproduction 2"
partitions.unfolded.t.neg.perm$Type3 <- "non-Fst"
dd.partitions <- rbind(partitions.adult.fst,partitions.adult.fst.perm2,partitions.reproductive.fst,partitions.reproductive.fst.perm,partitions.gametic.fst,partitions.gametic.fst.perm2,partitions.unfolded.fst.pos,partitions.unfolded.fst.pos.perm,partitions.unfolded.fst.neg,partitions.unfolded.fst.neg.perm,partitions.lst,partitions.lst.perm2,partitions.t,partitions.t.perm,partitions.unfolded.t.pos,partitions.unfolded.t.pos.perm,partitions.unfolded.t.neg,partitions.unfolded.t.neg.perm)

#Pick only relevant annotations
dd.partitions.s <- subset(dd.partitions,Category %in% c("Coding_UCSCL2_0","Intron_UCSCL2_0","UTR_3_UCSCL2_0","UTR_5_UCSCL2_0"))
#Clean up
dd.partitions.s$Category <- as.factor(as.character(dd.partitions.s$Category))
levels(dd.partitions.s$Category) <- c("Coding","Intron","3' UTR","5' UTR")
dd.partitions.s$Category <- factor(dd.partitions.s$Category,levels = c("Coding","3' UTR","5' UTR","Intron"))
dd.partitions.s$Metric <- as.factor(as.character(dd.partitions.s$Metric))
dd.partitions.s$Metric <- factor(dd.partitions.s$Metric,levels = c("Adult Fst", "Lst", "Reproductive Fst", "|t|", "Gametic Fst", "Unfolded Fst [Fst<0]","Unfolded Fst [Fst>0]", "Unfolded t [t<0]",  "Unfolded t [t>0]"))
levels(dd.partitions.s$Metric) <- c("Adult Fst", "Lst", "Reproductive Fst", "|t|", "Gametic Fst", "Unfolded Fst (Negative)", "Unfolded Fst (Positive)", "Unfolded t (Negative)", "Unfolded t (Positive)")
#Q values
dd.partitions.s$Enrichment_q <- NA
dd.partitions.s$Enrichment_q[dd.partitions.s$Category!="Intron" &dd.partitions.s$Type=="Observed"] <- p.adjust(dd.partitions.s$Enrichment_p[dd.partitions.s$Category!="Intron"&dd.partitions.s$Type=="Observed"],"BH") 

#Plot partitions####
ggplot(subset(dd.partitions.s,Type=="Observed" & Category!="Intron"),aes(x=Metric,y=Enrichment,col=Type2))+
  theme_classic()+
  geom_point(position=position_dodge(0.5),aes(shape=Type3,size=Type3))+
  scale_color_manual(values=c("Viability" = "darkorange", "Reproduction" = "forestgreen", "Total" = "purple", "Reproduction 2" = "#40B0A6"))+
  scale_shape_manual(values=c("Fst" = 19, "non-Fst" = 18))+
  scale_size_manual(values=c("Fst" = 3, "non-Fst" = 4))+
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,position=position_dodge(0.5),size=1)+
  theme(axis.title.x= element_blank(),axis.title.y= element_text(size=20),axis.text = element_text(size=15),legend.position="none",strip.background=element_blank(),strip.text=element_text(size=20),axis.text.x = element_text(angle = 45, hjust=1))+
  geom_hline(yintercept = 1,linetype="dashed")+
  facet_wrap(~Category,ncol = 1,scales="free_y")+
  coord_cartesian(ylim=c(-10,28))

  
#Import heritability####
  heritability.adult.fst <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.adultfst.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.adultfst.log"),skip=28,fill=T)[1,6]),"Adult Fst","Observed", "Viability", "Fst")
  heritability.adult.fst.perm2 <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.adultfst.perm2.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.adultfst.perm2.log"),skip=28,fill=T)[1,6]), "Adult Fst","Permuted", "Viability", "Fst")
  heritability.reproductive.fst <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.reproductivefst.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.reproductivefst.log"),skip=28,fill=T)[1,6]), "Reproductive Fst","Observed", "Reproduction", "Fst")
  heritability.reproductive.fst.perm <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.reproductivefst.perm.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.reproductivefst.perm.log"),skip=28,fill=T)[1,6]), "Reproductive Fst","Permuted", "Reproduction", "Fst")
  heritability.gametic.fst <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.gameticfst.log"),skip=28,fill=T)[1,5]), as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.gameticfst.log"),skip=28,fill=T)[1,6]),"Gametic Fst","Observed", "Total", "Fst")
  heritability.gametic.fst.perm2 <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.gameticfst.perm2.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.gameticfst.perm2.log"),skip=28,fill=T)[1,6]),"Gametic Fst","Permuted", "Total", "Fst")
  heritability.unfolded.fst.pos <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstpos.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstpos.log"),skip=28,fill=T)[1,6]), "Unfolded Fst","Observed", "Reproduction 2", "Fst")
  heritability.unfolded.fst.pos.perm <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstpos.perm.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstpos.perm.log"),skip=28,fill=T)[1,6]), "Unfolded Fst","Permuted", "Reproduction 2", "Fst")
  heritability.unfolded.fst.neg <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstneg.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstneg.log"),skip=28,fill=T)[1,6]),"Unfolded Fst","Observed", "Reproduction 2", "Fst")
  heritability.unfolded.fst.neg.perm <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstneg.perm.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedfstneg.perm.log"),skip=28,fill=T)[1,6]), "Unfolded Fst","Permuted", "Reproduction 2", "Fst")
  heritability.lst <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.lst.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.lst.log"),skip=28,fill=T)[1,6]), "Lst","Observed", "Viability", "non-Fst")
  heritability.lst.perm2 <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.lst.perm2.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.lst.perm2.log"),skip=28,fill=T)[1,6]), "Lst","Permuted", "Viability", "non-Fst")
  heritability.t <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.t.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.t.log"),skip=28,fill=T)[1,6]), "|t|","Observed", "Reproduction", "non-Fst")
  heritability.t.perm <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.t.perm.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.t.perm.log"),skip=28,fill=T)[1,6]),"|t|","Permuted", "Reproduction", "non-Fst")
  heritability.unfolded.t.pos <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtpos.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtpos.log"),skip=28,fill=T)[1,6]), "Unfolded t","Observed", "Reproduction 2", "non-Fst")
  heritability.unfolded.t.pos.perm <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtpos.perm.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtpos.perm.log"),skip=28,fill=T)[1,6]), "Unfolded t","Permuted", "Reproduction 2", "non-Fst")
  heritability.unfolded.t.neg <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtneg.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtneg.log"),skip=28,fill=T)[1,6]), "Unfolded t","Observed", "Reproduction 2", "non-Fst")
  heritability.unfolded.t.neg.perm <- c(as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtneg.perm.log"),skip=28,fill=T)[1,5]),as.character(read.table(paste0(dir,"/ukb_LDSC_v3/partitions.unfoldedtneg.perm.log"),skip=28,fill=T)[1,6]), "Unfolded t","Permuted", "Reproduction 2", "non-Fst")
  
  dd.heritability <- data.frame(rbind(heritability.adult.fst,heritability.adult.fst.perm2,heritability.reproductive.fst,heritability.reproductive.fst.perm,heritability.gametic.fst,heritability.gametic.fst.perm2,heritability.unfolded.fst.pos,heritability.unfolded.fst.pos.perm,heritability.unfolded.fst.neg,heritability.unfolded.fst.neg.perm,heritability.lst,heritability.lst.perm2,heritability.t,heritability.t.perm,heritability.unfolded.t.pos,heritability.unfolded.t.pos.perm,heritability.unfolded.t.neg,heritability.unfolded.t.neg.perm))
  names(dd.heritability) <- c("Estimate","SE","Metric","Type","Component","Type2")
  dd.heritability$SE <- as.numeric(substr(dd.heritability$SE,2,6))
  dd.heritability$Estimate <- as.numeric(as.character(dd.heritability$Estimate))
  dd.heritability$Metric <- factor(dd.heritability$Metric,levels = c("Adult Fst","Lst","Reproductive Fst","|t|","Gametic Fst","Unfolded Fst","Unfolded t"))
  
ggplot(subset(dd.heritability,Component!="Reproduction 2"),aes(x=Metric,y=Estimate,col=Component))+
    theme_classic()+
    geom_point(position=position_dodge(0.5),aes(shape=Type2,size=Type2))+
    scale_color_manual(values=c("Viability" = "darkorange", "Reproduction" = "forestgreen", "Total" = "purple", "Reproduction 2" = "#40B0A6"))+
    scale_shape_manual(values=c("Fst" = 19, "non-Fst" = 18))+
    scale_size_manual(values=c("Fst" = 3, "non-Fst" = 4))+
    geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=.2,position=position_dodge(0.5),size=1)+
    theme(axis.title.x= element_blank(),axis.title.y=element_text(size=20),axis.text = element_text(size=15),legend.position="none",strip.background=element_blank(),strip.text=element_text(size=20),axis.text.x = element_text(angle = 45, hjust=1))+
      facet_wrap(~Type)+
  ylab("Estimated SNP-heritability (LDSC)")+
    geom_hline(yintercept = 0,linetype="dashed")
 
