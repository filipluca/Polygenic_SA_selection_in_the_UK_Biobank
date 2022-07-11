#######################
## SNPEff #############
#######################

## Obtain SNPEff annotations for Biobank data 
#UK Biobank uses GRCh37/hg19 

## Autosomes

## Genotyped data

#Genomic positions for genotyped polymorphic sites
cat ukb52049_chr*.bim >  ~/Dropbox/dropbox_work/data/biobank/ukb_annotations_v3/ukb52049_chrAUTO.bim 

##CHROM, ID and POS from MAF/artefact-filtered genotyped SNPs
cut -f1,2,3,4 ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.Fst > ~/Dropbox/dropbox_work/data/biobank/ukb_annotations_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.4columns
#See R code for merging MAF-unfiltered .bim file with MAF/artefact-filtered .Fst file
#Add hash to first line of text file to make it a vcf

#Run SNPEff for the 1st time
java -Xmx6g -jar ~/Downloads/snpEff_latest_core/snpEff/snpEff.jar hg19 ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.vcf > ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ann.vcf
#Pick out sites that do not match reference genome, and invert reference and alternative allele order
awk -F'\t' '/WARNING_REF_DOES_NOT_MATCH_GENOME/ {t = $4; $4 = $5; $5 = t; print;}' OFS=$'\t' ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ann.vcf > ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.tmp.vcf;
#Combine these 'reference-flipped' sites with all sites (there should now be duplicate sites)
cat ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ann.vcf ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.tmp.vcf > ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.tmp2.vcf;
#Re-run SNPEff on all sites 
java -Xmx6g -jar ~/Downloads/snpEff_latest_core/snpEff/snpEff.jar hg19 ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.tmp2.vcf > ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.tmp2.ann.vcf;
#Remove erroneous sites
awk -F'\t' '!/WARNING_REF_DOES_NOT_MATCH_GENOME/' OFS=$'\t' ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.tmp2.ann.vcf > ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3_refcorrected.ann.vcf;
#Now Annotations should all be correct; extract SNP variant effects and gene names
java -Xmx6g -jar ~/Downloads/snpEff_latest_core/snpEff/SnpSift.jar extractFields ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3_refcorrected.ann.vcf "CHROM" "POS" "ID" "REF" "ALT" "ANN[0].EFFECT" "ANN[0].GENE" > ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3_refcorrected.ann.effect.gene.txt
#Remove temporary files
rm *tmp*; rm ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ann.vcf;



## Imputed data

#Genomic position for each polymorphic site (i.e. .bim file)
for i in {1..22}; do ~/Downloads/plink2 --bgen ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.bgen ref-first --sample ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.sample --make-just-bim --out ukb_imp_chr${i}_v3_sampleqc1_snpqc1b; done
cat ukb_imp_chr*_v3_sampleqc1_snpqc1b.bim > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b.bim

##CHROM, ID and POS from all imputed SNPs
cut -f1,2,3,4 ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst > ~/Dropbox/dropbox_work/data/biobank/ukb_annotations_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.4columns
#See R code for merging MAF-unfiltered .bim file with MAF-filtered .Fst file
#Add hash to first line of text file to make it a vcf


#Run SNPEff for the 1st time
java -Xmx12g -jar ~/Downloads/snpEff_latest_core/snpEff/snpEff.jar hg19 ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.vcf > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ann.vcf
#Pick out sites that do not match reference genome, and invert reference and alternative allele order
awk -F'\t' '/WARNING_REF_DOES_NOT_MATCH_GENOME/ {t = $4; $4 = $5; $5 = t; print;}' OFS=$'\t' ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ann.vcf > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.tmp.vcf;
#Combine these 'reference-flipped' sites with all sites (there should now be duplicate sites)
cat ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ann.vcf ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.tmp.vcf > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.tmp2.vcf;
#Re-run SNPEff on all sites 
java -Xmx12g -jar ~/Downloads/snpEff_latest_core/snpEff/snpEff.jar hg19 ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.tmp2.vcf > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.tmp2.ann.vcf;
#Remove erroneous sites
awk -F'\t' '!/WARNING_REF_DOES_NOT_MATCH_GENOME/' OFS=$'\t' ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.tmp2.ann.vcf > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_refcorrected.ann.vcf;
#Do any sites not match reference genome?
grep "WARNING_REF_DOES_NOT_MATCH_GENOME" ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_refcorrected.ann.vcf | wc -l
#Now Annotations should all be correct; extract SNP variant effects and gene names
java -Xmx12g -jar ~/Downloads/snpEff_latest_core/snpEff/SnpSift.jar extractFields ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_refcorrected.ann.vcf "CHROM" "POS" "ID" "REF" "ALT" "ANN[0].EFFECT" "ANN[0].GENE" > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_refcorrected.ann.effect.gene.txt
#Remove temporary files
rm *tmp*; rm ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ann.vcf;




#######################
## LDSC ###############
#######################


#Import LD score regression files

URL=https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE
#Download reference files
wget $URL/1000G_Phase3_frq.tgz
wget $URL/1000G_Phase3_plinkfiles.tgz
wget $URL/1000G_Phase3_weights_hm3_no_MHC.tgz
wget $URL/1000G_Phase3_ldscores.tgz
wget $URL/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
#Untar the files
tar zxvf 1000G_Phase3_frq.tgz
tar zxvf 1000G_Phase3_plinkfiles.tgz
tar zxvf 1000G_Phase3_weights_hm3_no_MHC.tgz
tar zxvf 1000G_Phase3_ldscores.tgz
tar zxvf 1000G_Phase3_baselineLD_v2.2_ldscores.tgz
rm *tgz
#Organize and rename directories
mv 1000G_EUR_Phase3_plink plink_files
mv 1000G_Phase3_frq/* plink_files; rm -rf 1000G_Phase3_frq
mv 1000G_Phase3_weights_hm3_no_MHC weights_hm3_no_MHC
mv LDscore ldscores
mkdir baselineLD_v2.2; mv baselineLD.* baselineLD_v2.2



##Identify compatible alleles between our data and w_hm3

#Header from .sumstats file
head -1 PASS_Adult_Fst.sumstats > header.sumstats
#Keep only relevant columns to generate a list of SNPs for our study
cut -f1,2,3 PASS_Adult_Fst.sumstats > ruzickaetal.snplist
#Make HapMap3 snplist space-delimited
tr -s '\t' ' ' < w_hm3.snplist > w_hm3.spacedelimited.snplist
#Generate HapMap3 snplist where A1 and A2 alleles have been flipped
awk 'BEGIN {FS=" "; OFS=" "} {print $1, $3, $2}' w_hm3.spacedelimited.snplist > w_hm3.flipped.spacedelimited.snplist

#Keep only compatible alleles (for original order of A1 and A2 in w_hm3)
comm -12 <(sort -k 1 ruzickaetal.snplist) <(sort -k 1 w_hm3.spacedelimited.snplist) > ruzickaetal.snplist.compatible.original
#Keep only compatible alleles (for flipped order of A1 and A2 in w_hm3)
comm -12 <(sort -k 1 ruzickaetal.snplist) <(sort -k 1 w_hm3.flipped.spacedelimited.snplist) > ruzickaetal.snplist.compatible.flipped
#Combine list of compatible SNPs
cat ruzickaetal.snplist.compatible.original ruzickaetal.snplist.compatible.flipped > ruzickaetal.snplist.compatible



## Calculate partitioned heritabilities using baseline model

#Adult Fst, observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Adult_Fst.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.adultfst;

#Adult Fst, permuted 
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Adult_Fst_perm2.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.adultfst.perm2;

#Reproductive Fst, observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Reproductive_Fst.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.reproductivefst;

#Reproductive Fst, permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Reproductive_Fst_perm.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.reproductivefst.perm;

#Gametic Fst, Observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Gametic_Fst.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.gameticfst;

# Gametic Fst,Permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Gametic_Fst_perm2.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.gameticfst.perm2;

#Unfolded Fst, Positive, Observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_Fst_positive.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedfstpos;

#Unfolded Fst, Positive, Permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_Fst_positive_perm.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedfstpos.perm;

#Unfolded Fst, Negative, Observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_Fst_negative.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedfstneg;

#Unfolded Fst, Negative, Permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_Fst_negative_perm.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedfstneg.perm;

#Lst, Observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Lst.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.lst;

#Lst Permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Lst_perm2.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.lst.perm2;

#t, Observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_t.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.t;

#t, Permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_t_perm.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.t.perm;

#Unfolded t, Positive, Observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_t_positive.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedtpos;

#Unfolded t, Positive, Permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_t_positive_perm.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedtpos.perm;

#Unfolded t, Negative, Observed
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_t_negative.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedtneg;

#Unfolded t, Negative, Permuted
~/Downloads/ldsc/ldsc.py --h2 sumstats/PASS_Unfolded_t_negative_perm.sumstats \
--ref-ld-chr baselineLD_v2.2/baselineLD. \
--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr plink_files/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--out partitions.unfoldedtneg.perm;


## Calculate heritability for all metrics
#~/Downloads/ldsc/ldsc.py \
#--h2 sumstats/PASS_Adult_Fst.sumstats \
#--ref-ld-chr ldscores/LDscore. \
#--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
#--out heritability.adultfst;


## Calculate genetic correlation between all traits
#~/Downloads/ldsc/ldsc.py \
#--rg sumstats/PASS_Unfolded_Fst_positive.sumstats,sumstats/PASS_unfolded_t_positive.sumstats \
#--ref-ld-chr ldscores/LDscore. \
#--w-ld-chr weights_hm3_no_MHC/weights.hm3_noMHC. \
#--out test.correlation



