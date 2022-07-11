##############################################
##############################################
## Sample/marker quality control -- genotype data
##############################################
##############################################

## Standardise prefix across bed/bim/fam files 

#Before standardising, prefix of files is e.g. ukb52049_cal_chr10_v2_s488282.fam, ukb_cal_chr10_v2.bed, ukb_snp_chr10_v2.bim
#After, it is ukb52049_chr1.fam/bed/bim, etc...

## First round of sample-level filtering

## See R code for details of sample-level filtering 

#Only keep individuals that meet filtering criteria
for i in {1..22}; do ~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chr${i} --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --make-bed --out ukb52049_chr${i}_sampleqc1; done
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrX --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --make-bed --out ukb52049_chrX_sampleqc1;
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrXY --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --make-bed --out ukb52049_chrXY_sampleqc1;
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrY --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --make-bed --out ukb52049_chrY_sampleqc1;
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrMT --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --make-bed --out ukb52049_chrMT_sampleqc1;

## Marker-level filtering
#Filter based on genotype call rate (>95%), MAF(>0.0001), HWE, remove indels
for i in {1..22}; do ~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chr${i}_sampleqc1 --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --make-bed --out ukb52049_chr${i}_sampleqc1_snpqc1b; rm ukb52049_chr${i}_sampleqc1.*; done
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrX_sampleqc1 --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only  --make-bed --out ukb52049_chrX_sampleqc1_snpqc1b; rm ukb52049_chrX_sampleqc1.*;
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrXY_sampleqc1 --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only  --make-bed --out ukb52049_chrXY_sampleqc1_snpqc1b; rm ukb52049_chrXY_sampleqc1.*;
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrY_sampleqc1 --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only  --make-bed --out ukb52049_chrY_sampleqc1_snpqc1b; rm ukb52049_chrY_sampleqc1.*;
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chrMT_sampleqc1 --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only  --make-bed --out ukb52049_chrMT_sampleqc1_snpqc1b; rm ukb52049_chrMT_sampleqc1.*;



##############################################
##############################################
## Sample/marker quality control -- imputed data
##############################################
##############################################

#Filter samples (population structure, non-White British, aneuploids, high missing rates, age<45) and snps (MAF, missingness, HWE)
~/Downloads/plink2 --bgen ukb_imp_chr22_v3.bgen ref-first --sample ukb52049_imp_chr22_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr22_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr21_v3.bgen ref-first --sample ukb52049_imp_chr21_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3  --out ukb_imp_chr21_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr20_v3.bgen ref-first --sample ukb52049_imp_chr20_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3  --out ukb_imp_chr20_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr19_v3.bgen ref-first --sample ukb52049_imp_chr19_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr19_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr18_v3.bgen ref-first --sample ukb52049_imp_chr18_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr18_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr17_v3.bgen ref-first --sample ukb52049_imp_chr17_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr17_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr16_v3.bgen ref-first --sample ukb52049_imp_chr16_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr16_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr15_v3.bgen ref-first --sample ukb52049_imp_chr15_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr15_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr14_v3.bgen ref-first --sample ukb52049_imp_chr14_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr14_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr13_v3.bgen ref-first --sample ukb52049_imp_chr13_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr13_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr12_v3.bgen ref-first --sample ukb52049_imp_chr12_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr12_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr11_v3.bgen ref-first --sample ukb52049_imp_chr11_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr11_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr10_v3.bgen ref-first --sample ukb52049_imp_chr10_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr10_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr9_v3.bgen ref-first --sample ukb52049_imp_chr9_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr9_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr8_v3.bgen ref-first --sample ukb52049_imp_chr8_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr8_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr7_v3.bgen ref-first --sample ukb52049_imp_chr7_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr7_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr6_v3.bgen ref-first --sample ukb52049_imp_chr6_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr6_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr5_v3.bgen ref-first --sample ukb52049_imp_chr5_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr5_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr4_v3.bgen ref-first --sample ukb52049_imp_chr4_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr4_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr3_v3.bgen ref-first --sample ukb52049_imp_chr3_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr3_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr2_v3.bgen ref-first --sample ukb52049_imp_chr2_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr2_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chr1_v3.bgen ref-first --sample ukb52049_imp_chr1_v3_s487297.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chr1_v3_sampleqc1_snpqc1b; 
#X chromosome
~/Downloads/plink2 --bgen ukb_imp_chrX_v3.bgen ref-first --sample ukb52049_imp_chrX_v3_s486646.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chrX_v3_sampleqc1_snpqc1b; 
~/Downloads/plink2 --bgen ukb_imp_chrXY_v3.bgen ref-first --sample ukb52049_imp_chrXY_v3_s486332.sample --keep ~/Documents/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.3 --out ukb_imp_chrXY_v3_sampleqc1_snpqc1b; 


## Make bgen files into 8-bit output, for more compact storage and usage in BOLT-LMM
~/Downloads/plink2 --bgen ukb_imp_chrXY_v3.bgen ref-first --sample ukb52049_imp_chrXY_v3_s486332.sample --keep ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.2 bits=8 --out ukb_imp_chrXY_v3_sampleqc1_snpqc1b_8bit; 
~/Downloads/plink2 --bgen ukb_imp_chr22_v3.bgen ref-first --sample ukb52049_imp_chr22_v3_s487297.sample --keep ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.2 bits=8 --out ukb_imp_chr22_v3_sampleqc1_snpqc1b_8bit; 
for i in {1..22}; do ~/Downloads/plink2 --bgen ukb_imp_chr${i}_v3.bgen ref-first --sample ukb52049_imp_chr${i}_v3_s487297.sample --keep ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v2.txt --maf 0.0001 --geno 0.05 --hwe 0.000001 --snps-only --export bgen-1.2 bits=8 --out ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit; done 


##############################################
##############################################
## LD pruning SNPs 
##############################################
##############################################

#Genotyped
#Concatenate all autosomal sites
~/Downloads/plink_mac_20190304/plink --bfile ukb52049_chr1_sampleqc1_snpqc1b --merge-list all_autosomal_bed_files.txt --make-bed --out ukb52049_chrAUTO_sampleqc1_snpqc1b 
#LD prune, excluding MAF<0.01 and artefacts
~/Downloads/plink2 --bfile ukb52049_chrAUTO_sampleqc1_snpqc1b --indep-pairwise 50 10 0.2 --extract ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids --remove ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes_v3/indiv_without_LRS_v3.txt --out ~/Dropbox/dropbox_work/data/biobank/ukb_ldpruning_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3_ldpruned

#Imputed
#LD prune, excluding MAF<0.01, INFO<0.8 and artefacts
#see plink2_ldprune_chrAUTO_v3.sh 
cat ukb_imp_chr*0.2.prune.in > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_ldpruned_0.2.prune.in
cat ukb_imp_chr*0.2.prune.out > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_ldpruned_0.2.prune.out


##############################################
##############################################
## Principal components
##############################################
##############################################

#see plink2_extract_top_20PCs.sh


##############################################
##############################################
## SNPs to remove when conducting GWAS
##############################################
##############################################



## File containing imputed SNPs with MAF<0.01, info score<0.8 and artefacts
# File containing all impiuted SNPs (no filtering at all)
cat ukb_imp_chr*_v3_sampleqc1_snpqc1b_no_info_filter_no_further_maf_filter.Fst > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_no_info_filter_no_further_maf_filter.Fst
cut -f2 ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_no_info_filter_no_further_maf_filter.Fst > cut -f2 ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_no_info_filter_no_further_maf_filter.ids
# File containing only imputed SNPs that have been filtered out (i.e., MAF<0.01, Info<0.8, etc)
comm -13 <(sort ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids) <(sort ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_no_info_filter_no_further_maf_filter.ids) > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_below0.8_maf_below0.01_artefacts_retained_v3.ids
