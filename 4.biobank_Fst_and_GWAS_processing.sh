##################################################
##################################################
## Pm, Pf, Pm', Pf' (& Adult & Gametic Fst), observed
##################################################
##################################################

## Genotyped data 

#Autosomes
##Genotype counts
for i in {1..22}; do ~/Downloads/plink2 --bfile ukb52049_chr${i}_sampleqc1_snpqc1b --pheno ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v3_SexChildren_cluster.txt --loop-cats SexChildren --freq --geno-counts --hardy --missing --out ~/Dropbox/dropbox_work/data/biobank/ukb_allele_freq_v3/split_by_SexChildren/ukb52049_chr${i}_sampleqc1_snpqc1b; done
# See R code for Fst calculations
## Concatenate all autosomal Fst genotypes (make sure X chromosome and others are in a different folder)
cat ukb52049_chr*_sampleqc1_snpqc1b_no_further_maf_filter.Fst > ukb52049_chrAUTO_sampleqc1_snpqc1b_no_further_maf_filter.Fst
#See R code for artefact filtering


## Imputed data

#Autosomes
##Genotype counts
for i in {1..22}; do ~/Downloads/plink2 --bgen ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.bgen ref-first --sample ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.sample --pheno ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v3_SexChildren_cluster.txt --loop-cats SexChildren --geno-counts --out ~/Dropbox/dropbox_work/data/biobank/ukb_allele_freq_v3/split_by_SexChildren/ukb_imp_chr${i}_v3_sampleqc1_snpqc1b; done
# See R code for Fst calculations
# See R code for removal of low-info sites
## Filter autosomal SNPs with low info scores (<0.8)
for i in {1..22}; do awk -F' ' 'NR==FNR{c[$1]++;next};c[$2] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_imp_mfi/ukb_imp_chr${i}_v3_info_above_0.8.ids ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_no_info_filter_no_further_maf_filter.Fst > ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_info_above0.8_no_further_maf_filter.Fst; done
## Concatenate autosomal SNPs (make sure X-linked SNPs are in different folder)
cat ukb_imp_chr*_v3_sampleqc1_snpqc1b_info_above0.8_no_further_maf_filter.Fst > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_no_further_maf_filter.Fst
#See R code for MAF filtering


##################################################
##################################################
## Pm' and Pf' (& Gametic Fst), permuted offspring (i.e., Perm / Perm1)
##################################################
##################################################

## Genotype data

#Autosomes
#Genotype counts
for i in {1..22}; do ~/Downloads/plink2 --bfile ukb52049_chr${i}_sampleqc1_snpqc1b --pheno ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v3_SexChildrenPerm_cluster.txt --loop-cats SexChildrenPerm --freq --geno-counts --missing --out ~/Dropbox/dropbox_work/data/biobank/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb52049_chr${i}_sampleqc1_snpqc1b_perm1; done
## Concatenate chromosomes (make sure X chromosome and others are in a different folder)
cat ukb52049_chr*_sampleqc1_snpqc1b_perm1_no_further_maf_filter.Fst > ukb52049_chrAUTO_sampleqc1_snpqc1b_perm1_no_further_maf_filter.Fst
## Filter autosomal SNPs with MAF<0.01 and remove artefacts 
awk -F' ' 'NR==FNR{c[$1]++;next};c[$2] > 0' ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids ukb52049_chrAUTO_sampleqc1_snpqc1b_perm1_no_further_maf_filter.Fst > ukb52049_chrAUTO_sampleqc1_snpqc1b_perm1_maf_above0.01_artefacts_removed_v3.Fst


## Imputed data

#Autosomes
#Genotype counts
for i in {1..22}; do ~/Downloads/plink2 --bgen ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.bgen ref-first --sample ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.sample --pheno ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v3_SexChildrenPerm_cluster.txt --loop-cats SexChildrenPerm --geno-counts --out ~/Dropbox/dropbox_work/data/biobank/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_perm1; done
## Concatenate chromosomes (make sure X-linked SNPs are in different folder)
cat ukb_imp_chr*_v3_sampleqc1_snpqc1b_perm1_no_info_filter_no_further_maf_filter.Fst > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm1_no_info_filter_no_further_maf_filter.Fst
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$2] > 0' ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm1_no_info_filter_no_further_maf_filter.Fst > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm1_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst


#Imputed data, filtered, 100 permutations
for i in {1..22}; do ~/Downloads/plink2 --bgen ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.bgen ref-first --sample ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.sample --extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids --pheno ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v3_SexChildrenPerm_cluster.txt --loop-cats SexChildrenPerm --geno-counts --out ~/Dropbox/dropbox_work/data/biobank/ukb_allele_freqperm_v3/split_by_SexChildrenPerm/ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_perm1; done


##################################################
##################################################
## Pm' and Pf' (& Gametic Fst), permuted sexes (Perm2)
##################################################
##################################################

## Genotype data

#Autosomes
#Genotype counts
for i in {1..22}; do ~/Downloads/plink2 --bfile ukb52049_chr${i}_sampleqc1_snpqc1b --pheno ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v3_SexChildrenPerm2_cluster.txt --loop-cats SexChildrenPerm2 --geno-counts --out ~/Dropbox/dropbox_work/data/biobank/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb52049_chr${i}_sampleqc1_snpqc1b_perm2; done
## Concatenate all autosomal Fst genotypes (make sure X chromosome and others are in a different folder)
cat ukb52049_chr*_sampleqc1_snpqc1b_perm2_no_further_maf_filter.Fst > ukb52049_chrAUTO_sampleqc1_snpqc1b_perm2_no_further_maf_filter.Fst
## Filter autosomal SNPs with MAF<0.01 and remove artefacts 
awk -F' ' 'NR==FNR{c[$1]++;next};c[$2] > 0' ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids ukb52049_chrAUTO_sampleqc1_snpqc1b_perm2_no_further_maf_filter.Fst > ukb52049_chrAUTO_sampleqc1_snpqc1b_perm2_maf_above0.01_artefacts_removed_v3.Fst


## Imputed data

#Autosomes
#Genotype counts
for i in {1..22}; do ~/Downloads/plink2 --bgen ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.bgen ref-first --sample ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_8bit.sample --pheno ~/Dropbox/dropbox_work/data/biobank/ukb_phenotypes/indiv_filtered_v3_SexChildrenPerm2_cluster.txt --loop-cats SexChildrenPerm2 --geno-counts --out ~/Dropbox/dropbox_work/data/biobank/ukb_allele_freqperm_v3/split_by_SexChildrenPerm2/ukb_imp_chr${i}_v3_sampleqc1_snpqc1b_perm2; done
## Concatenate autosomal SNPs (make sure X-linked SNPs are in different folder)
cat ukb_imp_chr*_v3_sampleqc1_snpqc1b_perm2_no_info_filter_no_further_maf_filter.Fst > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm2_no_info_filter_no_further_maf_filter.Fst
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$2] > 0' ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm2_no_info_filter_no_further_maf_filter.Fst > ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_perm2_info_above0.8_maf_above0.01_artefacts_removed_v3.Fst;




##################################################
##################################################
## Gwas of sex
##################################################
##################################################

## Observed
## Concatenate GWAS output across chromosomes
head -1 stats_gwas_mf_imputed_chr1.tab > header.gwas_mf;
awk 'FNR>1' stats_gwas_mf_imputed_chr*.tab > stats_gwas_mf_imputed_chrAUTO.tmp.tab;
cat header.gwas_mf stats_gwas_mf_imputed_chrAUTO.tmp.tab > stats_gwas_mf_imputed_chrAUTO.tab;
rm *tmp*;
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids stats_gwas_mf_imputed_chrAUTO.tab > stats_gwas_mf_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab;

##Permuted
## Concatenate GWAS output across chromosomes
awk 'FNR>1' stats_gwas_mf_perm2_imputed_chr*.tab > stats_gwas_mf_perm2_imputed_chrAUTO.tmp.tab;
cat header.gwas_mf stats_gwas_mf_perm2_imputed_chrAUTO.tmp.tab > stats_gwas_mf_perm2_imputed_chrAUTO.tab;
rm *tmp*;
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids stats_gwas_mf_perm2_imputed_chrAUTO.tab > stats_gwas_mf_perm2_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab;


##################################################
##################################################
## Gwas of LRS
##################################################
##################################################

## Females

## Observed
## Concatenate GWAS output across chromosomes 
awk 'FNR>1' stats_gwas_nchildren_female_imputed_chr*.tab > stats_gwas_nchildren_female_imputed_chrAUTO.tab; 
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids stats_gwas_nchildren_female_imputed_chrAUTO.tab > stats_gwas_nchildren_female_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab;

##Permuted
## Concatenate GWAS output across chromosomes 
awk 'FNR>1' stats_gwas_nchildren_female_perm_imputed_chr*.tab > stats_gwas_nchildren_female_perm_imputed_chrAUTO.tab; 
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids stats_gwas_nchildren_female_perm_imputed_chrAUTO.tab > stats_gwas_nchildren_female_perm_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab;

## Males
## Observed
## Concatenate GWAS output across chromosomes 
awk 'FNR>1' stats_gwas_nchildren_male_imputed_chr*.tab > stats_gwas_nchildren_male_imputed_chrAUTO.tab; 
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids stats_gwas_nchildren_male_imputed_chrAUTO.tab > stats_gwas_nchildren_male_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab;

##Permuted
## Concatenate GWAS output across chromosomes 
awk 'FNR>1' stats_gwas_nchildren_male_perm_imputed_chr*.tab > stats_gwas_nchildren_male_perm_imputed_chrAUTO.tab; 
## Filter out autosomal SNPs with low info scores (<0.8), MAF<0.01 and artefacts
awk -F' ' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids stats_gwas_nchildren_male_perm_imputed_chrAUTO.tab > stats_gwas_nchildren_male_perm_imputed_chrAUTO_info_above0.8_maf_above0.01_artefacts_removed_v3.tab;

