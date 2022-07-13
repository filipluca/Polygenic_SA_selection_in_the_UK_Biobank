#!/bin/bash

#SBATCH --job-name=plink2_ldprune_bgens_chrAUTO

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#SBATCH --account=rf45

### format: dd-hh:mm:ss
#SBATCH --time=01-06:00:00
#SBATCH --mem=128000

#SBATCH --partition=comp

module load plink/2.0-alpha

plink2 \
--bgen ukb_imp_chr1_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr1_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr1_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr2_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr2_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr2_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr3_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr3_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr3_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr4_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr4_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr4_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr5_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr5_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr5_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr6_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr6_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr6_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr7_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr7_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr7_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr8_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr8_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr8_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr9_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr9_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr9_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr10_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr10_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr10_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr11_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr11_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr11_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr12_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr12_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr12_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr13_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr13_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr13_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr14_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr14_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr14_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr15_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr15_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr15_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr16_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr16_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr16_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr17_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr17_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr17_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr18_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr18_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr18_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr19_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr19_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr19_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr20_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr20_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr20_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr21_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr21_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr21_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2

plink2 \
--bgen ukb_imp_chr22_v3_sampleqc1_snpqc1b_8bit.bgen \
--sample ukb_imp_chr22_v3_sampleqc1_snpqc1b_8bit.sample \
--indep-pairwise 50 10 0.2 \
--extract ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids \
--remove indiv_without_LRS_v3.txt \
--out ukb_imp_chr22_v3_sampleqc1_snpqc1b_8bit_ldpruned_0.2
