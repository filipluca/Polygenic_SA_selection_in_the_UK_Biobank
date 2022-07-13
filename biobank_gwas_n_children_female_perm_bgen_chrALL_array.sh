#!/bin/bash

#SBATCH --job-name=bolt_gwas_nchilren_perm_female_bgen_chrALL
#SBATCH --array=1-22%4

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#SBATCH --account=rf45
#SBATCH --array=1-22%4

### format: dd-hh:mm:ss
#SBATCH --time=0-12:00:00
#SBATCH --mem=64000

#SBATCH --partition=comp

module load bolt-lmm/2.3.4

bolt \
--bfile=ukb52049_chrAUTO_sampleqc1_snpqc1b \
--phenoFile=pheno_filtered_v3_perm1_for_boltlmm_nchildren.txt \
--phenoCol=Number_of_children_permuted \
--remove=list_of_male_ids_v3.txt \
--covarFile=pheno_filtered_v3_for_boltlmm_nchildren.txt \
--covarCol=f.54.0.0 \
--covarMaxLevels=22 \
--qCovarCol=f.21003.0.0 \
--qCovarCol=PC{1:20} \
--lmm \
--LDscoresFile=LDSCORE.1000G_EUR.tab.gz \
--statsFile=stats_gwas_nchildren_female_perm_genotyped_autosomes_when_running_imputed_chr${SLURM_ARRAY_TASK_ID}.tab \
--numThreads=8 \
--modelSnps=ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_ldpruned_0.2.prune.in \
--bgenFile=/home/fruzicka/rf45_scratch/ukb_imp_chr${SLURM_ARRAY_TASK_ID}_v3_sampleqc1_snpqc1b_8bit.bgen \
--sampleFile=/home/fruzicka/rf45_scratch/ukb_imp_chr${SLURM_ARRAY_TASK_ID}_v3_sampleqc1_snpqc1b_8bit.sample \
--remove=bolt.in_plink_but_not_imputed.FID_IID.7.txt \
--remove=indiv_without_LRS_v3.txt \
--exclude=ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_below0.8_maf_below0.01_artefacts_retained_v3.ids \
--statsFileBgenSnps=stats_gwas_nchildren_female_perm_imputed_chr${SLURM_ARRAY_TASK_ID}.tab
