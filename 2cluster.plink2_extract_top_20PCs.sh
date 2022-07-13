#!/bin/bash

#SBATCH --job-name=plink2_extract_PCs

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#SBATCH --account=rf45

### format: dd-hh:mm:ss
#SBATCH --time=1-00:00:00
#SBATCH --mem=256000

#SBATCH --partition=comp

module load plink/2.0-alpha

plink2 \
--bfile ukb52049_chrAUTO_sampleqc1_snpqc1b \
--extract ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3_ldpruned.prune.in \
--remove indiv_without_LRS_v3.txt \
--pca approx 20 \
--out ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3


