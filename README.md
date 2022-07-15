# Polygenic signals of sex differences in selection in humans from the UK Biobank


Welcome!

This page contains the code underlying analyses and figures presented in 'Polygenic signals of sex differences in selection in contemporary humans' (Ruzicka et al. 2022, bioRxiv). 

The code is organised into text files 1-8, which approximately follow the order of the analyses presented in the manuscript. I used a mixture of R and bash scripts, so it is usually necessary to go back and forth between these two sets of scripts to replicate a given analysis. The contents of each file are described in bullet-point form below. I have also annotated the code for clarity, but if it's insufficiently clear please e-mail me at filip.ruzicka [at] monash.edu. See also https://github.com/lukeholman/UKBB_LDSC for some further code/output pertaining to functional analyses in the manuscript

The 8th text file contains more detailed descriptions of the datafiles shared on zenodo (https://zenodo.org/record/6824671). These zenodo datafiles contain  metadata used to produce the figures in the manuscript, and will be the most useful objects for future users.


#### 1. Downloading UK Biobank data

-Download involves .bim and .bgen files

#### 2. Genotype quality control

-Implement basic sample-level and site-level quality controls (see relevant code from file #3)

-LD prune and extract top PCs (incl specific scripts for running the code on the cluster)

#### 3. Phenotypic filtering

-Basic sample-level quality controls for relatedness, aneuploidy, etc.

-Filter out individuals with missing and/or unreliable fitness data

#### 4. Fst and GWAS processing

-Calculate sex-specific allele frequencies, among adults and projected gametes, among observed and permuted data

-Filter out sites with MAF<0.01 and potential artefacts

-Calculate Fst

-Run GWAS analysis in BOLT-LMM (incl specific scripts for running the code on the cluster)

#### 5. Fst and GWAS data analysis

-Import Fst and GWAS results in R 

-Summary statistics and statistical tests for overall enrichment

-P-and Q-values

-Associations between Fst/GWAS metrics and MAF

-Plot results

#### 6. Functional analyses

-Prepare files for LD score regression functional enrichment analysis (i.e., make sure alleles are compatible with HapMap, polarise Fst and mixed-model statistics)

-Run Stratified LDSC to examine functional enrichments

-Plot functional enrichments

#### 7. Balancing selection analyses

-Process data on allele ages, between-population Fst, and previous candidates for balancing selection

-Run associations between each metric of balancing selection and metrics of sex-differential selection

-Plot results

#### 8. Data files

-Description of summary data files uploaded to zenodo

#### Supplementary

-Code used to generate Supplementary Figures (section D), on the covariance between Fst and MAF
