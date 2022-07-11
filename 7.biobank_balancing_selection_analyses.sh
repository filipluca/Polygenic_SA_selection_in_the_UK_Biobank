############################################
## Allele ages
############################################

## Concatenate allele age files
awk 'FNR>4' atlas.chr*.csv > atlas.chrAUTO.csv;
#head -4 atlas.chr1.csv > header.csv;

## Extract filtered SNP IDs

#Imputed
awk -F',' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids atlas.chrAUTO.csv > atlas.ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.csv
cat header.csv atlas.ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.csv > atlas.ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3_with_header.csv
rm atlas.ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.csv
#Genotyped
awk -F',' 'NR==FNR{c[$1]++;next};c[$1] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids atlas.chrAUTO.csv > atlas.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.csv
cat header.csv atlas.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.csv > atlas.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3_with_header.csv
rm atlas.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.csv


############################################
## 1000 genomes, allele frequencies of each subpopulation
############################################

## Estimate allele frequencies
for i in {1..22}; do ~/Downloads/bcftools/bcftools view -Ou -S sample_ids_YRI.txt ../../X_chromosome_purifying_selection/1000_genomes/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | ~/Downloads/bcftools/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN \n' > ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.YRI; done;
for i in {1..22}; do ~/Downloads/bcftools/bcftools view -Ou -S sample_ids_GIH.txt ../../X_chromosome_purifying_selection/1000_genomes/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | ~/Downloads/bcftools/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN \n' > ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.GIH; done;
#for i in {1..22}; do ~/Downloads/bcftools/bcftools view -Ou -S sample_ids_TSI.txt ../../X_chromosome_purifying_selection/1000_genomes/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | ~/Downloads/bcftools/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN \n' > ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.TSI; done;
#for i in {1..22}; do ~/Downloads/bcftools/bcftools view -Ou -S sample_ids_PUR.txt ../../X_chromosome_purifying_selection/1000_genomes/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | ~/Downloads/bcftools/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN \n' > ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.PUR; done;
#for i in {1..22}; do ~/Downloads/bcftools/bcftools view -Ou -S sample_ids_CHS.txt ../../X_chromosome_purifying_selection/1000_genomes/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | ~/Downloads/bcftools/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN \n' > ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.CHS; done;

##Concatenate chromosomes
cat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.YRI > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.YRI
cat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.GIH > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.GIH
#cat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.TSI > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.TSI
#cat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.PUR > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.PUR
#cat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.CHS > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.CHS

## Extract filtered SNP IDs
#Imputed
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.YRI > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.YRI;
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.GIH > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.GIH;
#awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.TSI > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.TSI;
#awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.PUR > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.PUR;
#awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb_imp_chrAUTO_v3_sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.CHS > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.sampleqc1_snpqc1b_info_above0.8_maf_above0.01_artefacts_removed_v3.CHS;


#Genotyped
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.YRI > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.YRI;
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.GIH > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.GIH;
#awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.TSI > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.TSI;
#awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.PUR > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.PUR;
#awk -F'\t' 'NR==FNR{c[$1]++;next};c[$3] > 0' ~/Dropbox/dropbox_work/data/biobank/ukb_fst_v3/ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.ids ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.CHS > ALL.chrAUTO.phase3_shapeit2_mvncall_integrated_v5a.20130502.ukb52049_chrAUTO_sampleqc1_snpqc1b_maf_above0.01_artefacts_removed_v3.CHS;


############################################
## 1000 genomes, Tajima's D
############################################

parallel -J 3 'vcftools --keep sample_ids_YRI.txt --max-missing 0.95 --gzvcf ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --TajimaD 10000 --out chr{.}.YRI' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22;
parallel -J 3 'vcftools --keep sample_ids_GIH.txt --max-missing 0.95 --gzvcf ALL.chr{.}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --TajimaD 10000 --out chr{.}.GIH' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22;

cat chr*.YRI.Tajima.D > chrAUTO.YRI.Tajima.D
cat chr*.GIH.Tajima.D > chrAUTO.GIH.Tajima.D
