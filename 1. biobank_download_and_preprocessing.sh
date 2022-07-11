####################################################
## BIOBANK DOWNLOAD ####
####################################################

## See guide for downloading Biobank datasets: http://biobank.ctsu.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.0.pdf

## Download helper programs: ukbdm5, ukbunpack, ukbconv, ukbfetch, ukblink, ukbgene, encoding.ukb
chmod 755 ukbmd5
chmod 755 ukbunpack
chmod 755 ukbconv
chmod 755 ukbfetch
chmod 755 ukblink
chmod 755 ukbgene

## Download main dataset (ukb37327.enc) using instructions from e-mail (Application ID 52049, ID 37323)

## Helper utilities are not MAC compatible, so do the following:
#update to Xcode 11, then
#brew install linux-noah/noah/noah
## Decrypt main data file
#noah ./ukbmd5 ukb37327.enc
#noah ./ukbunpack ukb37327.enc k52049.key
## Create data dictionary
#noah ./ukbconv ukb37327.enc_ukb docs
## Convert to R format
#noah ./ukbconv ukb37327.enc_ukb r
## Download BED file from chromosome 17
#noah ./ukbgene cal -c17 -ak52049.key 


#Guide on how to install Virtual Machine (Linux) on Mac 
#https://www.dev2qa.com/how-to-install-ubuntu-on-virtualbox-mac/
#(if installation fails, try assigning more ram or using non-LTS version)
#If menu bar disappears from Ubuntu, try ctrl+opt+t
#To increase VirtualBox HD size, navigate to folder with VirtualBox disk image (e.g. /Users/fruz0001/VirtualBox\ VMs/Ubuntu/), then
VBoxManage modifyhd YOUR_HARD_DISK.vdi --resize SIZE_IN_MB

##Within Biobank folder on Virtual Machine

## Decrypt main data file
./ukbmd5 ukb37327.enc
./ukbunpack ukb37327.enc k52049.key
## Create data dictionary
./ukbconv ukb37327.enc_ukb docs
## Convert to R format
./ukbconv ukb37327.enc_ukb r

##GENOTYPE DATA
#Download .fam files for all 26 chromosomes
#Make "chr_names.txt" file containing "-c1" on line 1, "-c2" on line 2, etc until line 26
for i in $(cat chr_names.txt); do ./ukbgene cal $i -m -ak52049.key; done 
##Download .bed files for each chromosome, a few chromosomes at a time (too large to loop across all chromosomes at once)
for i in $(cat chr_names_1_to_5.txt); do ./ukbgene cal $i -ak52049.key; done 
for i in $(cat chr_names_6_to_10.txt); do ./ukbgene cal $i -ak52049.key; done 
for i in $(cat chr_names_11_to_21.txt); do ./ukbgene cal $i -ak52049.key; done 
for i in $(cat chr_names_22_to_26.txt); do ./ukbgene cal $i -ak52049.key; done 

#Download relatedness data
./ukbgene rel -ak52049.keyls

#Download genetic data description
 wget  -nd  biobank.ndph.ox.ac.uk/showcase/showcase/docs/ukb_genetic_data_description.txt
 
#Download SNP Quality Control information
 wget  -nd  biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/ukb_snp_qc.txt

##IMPUTATION DATA

#Download '.bgen' files
#Navigate to shared folder in Virtual machine 
#see https://askubuntu.com/questions/161759/how-to-access-a-shared-folder-in-virtualbox

#Download '.sample' files
for i in {1..26}; do ./ukbgene imp -c${i} -m -ak52049.key; done

#Download .bgen files
./ukbgene imp -c22 -ak52049.key
./ukbgene imp -c21 -ak52049.key
#And so on





