#setwd("/tscc/projects/ps-gaultonlab/abhowell/scenic_data2")
#You: phg001854 (ARRAY)  
#phg001854.v1.nPOD_T1D.sample-info.MULTI.tar.gz ======= Sample_phg001854 
#phg001854.v1.nPOD_T1D.genotype-calls-matrixfmt.c1.GRU.tar.gz ======== Plink
setwd("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/250128_nPOD_redownload/Plink") 

#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)

#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

#Abby - read bim into R, use dnSNP to change coordiantes to 38
#plink file is bim/fam/bed - so plink file is a bim file 
#bed file for plink is binary file of genotypes - not a normal bed 


bim = read.csv("nPOD_UFDIchip_renamed.bim")


# Read in the supporting data
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
snpcount(snps)
seqinfo(snps)

# Look up SNP positions
my_snps <- snpsById(snps, mvp.cred$CS.SNP[str_detect(mvp.cred$CS.SNP, 'rs')], ifnotfound='warning')
my_snps




