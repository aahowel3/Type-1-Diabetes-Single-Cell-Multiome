#You: phg001854 (ARRAY)  
#phg001854.v1.nPOD_T1D.sample-info.MULTI.tar.gz ======= Sample_phg001854 
#phg001854.v1.nPOD_T1D.genotype-calls-matrixfmt.c1.GRU.tar.gz ======== Plink
library(data.table)
library(dplyr)
setwd("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/250128_nPOD_redownload/Plink") 

#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
snps_37 = SNPlocs.Hsapiens.dbSNP155.GRCh37
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
snps_38 = SNPlocs.Hsapiens.dbSNP155.GRCh38

#Abby - read bim into R, use dnSNP to change coordiantes to 38
#plink file is bim/fam/bed - so plink file is a bim file 
#bed file for plink is binary file of genotypes - not a normal bed 

library(genio)
bim = read_bim("nPOD_UFDIchip_renamed.bim", verbose = TRUE)

##checking for lines in bim file that have more than 1 rsid seperate by a ;
ch <- ";"
count <- sapply(as.character(bim$id), 
                function(x, letter = ch){
                  str <- strsplit(x, split = "")
                  sum(unlist(str) == letter)
                })
print ("Count of !")
table(count)

bim$count = count
#there are affy arrays with NO rsid (count = zero) so filter that as well
bim2 = bim[bim$count == 1,]

#if the max amount of rsIDs is 6 max split is 7
bim_split = bim2 %>% separate(id, paste0('X', c(1:2)), sep = ';', remove = T)

snps_37 = as.data.frame(snpsById(snps_37, bim_split$X2, ifnotfound="drop"))
#intersect the complete bim list wiht those that appear in 37
#cant just use the 37 list because its missing the affy prefix
names(snps_37)[names(snps_37) == 'RefSNP_id'] <- 'X2'
bim_split_37 = bim_split %>% 
  inner_join(snps_37, by="X2")



snps_38_return = as.data.frame(snpsById(snps_38, bim_split_37$X2, genome="GRCh38", ifnotfound="drop"))
#some lost in 38 - intersect 
names(snps_38_return)[names(snps_38_return) == 'RefSNP_id'] <- 'X2'


bim_final_38 = bim_split_37 %>% 
  inner_join(snps_38_return, by="X2")

bim_final_38$id = paste0(bim_final_38$X1,";",bim_final_38$X2)


#pull final columns to mimic bim file
bim_final_38 = bim_final_38[c("chr", "id","posg","pos","alt","ref")]
bim_final_38 = unique(bim_final_38)

write.table(bim_final_38$id, "selected_rsIDS.csv", row.names = FALSE, quote=FALSE)


write.table(bim_final_38, "nPOD_UFDIchip_renamed_filtered.bim", row.names = FALSE, quote=FALSE)
