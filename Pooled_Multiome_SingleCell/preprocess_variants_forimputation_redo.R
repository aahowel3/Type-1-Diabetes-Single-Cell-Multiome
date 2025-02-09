#setwd("/tscc/projects/ps-gaultonlab/abhowell/scenic_data2")
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
#awesome
#up to SIX different rsIDs per variant fuck
#if the max amount of rsIDs is 6 max split is 7
bim_split2 = bim %>% separate(id, paste0('X', c(1:7)), sep = ';', remove = T)
#well we can safely drop the rows where there is just no rsID
bim_split2[is.na(bim_split2)] = "rs0"


#do this 6 times, glue those dataframes together rowwise
snps_37col1 = as.data.frame(snpsById(snps_37, bim_split2$X2, ifnotfound="drop"))
names(snps_37col1)[names(snps_37col1) == 'RefSNP_id'] <- 'X2'
bim_split2_col1 = bim_split2 %>% 
  inner_join(snps_37col1, by="X2")


#do this 6 times, glue those dataframes together rowwise
snps_37col2 = as.data.frame(snpsById(snps_37, bim_split2$X3, ifnotfound="drop"))
names(snps_37col2)[names(snps_37col2) == 'RefSNP_id'] <- 'X3'
bim_split2_col2 = bim_split2 %>% 
  inner_join(snps_37col2, by="X3")


#do this 6 times, glue those dataframes together rowwise
snps_37col3 = as.data.frame(snpsById(snps_37, bim_split2$X4, ifnotfound="drop"))
names(snps_37col3)[names(snps_37col3) == 'RefSNP_id'] <- 'X4'
bim_split2_col3 = bim_split2 %>% 
  inner_join(snps_37col3, by="X4")


#do this 6 times, glue those dataframes together rowwise
snps_37col4 = as.data.frame(snpsById(snps_37, bim_split2$X5, ifnotfound="drop"))
names(snps_37col4)[names(snps_37col4) == 'RefSNP_id'] <- 'X5'
bim_split2_col4 = bim_split2 %>% 
  inner_join(snps_37col4, by="X5")

#do this 6 times, glue those dataframes together rowwise
snps_37col5 = as.data.frame(snpsById(snps_37, bim_split2$X6, ifnotfound="drop"))
names(snps_37col5)[names(snps_37col5) == 'RefSNP_id'] <- 'X6'
bim_split2_col5 = bim_split2 %>% 
  inner_join(snps_37col5, by="X6")


#do this 6 times, glue those dataframes together rowwise
snps_37col6 = as.data.frame(snpsById(snps_37, bim_split2$X7, ifnotfound="drop"))
names(snps_37col6)[names(snps_37col6) == 'RefSNP_id'] <- 'X7'
bim_split2_col6 = bim_split2 %>% 
  inner_join(snps_37col6, by="X7")

#drop duplicates removes the double return where the first AND second rsID returned a "keep" 
bim_filter = rbind(bim_split2_col1,bim_split2_col2,bim_split2_col3,bim_split2_col4,bim_split2_col5,bim_split2_col6)
bim_filter_nodup = bim_filter[!duplicated(bim_filter), ]

#some variants that appear in 37 do not appear in 38 - can't filter bed file from the 37 variants that exist and then try to use 
#it with the 38 - lenghts will be differetn
#this one should not change anything - no instances where the first rsID was a hit but the second wasnt
bim_filter_nodup = bim_filter_nodup[!bim_filter_nodup$X2 == "rs0", ]
#here we want to keep the missing entries in rsiD col 2, don't want to have to bother with the double mapping
bim_filter_nodup = bim_filter_nodup[bim_filter_nodup$X3 == "rs0", ]


snps_38_return = as.data.frame(snpsById(snps_38, bim_filter_nodup$X2, genome="GRCh38", ifnotfound="drop"))
#some lost in 38 - intersect 
names(snps_38_return)[names(snps_38_return) == 'RefSNP_id'] <- 'X2'
bim_split2 = bim_filter_nodup %>% 
  inner_join(snps_38_return, by="X2")

check=bim_split2[bim_split2$pos.x != bim_split2$pos.y,]


bim_split2$id = paste0(bim_split2$X1,";",bim_split2$X2)
write.table(bim_filter_nodup$id, "selected_rsIDS.csv", row.names = FALSE, quote=FALSE)

#conda activate
#plink --bfile nPOD_UFDIchip_renamed --extract  selected_rsIDS.csv --keep-allele-order --make-bed --out nPOD_UFDIchip_renamed_filtered



seqlevelsStyle(check38) <- "UCSC"
inferRefAndAltAlleles(check38, BSgenome.Hsapiens.UCSC.hg38)


#plink --bfile nPOD_UFDIchip_renamed --extract  selected_rsIDS.csv --keep-allele-order --make-bed --out nPOD_UFDIchip_renamed_filtered

#bim_final_38 = bim_split2[c("seqnames", "id","posg","pos.y","alt_alleles","ref_allele")]
#instead of taking new ref/alt alleles take the old ones because they don't have multiple alleles
bim_final_38 = bim_split2[c("seqnames", "id","posg","pos.y","alt","ref")]

names(bim_final_38) = c("chr", "id","posg","pos","alt","ref")





#/tscc/nfs/home/abhowell/imputeInversion/HRC-1000G-check-bim.pl -b nPOD_UFDIchip_renamed_filtered_step4.bim -f nPOD_UFDIchip_renamed_filtered_step4.frq -r PASS.Variantsbravo-dbsnp-all.tab.gz -h







#file to prepare for bim/bed file modify command
#plink1.9 --bfile PREFIX_OF_PLINK_FILES --extract VARIANT_LIST --keep-allele-order --make-bed --out NEW_PLINK_PREFIX
write.table(bim_filter_nodup$id, "selected_rsIDS_redo.csv", row.names = FALSE, quote=FALSE)

#conda activate plink
#plink --bfile nPOD_UFDIchip_renamed --extract  selected_rsIDS_redo.csv --keep-allele-order --make-bed --out nPOD_UFDIchip_renamed_filtered
bim2 = read_bim("nPOD_UFDIchip_renamed_filtered_old.bim", verbose = TRUE)
##checking for lines in bim file that have more than 1 rsid seperate by a ;
ch <- ";"
count <- sapply(as.character(bim2$id), 
                function(x, letter = ch){
                  str <- strsplit(x, split = "")
                  sum(unlist(str) == letter)
                })
print ("Count of !")
table(count)
bim2$count = count
#68/900,000 drop the double entry ones
bim2 = bim2[bim2$count == 1, ]
#drop the ones with two entries - 

bim_split = bim2 %>% separate(id, paste0('X', c(1:2)), sep = ';', remove = T)
snps_38_return = as.data.frame(snpsById(snps_38, bim_split$X2, genome="GRCh38", ifnotfound="drop"))
#some lost in 38 - intersect 
names(snps_38_return)[names(snps_38_return) == 'RefSNP_id'] <- 'X2'
bim_split2 = bim_split %>% 
  inner_join(snps_38_return, by="X2")

bim_split2$id = paste0(bim_split2$X1,";",bim_split2$X2)
#create new GH38 bim using snps_38_return info 

check38 = snpsById(snps_38, bim_split$X2, ifnotfound="drop")
seqlevelsStyle(check38) <- "UCSC"
inferRefAndAltAlleles(check38, BSgenome.Hsapiens.UCSC.hg38)


#bim_final_38 = bim_split2[c("seqnames", "id","posg","pos.y","alt_alleles","ref_allele")]
#instead of taking new ref/alt alleles take the old ones because they don't have multiple alleles
bim_final_38 = bim_split2[c("seqnames", "id","posg","pos.y","alt","ref")]

names(bim_final_38) = c("chr", "id","posg","pos","alt","ref")


#head_bim = head(bim_final_38,10)


bim_final_38 = bim_final_38 %>% rowwise() %>% 
  mutate(alt = paste(alt, collapse=',')) %>%
  ungroup()

#bim_final_38_nomulti = bim_final_38[- grep(",", bim_final_38$alt),]


#in tscc mv nPOD_UFDIchip_renamed_filtered.bim nPOD_UFDIchip_renamed_filtered_old.bim
write.table(bim_final_38, "nPOD_UFDIchip_renamed_filtered.bim", row.names = FALSE, quote=FALSE, sep="\t")


write.table(bim_final_38$id, "selected_rsIDS_secondfilter.csv", row.names = FALSE, quote=FALSE)

#plink --bfile nPOD_UFDIchip_renamed --extract  selected_rsIDS_secondfilter.csv --keep-allele-order --make-bed --out nPOD_UFDIchip_renamed_filtered



#checking something else
library(data.table)
library(dplyr)
setwd("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/250128_nPOD_redownload/Plink") 


library(genio)
bim_myedits = read_bim("nPOD_UFDIchip_renamed_filtered_step4-updated-chr6.bim", verbose = TRUE)
bim_noedits = read_bim("nPOD_UFDIchip_renamed-updated-chr6.bim", verbose = TRUE)


bim_split= bim_noedits %>% 
  inner_join(bim_myedits, by="id")


bim_split[bim_split$alt.x != bim_split$alt.y,]
bim_split[bim_split$ref.x != bim_split$ref.y,]
bim_split[bim_split$pos.x != bim_split$pos.y,]

#back up, did any of my changes made in R stick? Or did it just drop those?

bim_myedits = read_bim("nPOD_UFDIchip_renamed_filtered_step4.bim", verbose = TRUE)
bim_noedits = read_bim("nPOD_UFDIchip_renamed.bim", verbose = TRUE)

bim_split= bim_noedits %>% 
  inner_join(bim_myedits, by="id")


bim_split[bim_split$pos.x != bim_split$pos.y,]



