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


library(genio)
bim = read_bim("nPOD_UFDIchip_renamed.bim", verbose = TRUE)
bim_split <- data.frame(do.call('rbind', strsplit(as.character(bim$id),';',fixed=TRUE)))
#col split is odd bc some rows have 1-3 rsids
varlist = as.data.frame(unlist(c(bim_split$X1, bim_split$X2, bim_split$X3,bim_split$X4,bim_split$X5,bim_split$X6,bim_split$X7)))
varlist_un=unique(varlist)
varlist_final = dplyr::filter(varlist_un, !grepl("Aff", `unlist(c(bim_split$X1, bim_split$X2, bim_split$X3, bim_split$X4, bim_split$X5, bim_split$X6, bim_split$X7))`))

return_37 = snpsById(snps_37, varlist_final$`unlist(c(bim_split$X1, bim_split$X2, bim_split$X3, bim_split$X4, bim_split$X5, bim_split$X6, bim_split$X7))`, ifnotfound="drop")
return_37= as.data.frame(return_37)

return_38 = snpsById(snps_38, varlist_final$`unlist(c(bim_split$X1, bim_split$X2, bim_split$X3, bim_split$X4, bim_split$X5, bim_split$X6, bim_split$X7))`, ifnotfound="drop")
return_38= as.data.frame(return_38)


##checking for lines in bim file that have more than 1 rsid seperate by a ;
ch <- ";"
count <- sapply(as.character(bim$id), 
                function(x, letter = ch){
                  str <- strsplit(x, split = "")
                  sum(unlist(str) == letter)
                })
print ("Count of !")
# returning the number of occurrences
bim$count = count
bim_issue = bim[bim$count > 1, ]
#checking for alt/ref columns that are INDELs based on character count
bim_issue$alt_coutn = nchar(bim_issue$alt)
bim_issue$ref_coutn = nchar(bim_issue$ref)
bim_issue$combined = bim_issue$alt_coutn + bim_issue$ref_coutn
bim_sevissue = bim_issue[bim_issue$combined == 2, ]

#file to prepare for bim/bed file modify command
#plink1.9 --bfile PREFIX_OF_PLINK_FILES --extract VARIANT_LIST --keep-allele-order --make-bed --out NEW_PLINK_PREFIX
write.table(return_37$RefSNP_id, "selected_rsIDS.csv", row.names = FALSE, quote=FALSE)


