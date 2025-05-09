library(Seurat)
library(irGSEA)
library(ComplexHeatmap)
library(dplyr)
library(Matrix)

#extract necessary atac_matrix.mtx, barcodes.tsv, and regions.tsv for the create_cistopic_object step
setwd("/tscc/projects/ps-gaultonlab/projects/npod1/cellranger_output/Multiome/Combined/")

#this is ALL of npod1 - includes the 8 multiome samples and others 
#the lenghts of the RNA and ATAC don't match because there are some that had more single modality than others 
npod1RNA <- readRDS("/tscc/projects/ps-gaultonlab/abhowell/npod1_data/RNA_final_20230515_rawMTxAssayIncl.rds")
npod1ATAC = readRDS("/tscc/projects/ps-gaultonlab/abhowell/npod1_data/20230515_snATAC_finalandFixedPeaksIncl_correctedFixPeak_wWindows.rds") 
dim(npod1RNA@meta.data)
dim(npod1ATAC@meta.data)

#this is the 8 multiome only object
#the rub is that there are no celltype annotations on this
npod1_multiome <- readRDS('/tscc/projects/ps-gaultonlab/rlmelton/megapanc/040725/240210_WE_nPOD_Multiome_hvw_window_adata_clustered_2:20_remove_samples.RDS')
npod1_multiome_meta=as.data.frame(npod1_multiome@meta.data)
dim(npod1_multiome)
#cleanest way to do this is subest npod1RNA and npod1ATAC by the correct 8 samples
#make sure all barcodes overlap 
npod1RNA_multome = subset(x = npod1RNA, subset = samples %in% c("MM_667","MM_666", "MM_662", "MM_665",
                                                                "MM_661", "MM_510", "MM_660", "MM_664"))

npod1ATAC_multome = subset(x = npod1ATAC, subset = library2 %in% c("MM_667","MM_666", "MM_662", "MM_665",
                                                                "MM_661", "MM_510", "MM_660", "MM_664"))


#still not same amount of barcodes will have to intersect them
dim(npod1RNA_multome@meta.data)
dim(npod1ATAC_multome@meta.data)
dim(npod1_multiome@meta.data)

#make sure all samples are included
table(npod1RNA_multome@meta.data$samples)
table(npod1ATAC_multome@meta.data$library2)

#rebecca agrees easiest is to combine ATAC and RNA - not to map celltype annotations onto already multiome
df1 = npod1RNA_multome@meta.data
df2 = npod1ATAC_multome@meta.data
df1check = df1[grep("AAACATGCAGTAATAG-1", rownames(df1)),]
df2check = df2[grep("AAACATGCAGTAATAG-1", rownames(df2)),]


library(dplyr)
#pull cell names
cells1 <- colnames(npod1RNA_multome)
cells2 <- colnames(npod1ATAC_multome)
# Identify shared cells between both objects
shared_cells <- intersect(x = cells1, y = cells2)
# Create Obj1 minus shared cells
npod1ATAC_multome_shared <- subset(npod1ATAC_multome, cells = shared_cells)
npod1RNA_multome_shared <- subset(npod1RNA_multome, cells = shared_cells)

atacmeta=(npod1ATAC_multome_shared@meta.data)
table(atacmeta$celltype)
table(atacmeta$predicted.id)

dim(npod1ATAC_multome_shared@meta.data)
dim(npod1RNA_multome_shared@meta.data)

#add atac onto RNA and check metadata dim lenght
npod1RNA_multome_shared[["ATAC"]] <- as(npod1ATAC_multome_shared[["ATAC"]], Class = "Assay5")
dim(npod1RNA_multome_shared@meta.data)
#rename object
npod1_multome_shared = npod1RNA_multome_shared


writeMM(npod1_multome_shared@assays$RNA$counts, "npod1_rna_counts.mtx")
writeMM(npod1_multome_shared@assays$ATAC$counts, "npod1_atac_counts.mtx")

write.table(
  npod1_multome_shared@meta.data, "npod1_meta_data.tsv", 
  col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(
  colnames(npod1_multome_shared), "npod1_barcodes.tsv", 
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(
  rownames(npod1_multome_shared[["RNA"]]), "npod1_gene_names.tsv", 
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(
  rownames(npod1_multome_shared[["ATAC"]]), "npod1_region_names.tsv", 
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")



