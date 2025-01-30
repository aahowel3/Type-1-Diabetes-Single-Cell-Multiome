library(Seurat)
library(irGSEA)
library(ComplexHeatmap)
library(dplyr)
setwd("/tscc/projects/ps-gaultonlab/abhowell/scenic_data/")

#Kegg pathways
npod2 <- readRDS("031424_PancGRS_8lanes_merged_filt5perMito_100katacRNAcounts_peakCallLieden_celltypesLabelled_wlymph_rmLowQualAcinar_peakCallCelltype_acinarSubtypes_usingNPOD1refmap_newWindowsAssayNoCHR_finalPeakCalls_wGeneAct.rds")

npod2_meta=as.data.frame(npod2@meta.data)


#npod2_meta$newrow=paste(rownames(npod2_meta), npod2_meta$orig.ident, sep = "-")
#use the actual donor IDs from demuxlet to label your metadata
npod2_meta$newrow=paste(rownames(npod2_meta), npod2_meta$donor_demux, sep = "-")


rownames(npod2_meta) = npod2_meta$newrow


write.table(npod2_meta, file = "npod2_metadata_donordemuxids.tsv", row.names=TRUE, sep="\t")

