library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)
#https://www.10xgenomics.com/analysis-guides/background-removal-guidance-for-single-cell-gene-expression-datasets-using-third-party-tools
#package to read in cellbender to Seurat
#https://github.com/broadinstitute/CellBender/issues/315

data = Read_CellBender_h5_Mat("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCC/cellbender/HPAP-101-hg38.cellbender_FPR_0.01_filtered.h5")
obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
cbres <- NormalizeData(obj)
cbres <- FindVariableFeatures(cbres, selection.method = "vst", nfeatures = 2000)
cbres <- ScaleData(cbres)
cbres <- FindVariableFeatures(cbres)
cbres = RunPCA(cbres)
cbres <- FindNeighbors(cbres, dims = 1:10)
cbres <- FindClusters(algorithm = 1, cbres, resolution = 0.5)
cbres <- RunUMAP(cbres, dims = 1:9)
DimPlot(cbres, reduction = "umap")


#once you get which ones are DE between clusters need to use a list of known ones to then say what cell type it is
#from Wes 
markers <- list()
markers[["Islets"]] <- c('INS','GCG','MKI67','SST','PPY','GHRL')
markers[["Acinar"]] <- c('CTRB2','PRSS2','CPA1','REG1A') 
markers[["Ductal_muc5b"]] <- c('CFTR','MUC5B')
markers[["stellate"]] <- c('COL6A3','PDGFRB','IGA1','SPARCL1','SPARC','COL1A1')
markers[["Macrophage"]] <- c('C1QA','C1QB','C1QC')
markers[["T_cells"]] <- c('CD3D','CD4','CD8B','CXCR4')


# Create dotplot based on RNA expression
DotPlot(cbres, markers, assay="RNA") +
  RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))

#prev was from known markers, next do differentially expressed markers per cluster
markers <- FindAllMarkers(object = cbres, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

write.csv(markers, 
          file = "all_markers_annotated.csv", 
          quote = FALSE, 
          row.names = FALSE)

markers = read.csv("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCC/all_markers_annotated.csv")
top5 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, 
        wt = avg_log2FC)

#convert to correct list format that dotplot will take
library(tidyverse)
list2 = top5 %>%
  group_by(cluster) %>%
  summarise(named_vec = list(gene)) %>%
  deframe()

DotPlot(cbres, list2, assay="RNA") +
  RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))

#wes helped you ID a lot of your clusters - now you can properly name them
seurat_integrated <- RenameIdents(object = cbres, 
                                  "0" = "Beta",
                                  "1" = "Beta",
                                  "2" = "Alpha",
                                  "3" = "Alpha",
                                  "4" = "Weird",
                                  "5" = "Delta",
                                  "6" = "Beta",
                                  "7" = "Acinar",
                                  "8" = "Activated Stellate",
                                  "9" = "Beta",
                                  "10" = "Quiescent Stellate",
                                  "11" = "Beta with ambient contamination",
                                  "12" = "Ductal",
                                  "13" = "Gamma",
                                  "14" = "Potential Doublets",
                                  "15" = "Potential Mast Cells based on TPSB2 DE Marker",
                                  "16" = "Potential Doublets",
                                  "17" = "Potential Doublets",
                                  "18" = "Quiescent Stellate")

DimPlot(seurat_integrated, reduction = "umap", label=TRUE)


