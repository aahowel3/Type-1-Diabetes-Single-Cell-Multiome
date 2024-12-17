library(Seurat)
library(SeuratData)
library(SeuratDisk)

setwd("/tscc/nfs/home/abhowell/hubmappipeline/salmon-rnaseq/")
Convert("secondary_analysis.h5ad", dest = "h5seurat", overwrite = TRUE)
hubmap <- LoadH5Seurat("secondary_analysis.h5seurat")

#here this objects features are called n_genes, n_counts, leiden, and umap density 
#VlnPlot(hubmap, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#check by doing hubmap$
#because there is already a leiden and umap column it has already gone through normalization/PCA etc
#remember your cellbender object started out as "obj" with 3 features (original barcode/gene count/gene name files)
#and then became 5 features as "cbres" because you did the clustering 
#to figure out what went here looked a cbres object - umap is not a $ feature it is a @reductions$
hubmap@reductions$umap
DimPlot(hubmap, reduction = "umap")
#hmm so why isn't it colored in nicely? 
#change ident to leiden clusters,  
hubmap@active.ident = hubmap$leiden 

#once you get which ones are DE between clusters need to use a list of known ones to then say what cell type it is
#from Wes 
markers <- list()
markers[["Islets"]] <- c('INS','GCG','MKI67','SST','PPY','GHRL')
markers[["Acinar"]] <- c('CTRB2','PRSS2','CPA1','REG1A') 
markers[["Ductal_muc5b"]] <- c('CFTR','MUC5B')
markers[["stellate"]] <- c('COL6A3','PDGFRB','IGA1','SPARCL1','SPARC','COL1A1')
markers[["Macrophage"]] <- c('C1QA','C1QB','C1QC')
markers[["T_cells"]] <- c('CD3D','CD4','CD8B','CXCR4')


#cannot convert gene names after seurat object made absoluute pain
#convert markers to ENSG and then reconvert from there
#fix in hubmap pipeline later
library("AnnotationDbi")
library("org.Hs.eg.db")
markers2 = lapply(markers, function(x) mapIds(org.Hs.eg.db, keys=x, column="ENSEMBL", keytype="SYMBOL", multiVals="first"))

#STILL NOT RIGHT 
#salmon used ENSGs with a version number
#screw this grab list of originail names and glue extensions onto translated ensembls
library(tidyverse) 
ver_ensg=read.csv("salmon_out/alevin/quants_mat_cols.txt", header=FALSE)
ver_ensg2 = ver_ensg %>% separate_wider_delim(V1, ".", names = c("A", "B"))
ver_ensg2$C = paste(ver_ensg2$A,ver_ensg2$B, sep=".")


markers2df$ensg = ver_ensg2$C[match(markers2df$value, ver_ensg2$A)]
list2 = markers2df %>%
  group_by(name) %>%
  summarise(named_vec = list(ensg)) %>%
  deframe()

#how to check remaining genes 
gene_list <- rownames(hubmap)
gene_list = as.data.frame(gene_list)

# Create dotplot based on RNA expression
DotPlot(hubmap, list2, assay="RNA") +
  RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7)) +
  scale_x_discrete(labels=markers)

#wes helped you ID a lot of your clusters - now you can properly name them
seurat_integrated <- RenameIdents(object = hubmap, 
                                  "0" = "Beta",
                                  "1" = "Alpha",
                                  "2" = "Beta",
                                  "3" = "Beta",
                                  "4" = "Activated Stellate",
                                  "5" = "T Cells (?)",
                                  "6" = "Acinar",
                                  "7" = "Quiescent Stellate",
                                  "8" = "Potential Doublets")

DimPlot(seurat_integrated, reduction = "umap", label=TRUE)

