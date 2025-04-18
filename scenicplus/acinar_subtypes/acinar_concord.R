pycis = read.csv("/tscc/projects/ps-gaultonlab/projects/npod1/cellranger_output/Multiome/Combined/cistopic_obj_celldata.csv")
og = as.data.frame(npod1_multome_shared@meta.data) 

updated_og = pycis %>% 
  inner_join(og, by="nPOD_barcode")

combine = paste(updated_og$CurrentObjLabels.x, updated_og$pycisTopic_leiden_10_0.6)
combine1 = paste(updated_og$CurrentObjLabels.x, updated_og$pycisTopic_leiden_10_1.2)
combine2 = paste(updated_og$CurrentObjLabels.x, updated_og$pycisTopic_leiden_10_3)

length(table(combine))
length(table(combine1))
length(table(combine2))

DimPlot(npod1_multome_shared, reduction = "umap", label=TRUE)

atacmeta$nPOD_barcode = atacmeta$X
atac_pycis = pycis %>% 
  inner_join(atacmeta, by="nPOD_barcode")

combine = paste(atac_pycis$predicted.id, atac_pycis$pycisTopic_leiden_10_3)
length(table(combine))










