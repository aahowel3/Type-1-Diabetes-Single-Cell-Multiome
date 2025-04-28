ereg_combined = read.csv("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCC/eRegulon_direct_Combined.tsv", sep ="\t", header = TRUE)
ereg_T1D_early = read.csv("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCC/eRegulon_direct_T1D_early.tsv", sep ="\t", header = TRUE)



length(unique(ereg_combined$TF))
length(unique(ereg_T1D_early$TF))
length(intersect(unique(ereg_combined$TF), unique(ereg_T1D_early$TF)))

intersected = intersect(unique(ereg_combined$TF), unique(ereg_T1D_early$TF))

bach1_genes_combined = ereg_combined[ereg_combined$TF == "FLI1",]
bach1_genes_combined = bach1_genes_combined$Gene


bach1_genes_t1d = ereg_T1D_early[ereg_T1D_early$TF == "FLI1",]
bach1_genes_t1d = bach1_genes_t1d$Gene

length(bach1_genes_combined)
length(bach1_genes_t1d)
length(intersect(bach1_genes_combined, bach1_genes_t1d))

for (x in intersected)  {
  
  bach1_genes_t1d = ereg_T1D_early[ereg_T1D_early$TF == x,]
  bach1_genes_t1d = bach1_genes_t1d$Gene
  
  bach1_genes_combined = ereg_combined[ereg_combined$TF == x,]
  bach1_genes_combined = bach1_genes_combined$Gene
  
  print(cat(x, length(intersect(bach1_genes_combined, bach1_genes_t1d)), length(bach1_genes_combined), length(bach1_genes_t1d) ))
  
  
}