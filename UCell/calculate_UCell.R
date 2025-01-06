library(Seurat)
library(irGSEA)
library(ComplexHeatmap)
library(clusterProfiler)
library(tidyverse)
library(DBI)
library(RSQLite)
setwd("/tscc/projects/ps-gaultonlab/abhowell/npod1_data/")

kk <- clusterProfiler::gson_KEGG(species = "hsa")
gson::write.gson(kk, file = "./KEGG_20231128.gson")

# read gson file
kk2 <- gson::read.gson("./KEGG_20231128.gson")
# Convert to a data frame
kegg.list <- dplyr::left_join(kk2@gsid2name,
                              kk2@gsid2gene,
                              by = "gsid")
head(kegg.list)
# gsid               name      gene
# 1 hsa01100 Metabolic pathways        10
# 2 hsa01100 Metabolic pathways       100
# 3 hsa01100 Metabolic pathways     10005
# 4 hsa01100 Metabolic pathways     10007
# 5 hsa01100 Metabolic pathways 100137049
# 6 hsa01100 Metabolic pathways     10020

# Convert gene ID to gene symbol
gene_name <- clusterProfiler::bitr(kegg.list$gene, 
                                   fromType = "ENTREZID", 
                                   toType = "SYMBOL", 
                                   OrgDb = "org.Hs.eg.db")
kegg.list <- dplyr::full_join(kegg.list,
                              gene_name,
                              by = c("gene"="ENTREZID"))
# remove NA value if exist
kegg.list <- kegg.list[complete.cases(kegg.list[, c("gene", "SYMBOL")]), ]
head(kegg.list)
# gsid               name      gene  SYMBOL
# 1 hsa01100 Metabolic pathways        10    NAT2
# 2 hsa01100 Metabolic pathways       100     ADA
# 3 hsa01100 Metabolic pathways     10005   ACOT8
# 4 hsa01100 Metabolic pathways     10007  GNPDA1
# 5 hsa01100 Metabolic pathways 100137049 PLA2G4B
# 6 hsa01100 Metabolic pathways     10020     GNE

# convert to list required by irGSEA package
kegg.list$name <- factor(kegg.list$name)
kegg.list <- kegg.list %>% 
  dplyr::group_split(name, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(SYMBOL) %>% unique(.)) %>%
  purrr::set_names(levels(kegg.list$name))
head(kegg.list)

### go bp ###
# download go bp (human) and write as gson file
go <- clusterProfiler::gson_GO(OrgDb = "org.Hs.eg.db", ont = "BP")
gson::write.gson(go, file = "./go_20231128.gson")

# read gson file
go2 <- gson::read.gson("./go_20231128.gson")

# Convert to a data frame
go.list <- dplyr::left_join(go2@gsid2name,
                            go2@gsid2gene,
                            by = "gsid")
head(go.list)
# gsid                             name gene
# 1 GO:0000001        mitochondrion inheritance <NA>
#   2 GO:0000002 mitochondrial genome maintenance  142
# 3 GO:0000002 mitochondrial genome maintenance  291
# 4 GO:0000002 mitochondrial genome maintenance 1763
# 5 GO:0000002 mitochondrial genome maintenance 1890
# 6 GO:0000002 mitochondrial genome maintenance 2021

# Convert gene ID to gene symbol
go.list <- dplyr::full_join(go.list,
                            go2@gene2name,
                            by = c("gene"="ENTREZID"))
# remove NA value if exist
go.list <- go.list[complete.cases(go.list[, c("gene", "SYMBOL")]), ]
head(go.list)
# gsid                             name gene  SYMBOL
# 2 GO:0000002 mitochondrial genome maintenance  142   PARP1
# 3 GO:0000002 mitochondrial genome maintenance  291 SLC25A4
# 4 GO:0000002 mitochondrial genome maintenance 1763    DNA2
# 5 GO:0000002 mitochondrial genome maintenance 1890    TYMP
# 6 GO:0000002 mitochondrial genome maintenance 2021   ENDOG
# 7 GO:0000002 mitochondrial genome maintenance 3980    LIG3

# convert to list required by irGSEA package
go.list$name <- factor(go.list$name)
go.list <- go.list %>% 
  dplyr::group_split(name, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(SYMBOL) %>% unique(.)) %>%
  purrr::set_names(levels(go.list$name))
head(go.list)


# code modified by https://rdrr.io/github/cashoes/sear/src/data-raw/1_parse_msigdb_sqlite.r
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = './msigdb_v2024.1.Hs.db')
DBI::dbListTables(con)


geneset_db <- dplyr::tbl(con, 'gene_set')                                              # standard_name, collection_name
details_db <- dplyr::tbl(con, 'gene_set_details')                                      # description_brief, description_full
geneset_genesymbol_db <- dplyr::tbl(con, 'gene_set_gene_symbol')                       # meat and potatoes
genesymbol_db <- dplyr::tbl(con, 'gene_symbol')                                        # mapping from ids to gene symbols
collection_db <- dplyr::tbl(con, 'collection') %>% dplyr::select(collection_name, full_name)  # collection metadata

# join tables
msigdb <- geneset_db %>%
  dplyr::left_join(details_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(collection_db, by = 'collection_name') %>%
  dplyr::left_join(geneset_genesymbol_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(genesymbol_db, by = c('gene_symbol_id' = 'id')) %>%
  dplyr::select(collection = collection_name, subcollection = full_name, geneset = standard_name, description = description_brief, symbol) %>%
  dplyr::as_tibble() 

# clean up
DBI::dbDisconnect(con)


unique(msigdb$collection)
# [1] "C1"                 "C2:CGP"             "C2:CP:BIOCARTA"    
# [4] "C2:CP:KEGG_LEGACY"  "C2:CP:PID"          "C3:MIR:MIRDB"      
# [7] "C3:MIR:MIR_LEGACY"  "C3:TFT:GTRD"        "C3:TFT:TFT_LEGACY" 
# [10] "C4:3CA"             "C4:CGN"             "C4:CM"             
# [13] "C6"                 "C7:IMMUNESIGDB"     "C7:VAX"            
# [16] "C8"                 "C5:GO:BP"           "C5:GO:CC"          
# [19] "C5:GO:MF"           "H"                  "C5:HPO"            
# [22] "C2:CP:KEGG_MEDICUS" "C2:CP:REACTOME"     "C2:CP:WIKIPATHWAYS"
# [25] "C2:CP" 
unique(msigdb$subcollection)
# [1] "C1"                 "C2:CGP"             "C2:CP:BIOCARTA"    
# [4] "C2:CP:KEGG_LEGACY"  "C2:CP:PID"          "C3:MIR:MIRDB"      
# [7] "C3:MIR:MIR_LEGACY"  "C3:TFT:GTRD"        "C3:TFT:TFT_LEGACY" 
# [10] "C4:3CA"             "C4:CGN"             "C4:CM"             
# [13] "C6"                 "C7:IMMUNESIGDB"     "C7:VAX"            
# [16] "C8"                 "C5:GO:BP"           "C5:GO:CC"          
# [19] "C5:GO:MF"           "H"                  "C5:HPO"            
# [22] "C2:CP:KEGG_MEDICUS" "C2:CP:REACTOME"     "C2:CP:WIKIPATHWAYS"
# [25] "C2:CP" 

# convert to list[hallmarker] required by irGSEA package
msigdb.h <- msigdb %>% 
  dplyr::filter(collection=="H") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.h$geneset <- factor(msigdb.h$geneset)
msigdb.h <- msigdb.h %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.h$geneset))

# convert to list[go bp] required by irGSEA package
msigdb.go.bp <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:BP") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.bp$geneset <- factor(msigdb.go.bp$geneset)
msigdb.go.bp <- msigdb.go.bp %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.bp$geneset))

# convert to list[KEGG] required by irGSEA package
msigdb.kegg <- msigdb %>% 
  dplyr::filter(collection=="C2:CP:KEGG_MEDICUS") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.kegg$geneset <- factor(msigdb.kegg$geneset)
msigdb.kegg <- msigdb.kegg %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.kegg$geneset))


msigdb.react <- msigdb %>% 
  dplyr::filter(collection=="C2:CP:REACTOME") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.react$geneset <- factor(msigdb.react$geneset)
msigdb.react <- msigdb.react %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.react$geneset))

# Look for the gene sets associated with angiogenesis from gene sets names and 
# gene sets descriptions

category <- c("angiogenesis", "vessel")

msigdb.vessel <- list()
for (i in category) {
  # Ignore case matching
  find.index.description <- stringr::str_detect(msigdb$description, pattern = regex(all_of(i), ignore_case=TRUE))
  find.index.name <- stringr::str_detect(msigdb$geneset, pattern = regex(all_of(i), ignore_case=TRUE))
  msigdb.vessel[[i]] <- msigdb[find.index.description | find.index.name, ] %>% mutate(category = i)
  
}
msigdb.vessel <- do.call(rbind, msigdb.vessel)

head(msigdb.vessel)
# # A tibble: 6 × 6
# collection subcollection                      geneset            description    symbol category
# <chr>      <chr>                              <chr>              <chr>          <chr>  <chr>   
#   1 C2:CGP     Chemical and Genetic Perturbations HU_ANGIOGENESIS_UP Up-regulated … HECW1  angioge…
# 2 C2:CGP     Chemical and Genetic Perturbations HU_ANGIOGENESIS_UP Up-regulated … JADE2  angioge…
# 3 C2:CGP     Chemical and Genetic Perturbations HU_ANGIOGENESIS_UP Up-regulated … SEMA3C angioge…
# 4 C2:CGP     Chemical and Genetic Perturbations HU_ANGIOGENESIS_UP Up-regulated … STUB1  angioge…
# 5 C2:CGP     Chemical and Genetic Perturbations HU_ANGIOGENESIS_UP Up-regulated … FAH    angioge…
# 6 C2:CGP     Chemical and Genetic Perturbations HU_ANGIOGENESIS_UP Up-regulated … COL7A1 angioge…

length(unique(msigdb.vessel$geneset))
# [1] 112

# convert gene sets associate



msigdb <- clusterProfiler::read.gmt("./msigdb.v2023.2.Hs.symbols.gmt")

# convert to list[hallmarker] required by irGSEA package
msigdb.h <- msigdb %>% 
  dplyr::filter(str_detect(term, pattern = regex("HALLMARK_", ignore_case=TRUE)))
msigdb.h$term <- factor(msigdb.h$term)
msigdb.h <- msigdb.h %>% 
  dplyr::group_split(term, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(gene) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.h$term))

# convert to list[go bp] required by irGSEA package

msigdb.go.bp <- msigdb %>% 
  dplyr::filter(str_detect(term, pattern = regex("GOBP_", ignore_case=TRUE)))
msigdb.go.bp$term <- factor(msigdb.go.bp$term)
msigdb.go.bp <- msigdb.go.bp %>% 
  dplyr::group_split(term, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(gene) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.bp$term))

# convert to list[KEGG] required by irGSEA package
msigdb.kegg <- msigdb %>% 
  dplyr::filter(str_detect(term, pattern = regex("KEGG_", ignore_case=TRUE)))
msigdb.kegg$term <- factor(msigdb.kegg$term)
msigdb.kegg <- msigdb.kegg %>% 
  dplyr::group_split(term, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(gene) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.kegg$term))









######################################
#######################################
########################################
#after setting up the irgsea databases load in the NPOD1 data and the HPAP data to Ucell scores for each 
npod1 <- readRDS("RNA_final_20230515_rawMTxAssayIncl.rds")
#need to assing subtypes via sep metadata file and then subset
npod1_meta=read.csv("nPOD_snRNA_cellCount_subtypes_wProportions_wMeanGeneRAW_wMT_20230515.txt", sep='\t')
npod1$nPOD_ID_2 = gsub("multi_", "", npod1$nPOD_ID)
npod1$condition_subtype = npod1_meta$Condition_subtype[match(npod1$nPOD_ID_2, npod1_meta$Sample.ID)]

npod1_kegg <- irGSEA.score(object = npod1, assay = "RNA", 
                              slot = "data", seeds = 123, ncores = 4,
                              min.cells = 3, min.feature = 0,
                              maxGSSize = 500,
                              custom = T, geneset = msigdb.kegg, 
                              method = c("UCell"),
                              aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                              kcdf = 'Gaussian')

SaveSeuratRds(npod1_kegg, file = "npod1_kegg_UCell.rds")

npod1_react <- irGSEA.score(object = npod1, assay = "RNA", 
                           slot = "data", seeds = 123, ncores = 4,
                           min.cells = 3, min.feature = 0,
                           maxGSSize = 500,
                           custom = T, geneset = msigdb.react, 
                           method = c("UCell"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian')

SaveSeuratRds(npod1_react, file = "npod1_react_UCell.rds")



#read in HPAP
hpap <- readRDS("2023_01_27_hpap_res0p3.rds")
#use HPAP metadatafile now too
hpap_meta= read.csv("hpap_metadata2_rebecca.tsv", sep="\t")
#apply labels from metadata to the seurat object
hpap$condition_subtype = hpap_meta$diabetes_status[match(hpap$library, hpap_meta$donor_ID)]


hpap_kegg <- irGSEA.score(object = hpap, assay = "RNA", 
                                     slot = "data", seeds = 123, ncores = 4,
                                     min.cells = 3, min.feature = 0,
                                     maxGSSize = 500,
                                     custom = T, geneset = msigdb.kegg, 
                                     method = c("UCell"),
                                     aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                                     kcdf = 'Gaussian')

SaveSeuratRds(hpap_kegg, file = "hpap_kegg_UCell.rds")


hpap_react <- irGSEA.score(object = hpap, assay = "RNA", 
                          slot = "data", seeds = 123, ncores = 4,
                          min.cells = 3, min.feature = 0,
                          maxGSSize = 500,
                          custom = T, geneset = msigdb.react, 
                          method = c("AUCell", "UCell"),
                          aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

SaveSeuratRds(hpap_react, file = "hpap_react_UCell.rds")

#in "/tscc/projects/ps-gaultonlab/abhowell/npod1_data"
npod1_kegg <- readRDS("npod1_kegg_UCell.rds")
npod1_kegg_Ucell=as.data.frame(npod1_kegg@assays$UCell$scale.data)
npod1_kegg_Ucell <- data.frame(t(npod1_kegg_Ucell))
npod1_kegg_Ucell$cell = rownames(npod1_kegg_Ucell)
#need to assign each cell barcode their cluster
npod1_kegg_Ucell$celltype = as.data.frame(npod1_kegg$celltype_assignment2)
npod1_kegg_Ucell$condition = npod1_kegg@meta.data$condition_subtype
npod1_kegg_Ucell$donor = npod1_kegg@meta.data$samples

#read in npod react and combine
npod1_react <- readRDS("npod1_react_UCell.rds")
npod1_react_Ucell=as.data.frame(npod1_react@assays$UCell$scale.data)
npod1_react_Ucell <- data.frame(t(npod1_react_Ucell))
npod1_react_Ucell$cell = rownames(npod1_react_Ucell)
#need to assign each cell barcode their cluster
npod1_react_Ucell$celltype = as.data.frame(npod1_react$celltype_assignment2)
npod1_react_Ucell$condition = npod1_react@meta.data$condition_subtype
npod1_react_Ucell$donor = npod1_react@meta.data$samples

npod1_paths=cbind(npod1_kegg_Ucell, npod1_react_Ucell)

#that's right for these we just dropped individuals we knew to be single AAB in NPOD
#we did not have specific single AAB or multi AAB labels
npod1_paths=npod1_paths[!npod1_paths$donor == "MM_403",]
npod1_paths=npod1_paths[!npod1_paths$donor == "MM_401",]

#select only beta cell
npod1_paths=npod1_paths[npod1_paths$celltype$`npod1_kegg$celltype_assignment2` == 'Beta',]

#fix odd column name
npod1_paths$celltype = npod1_paths$celltype$`npod1_kegg$celltype_assignment2`

#use only pathways that were signifigant
kegglist = read.csv("kegglist.csv", header = FALSE)
kegglist$V1 = gsub("_", ".", kegglist$V1)

#second option for reactome paths
reactlist = read.csv("reactome_list_drop2.csv", header = FALSE)
reactlist$V1 = gsub("_", ".", reactlist$V1)
reactlist_final = reactlist$V1[reactlist$V1 %in% colnames(npod1_paths)]

npod1_paths=npod1_paths[,c(kegglist$V1,reactlist_final,"cell","donor","celltype","condition")]

npod1_paths$database= "npod"

#in "/tscc/projects/ps-gaultonlab/abhowell/npod1_data"
hpap_kegg <- readRDS("hpap_kegg_UCell.rds")
hpap_kegg_Ucell=as.data.frame(hpap_kegg@assays$UCell$scale.data)
hpap_kegg_Ucell <- data.frame(t(hpap_kegg_Ucell))
hpap_kegg_Ucell$cell = rownames(hpap_kegg_Ucell)
#need to assign each cell barcode their cluster
hpap_kegg_Ucell$celltype = as.data.frame(hpap_kegg$cell_type)
hpap_kegg_Ucell$condition = hpap_kegg@meta.data$condition_subtype
hpap_kegg_Ucell$donor = hpap_kegg@meta.data$library



#read in npod react and combine
hpap_react <- readRDS("hpap_react_UCell.rds")
hpap_react_Ucell=as.data.frame(hpap_react@assays$UCell$scale.data)
hpap_react_Ucell <- data.frame(t(hpap_react_Ucell))
hpap_react_Ucell$cell = rownames(hpap_react_Ucell)
#need to assign each cell barcode their cluster
hpap_react_Ucell$celltype = as.data.frame(hpap_react$cell_type)
hpap_react_Ucell$condition = hpap_react@meta.data$condition_subtype
hpap_react_Ucell$donor = hpap_react@meta.data$library

hpap_paths=cbind(hpap_kegg_Ucell, hpap_react_Ucell)

#drop few multi AAbs in HPAP
#NAs were T2Ds drop those too
hpap_paths=hpap_paths[!hpap_paths$condition == "Multiple",]
hpap_paths=hpap_paths[!is.na(hpap_paths$condition),]


#select only beta cell
hpap_paths=hpap_paths[hpap_paths$celltype$`hpap_kegg$cell_type` == 'Beta',]

#fix odd column name
hpap_paths$celltype = hpap_paths$celltype$`hpap_kegg$cell_type`

#use only pathways that were signifigant
kegglist = read.csv("kegglist.csv", header = FALSE)
kegglist$V1 = gsub("_", ".", kegglist$V1)

#second option for reactome paths
reactlist = read.csv("reactome_list_drop2.csv", header = FALSE)
reactlist$V1 = gsub("_", ".", reactlist$V1)
reactlist_final = reactlist$V1[reactlist$V1 %in% colnames(hpap_paths)]

hpap_paths=hpap_paths[,c(kegglist$V1,reactlist_final,"cell","donor","celltype","condition")]

hpap_paths$database= "hpap"

combined_data=rbind(npod1_paths,hpap_paths)

#encodings for disease state are not the same - fix that
unique(combined_data$condition)
library(dplyr)
combined_data = combined_data %>%
  mutate(reclassify = case_when(condition == "ND" ~ 'ND',
                                condition == "T1D_early" ~ 'Early',
                                condition == "Early" ~ 'Early',
                                condition == "T1D_late" ~ 'Late',
                                condition == "Late" ~ 'Late',
                                condition == "Aab" ~ 'MultiAab',
                                condition == "One" ~ 'SingleAab'))

write.csv(combined_data, "combined_npod_hpap_kegg_react.csv")


