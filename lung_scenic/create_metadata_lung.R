library(dplyr)
library(stringr)
library(tibble)
library(Seurat)
library(dplyr)
library(Signac)
library(tibble)
library(Matrix)
setwd("/scratch/aahowel3/lung_scenic/")

hpap_full <- readRDS("92samp_UnionPeaks_UpdateFragments_newPeaks_removeDups_600k.rds")
#make sure these are in chr- notation
head(rownames(hpap_full@assays$Final_Peaks$counts))

meta = hpap_full@meta.data
table(meta$new.final.celltype)

length(table(meta$donor_demux))
meta$celltypebydonor = paste0(meta$new.final.celltype,meta$donor_demux)
#make sure to add this column to Seurat object so cells are parsable by combined donor+celltype
hpap_full@meta.data$celltypebydonor = paste0(meta$new.final.celltype,meta$donor_demux)
#output table to get downsampling schema
celltypebydonor_table = table(meta$celltypebydonor)
celltypebydonor_table = as.data.frame(celltypebydonor_table)

dat <- celltypebydonor_table
# list of known cell types from your screenshot
celltypes <- c(
  "Aberrant.Basaloid", "AT.trans", "Basal", "Goblet", "Multiciliated",
  "Platelet", "VSMC", "AEC", "AT0", "CAP1", "IM", "Myeloid.p", "PVEC",
  "AF1", "AT1", "CAP2", "LEC", "NK", "RAS",
  "AF2", "AT2", "AT2.p", "Club", "Mast", "Peribronchial.FB", "SCMF",
  "AM", "DC", "Mesothelium", "Pericyte", "SVEC",
  "ASMC", "B.Cells", "Fibrotic.FB", "MON", "Plasma", "T.Cells"
)

# Build regex pattern of valid celltypes
pattern <- paste0("^(", paste(celltypes, collapse="|"), ")")
# Extract celltype and donor
dat2 <- dat %>%
  mutate(
    celltype = str_extract(Var1, pattern),
    donor    = str_remove(Var1, pattern)
  )
write.csv(dat2,"celltypebydonor_table2.csv")


#downsample scehma
schema <- read.csv("downsample_schema.csv", stringsAsFactors = FALSE)
# Convert to named vector
celltype_counts <- setNames(schema$new_count, schema$celltypebydonor)


meta <- hpap_full@meta.data %>%
  rownames_to_column("cell_name")

sampled_cells <- meta %>%
  filter(celltypebydonor %in% names(celltype_counts)) %>%
  group_by(celltypebydonor) %>%
  group_map(~ {
    n_cells <- celltype_counts[[.y$celltypebydonor]]
    slice_sample(.x, n = min(n_cells, nrow(.x)))
  }) %>%
  bind_rows() %>%
  pull(cell_name)

# 5. Subset the Seurat object
hpap2 <- subset(hpap_full, cells = sampled_cells)

# 6. Confirm result
hpap_full = hpap2
hpap = hpap_full

hpap[["RNA"]]$counts <- as(object = hpap[["RNA"]]$counts, Class = "dgCMatrix")
hpap[["Final_Peaks"]]$counts <- as(object = hpap[["Final_Peaks"]]$counts, Class = "dgCMatrix")

writeMM(hpap[["RNA"]]$counts, "lung_RNA_counts.mtx")
writeMM(hpap[["Final_Peaks"]]$counts, "lung_atac_counts.mtx")

meta = hpap@meta.data
write.table(
  meta, "lung_meta_data.tsv",
  col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(
  colnames(hpap), "lung_barcodes.tsv",
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(
  rownames(hpap[["RNA"]]), "lung_gene_names.tsv",
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(
  rownames(hpap[["Final_Peaks"]]), "lung_region_names.tsv",
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

######
##downsampling schema
######
library(dplyr)

### -----------------------------
### Load original table
### -----------------------------
df <- read.csv("celltypebydonor_table2.csv", stringsAsFactors = FALSE)

### df must contain:
### Var1, celltype, donor, Freq

### -----------------------------
### Ranking-preserved target totals (from your final schema)
### -----------------------------
target_totals <- c(
  "AT2"=8168, "AT1"=4690, "AF1"=2772, "AM"=2771, "AT.trans"=2770,
  "Multiciliated"=2769, "CAP1"=2768, "T.Cells"=2767, "RAS"=2766,
  "MON"=2765, "IM"=2764, "CAP2"=2763, "Goblet"=2762, "NK"=2761,
  "AF2"=2760, "Fibrotic.FB"=2760, "LEC"=2758, "AEC"=2654, "AT0"=2302,
  "Club"=2172, "Peribronchial.FB"=2136, "PVEC"=2057, "SVEC"=1859,
  "SCMF"=1660, "DC"=1638, "Pericyte"=1209, "Basal"=1206, "VSMC"=1200,
  "B.Cells"=1121, "Plasma"=1120, "Mast"=1075, "Aberrant.Basaloid"=751,
  "Platelet"=556, "Mesothelium"=469, "ASMC"=333, "Myeloid.p"=150
)

### -----------------------------
### Compute original totals (sorted descending)
### -----------------------------
orig_totals <- df %>%
  group_by(celltype) %>%
  summarize(original_total = sum(Freq)) %>%
  arrange(desc(original_total))

group_order <- orig_totals$celltype
alpha <- 0.33  # cube-root smoothing exponent

df$new <- 0
prev_group_total <- 1e12   # large initial value

### -----------------------------
### MAIN LOOP: apply smoothing, capacity capping, ranking preservation
### -----------------------------
for (g in group_order) {
  
  sub_idx <- which(df$celltype == g)
  sub <- df[sub_idx, ]
  
  orig_total <- sum(sub$Freq)
  
  # Determine group target
  target <- ifelse(g %in% names(target_totals), target_totals[g], orig_total)
  
  # Cannot exceed original total
  target <- min(target, orig_total)
  
  # Ranking constraint: must not exceed previous group
  target <- min(target, prev_group_total)
  
  old <- sub$Freq
  
  if (sum(old) == 0) {
    df$new[sub_idx] <- 0
    prev_group_total <- 0
    next
  }
  
  ### Smoothing weights
  w <- old^alpha
  if (sum(w) == 0) w <- rep(1, length(old))
  
  raw_new <- w / sum(w) * target
  
  assigned <- floor(raw_new)
  
  ### Capacity enforcement: new ≤ original
  assigned <- pmin(assigned, old)
  
  ### Minimum of 1 for donors with original > 0
  assigned[old > 0 & assigned < 1] <- 1
  
  ### Adjust totals
  diff <- target - sum(assigned)
  
  # If we need to ADD counts
  if (diff > 0) {
    capacity <- old - assigned
    ord <- order(capacity, decreasing = TRUE)
    for (i in ord) {
      if (diff == 0) break
      cap <- capacity[i]
      if (cap > 0) {
        add_val <- min(cap, diff)
        assigned[i] <- assigned[i] + add_val
        diff <- diff - add_val
      }
    }
  }
  
  # If we need to REMOVE counts
  if (diff < 0) {
    remov <- -diff
    ord <- order(assigned, decreasing = TRUE)
    for (i in ord) {
      if (remov == 0) break
      max_rem <- assigned[i] - (ifelse(old[i] > 0, 1, 0))
      if (max_rem > 0) {
        dec <- min(max_rem, remov)
        assigned[i] <- assigned[i] - dec
        remov <- remov - dec
      }
    }
  }
  
  ### If still cannot meet target (capacity exhausted) — Option A:
  if (sum(assigned) > prev_group_total) {
    excess <- sum(assigned) - prev_group_total
    ord <- order(assigned, decreasing = TRUE)
    for (i in ord) {
      if (excess == 0) break
      max_rem <- assigned[i] - (ifelse(old[i] > 0, 1, 0))
      if (max_rem > 0) {
        dec <- min(max_rem, excess)
        assigned[i] <- assigned[i] - dec
        excess <- excess - dec
      }
    }
  }
  
  ### Save
  df$new[sub_idx] <- assigned
  prev_group_total <- sum(assigned)
}

### -----------------------------
### Write before/after output
### -----------------------------
write.csv(
  df %>% select(Var1, celltype, donor, before = Freq, after = new),
  "before_after_perdonor.csv",
  row.names = FALSE
)

write.table(
  paste0("\"", df$Var1, "\" = ", df$new),
  file = "celltype_counts_Rvector.txt",
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
