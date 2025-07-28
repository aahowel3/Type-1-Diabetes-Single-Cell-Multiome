setwd("/tscc/projects/ps-gaultonlab/abhowell/scenic_hic/")

library(rtracklayer)
library(GenomicRanges)
library(R.utils)

input_file <- "NIHMS1834333-supplement-MMC2.csv"
df <- read.csv(input_file)

# === Step 2: Create GRanges for chr_a/start_a/end_a and chr_b/start_b/end_b ===
gr_a <- GRanges(seqnames = df$chr_a, ranges = IRanges(start = df$start_a, end = df$end_a))
gr_b <- GRanges(seqnames = df$chr_b, ranges = IRanges(start = df$start_b, end = df$end_b))

# === Step 3: Download and import UCSC liftover chain file (hg19 → hg38) ===
chain_file <- "hg19ToHg38.over.chain"
if (!file.exists(chain_file)) {
  download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "chain.gz")
  gunzip("chain.gz", destname = chain_file)
}
chain <- import.chain(chain_file)

# === Step 4: Perform liftover ===
lifted_a <- liftOver(gr_a, chain)
lifted_b <- liftOver(gr_b, chain)

# === Step 5: Filter for rows with successful 1-to-1 mapping in both A and B ===
mapped_a <- sapply(lifted_a, length) == 1
mapped_b <- sapply(lifted_b, length) == 1
mapped_both <- mapped_a & mapped_b

lifted_a_clean <- unlist(lifted_a[mapped_both])
lifted_b_clean <- unlist(lifted_b[mapped_both])
df_clean <- df[mapped_both, ]

# === Step 6: Replace original columns with lifted coordinates ===
df_clean$chr_a <- as.character(seqnames(lifted_a_clean))
df_clean$start_a <- start(lifted_a_clean)
df_clean$end_a <- end(lifted_a_clean)

df_clean$chr_b <- as.character(seqnames(lifted_b_clean))
df_clean$start_b <- start(lifted_b_clean)
df_clean$end_b <- end(lifted_b_clean)

# === Step 7: Save output CSV ===
output_file <- "liftover_chr_a_and_chr_b_to_hg38.csv"
write.csv(df_clean, output_file, row.names = FALSE)

cat("Liftover complete. Output saved to:", output_file, "\n")

#spotchecked, I trust it
#https://genebe.net/tools/liftover
#to sort GRN gene-region connections to these 3 cell types 
#will need to add in RSS scores 
library(data.table)
library(tidyverse)
df2 = read.csv("eRegulon_direct_npod2.csv", sep =",", header = TRUE)
df1 = fread("npod2_celltype_gene.csv", sep =",", header = TRUE)

#add RSSs in so you can pull out gene-region relationships relevant to the celltypes in HiC
#acinar, alpha,beta
# Step 1: Melt the expression dataframe
df1_long <- df1 %>%
  pivot_longer(
    cols = -V1,
    names_to = "Gene_signature_name",
    values_to = "Expression"
  ) %>%
  rename(Celltype = V1)

# Step 2: Pivot wider so each Gene_signature_name has one row with each cell type as a column
df1_wide <- df1_long %>%
  pivot_wider(
    names_from = Celltype,
    values_from = Expression
  )

# Step 3: Rename the cell type columns to include "_gene"
# (Keep the Gene_signature_name column as-is)
df1_wide <- df1_wide %>%
  rename_with(.fn = ~ paste0(.x, "_gene"), .cols = -Gene_signature_name)

# Step 4: Merge with the second dataframe
merged_df <- left_join(df2, df1_wide, by = "Gene_signature_name")

#repeat for region scores
df1 = fread("npod2_celltype_region.csv", sep =",", header = TRUE)
# Step 1: Melt the expression dataframe
df1_long <- df1 %>%
  pivot_longer(
    cols = -V1,
    names_to = "Region_signature_name",
    values_to = "Expression"
  ) %>%
  rename(Celltype = V1)

# Step 2: Pivot wider so each Gene_signature_name has one row with each cell type as a column
df1_wide <- df1_long %>%
  pivot_wider(
    names_from = Celltype,
    values_from = Expression
  )

# Step 3: Rename the cell type columns to include "_region"
# (Keep the Gene_signature_name column as-is)
df1_wide <- df1_wide %>%
  rename_with(.fn = ~ paste0(.x, "_region"), .cols = -Region_signature_name)

npod2 <- left_join(merged_df, df1_wide, by = "Region_signature_name")

#hello new issue! You’re evaluating the region-gene intersection with HiC data, 
#but you’re double counting region-genes if theres the same region-gene controlled by a diff TF
npod2$region_gene = paste(npod2$Region, npod2$Gene)
npod2 <- npod2 %>%
  distinct(region_gene, .keep_all = TRUE)


#now were good to convert the gene names to a chromosme region 
library(biomaRt)
my_genes = npod2$Gene
length(unique(my_genes))

m <- useMart('ensembl', dataset='hsapiens_gene_ensembl') # create a mart object
df <- getBM(mart=m, attributes=c('hgnc_symbol', 'description', 'chromosome_name',
                                 'start_position', 'end_position', 'strand',
                                 'ensembl_gene_id'),
            filters='hgnc_symbol', values=my_genes)


#some of these are duplicate entries of the same chr, gene, id,drop
df <- df[!duplicated(df), ]

valid_chromosomes <- c(as.character(1:22), "X", "Y")
# Convert column to character if it's not already
df$chromosome_name <- as.character(df$chromosome_name)
# Filter the rows
df <- df[df$chromosome_name %in% valid_chromosomes, ]


#after all that these are your problem childrne 
dup_symbols <- df$hgnc_symbol[duplicated(df$hgnc_symbol) | duplicated(df$hgnc_symbol, fromLast = TRUE)]
# Filter the dataframe to keep only duplicated rows
df_duplicates <- df[df$hgnc_symbol %in% dup_symbols, ]
# View result
df_duplicates

#unique are 5 6262-5 = 6257 i.e the number of unique genes in the npod2 GRN
#this will add the chr coordinates to npod2 - duplicating the rows where a gene has multiple matching sites from ensembl
colnames(df)[colnames(df) == "hgnc_symbol"] <- "Gene"
# Merge npod2 with df by Gene (many-to-one or many-to-many join)
npod2_expanded <- merge(npod2, df[, c("Gene", "chromosome_name", "start_position", "end_position")], by = "Gene")

library(GenomicRanges)
library(dplyr)
library(tidyverse)
###
#look at the intersection of scneic gene-region and hic region-region
#first need to convert chr:X-Y into the same comparable format
# Step 1: Parse Region from npod2_expanded
npod2_expanded <- npod2_expanded %>%
  mutate(region_chr = str_extract(Region, "(?<=chr)\\w+"),
         region_start = as.numeric(str_extract(Region, "(?<=:)[0-9]+")),
         region_end = as.numeric(str_extract(Region, "(?<=-)\\d+")))

# Step 2: Create GRanges for npod2_expanded (both region-based and gene-locus-based)
gr_region <- GRanges(
  seqnames = npod2_expanded$region_chr,
  ranges = IRanges(start = npod2_expanded$region_start, end = npod2_expanded$region_end),
  npod_index = 1:nrow(npod2_expanded),
  source = "region"
)

gr_gene <- GRanges(
  seqnames = npod2_expanded$chromosome_name,
  ranges = IRanges(start = npod2_expanded$start_position, end = npod2_expanded$end_position),
  npod_index = 1:nrow(npod2_expanded),
  source = "gene"
)

# Step 3: Create GRanges for df_clean's chr_a and chr_b (with ±5kb buffer)
df_clean_expanded <- df_clean %>%
  mutate(chr_a = gsub("chr", "", chr_a),
         chr_b = gsub("chr", "", chr_b))

gr_a <- GRanges(
  seqnames = df_clean_expanded$chr_a,
  ranges = IRanges(start = df_clean_expanded$start_a - 500, end = df_clean_expanded$end_a + 500),
  df_index = 1:nrow(df_clean_expanded),
  anchor = "A"
)

gr_b <- GRanges(
  seqnames = df_clean_expanded$chr_b,
  ranges = IRanges(start = df_clean_expanded$start_b - 500, end = df_clean_expanded$end_b + 500),
  df_index = 1:nrow(df_clean_expanded),
  anchor = "B"
)

# Step 4: Find overlaps
hits_a_region <- findOverlaps(gr_region, gr_a)
hits_b_region <- findOverlaps(gr_region, gr_b)
hits_a_gene   <- findOverlaps(gr_gene, gr_a)
hits_b_gene   <- findOverlaps(gr_gene, gr_b)

# Step 5: Combine and deduplicate hits
get_overlap_df <- function(hits, label) {
  data.frame(
    npod_index = queryHits(hits),
    df_index = subjectHits(hits),
    anchor = label
  )
}

overlaps <- bind_rows(
  get_overlap_df(hits_a_region, "A"),
  get_overlap_df(hits_b_region, "B"),
  get_overlap_df(hits_a_gene, "A"),
  get_overlap_df(hits_b_gene, "B")
) %>%
  distinct() %>%
  group_by(npod_index, df_index) %>%
  summarise(overlap_type = paste(sort(unique(anchor)), collapse = "+"), .groups = "drop")

# Step 6: Merge overlapping hits
final_overlap_df <- overlaps %>%
  left_join(npod2_expanded %>% mutate(npod_index = row_number()), by = "npod_index") %>%
  left_join(df_clean_expanded %>% mutate(df_index = row_number()), by = "df_index")

# Step 7: Identify non-overlapping npod2 rows
all_npod_indices <- 1:nrow(npod2_expanded)
overlapping_npod_indices <- unique(overlaps$npod_index)
non_overlapping_indices <- setdiff(all_npod_indices, overlapping_npod_indices)

# Step 8: Create non-overlapping records
non_overlap_df <- npod2_expanded[non_overlapping_indices, ] %>%
  mutate(npod_index = non_overlapping_indices,
         overlap_type = "none")

# Get placeholder NA columns for df_clean_expanded (same structure, repeated to match rows)
na_df_clean <- df_clean_expanded[1, ] %>% slice(0)  # keep structure, zero rows
na_columns <- as.data.frame(matrix(NA, nrow = nrow(non_overlap_df), ncol = ncol(na_df_clean)))
colnames(na_columns) <- colnames(na_df_clean)

# Combine npod2 + NA df_clean
non_overlap_df <- bind_cols(non_overlap_df, na_columns)

# Step 9: Combine both
final_overlap_df_full <- bind_rows(final_overlap_df, non_overlap_df)


###because HiC is only on alpha beta acinar - filter by RSS score for fair comparison 
#Now remove connectiosn thar are not part of alpha,beta, acinar GRNs by RSS score
score_cols <- c("acinar_gene", "acinar_region", "beta_gene", "beta_region", "alpha_gene", "alpha_region")

# Combine unique top 10 eRegulon_names from each column
top_eregulons <- score_cols %>%
  map(~ npod2 %>%
        group_by(eRegulon_name) %>%
        slice_max(.data[[.x]], n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(desc(.data[[.x]])) %>%
        slice_head(n = 10) %>%
        pull(eRegulon_name)) %>%
  unlist() %>%
  unique()

# View result
top_eregulons
final_final_overlap_df = final_overlap_df_full[final_overlap_df_full$eRegulon_name %in% top_eregulons,]

#how many resolustion at what bins?
final_final_overlap_df$type_resolution = paste(final_final_overlap_df$overlap_type, final_final_overlap_df$resol)

plot_data <- final_final_overlap_df %>%
  separate(type_resolution, into = c("type", "resolution"), sep = " ") %>%
  mutate(resolution = factor(resolution, levels = c("1000", "2000", "4000")))

# Step 3: HiC baseline dataset
hic_baseline <- data.frame(
  type = "HiC_baseline",
  resolution = factor(c("1000", "2000", "4000"), levels = c("1000", "2000", "4000")),
  count = c(2, 104, 10728)
)

# Step 4: Compute counts from plot_data and bind with baseline
plot_counts <- plot_data %>%
  count(type, resolution, name = "count")

combined_data <- bind_rows(plot_counts, hic_baseline)

# Step 5: Plot as side-by-side bars with counts and log scale
ggplot(combined_data, aes(x = resolution, y = count, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.9),
            vjust = -0.2, size = 3) +
  scale_y_log10() +
  scale_fill_manual(values = c(
    "A" = "#1f77b4",
    "A+B" = "#2ca02c",
    "B" = "#ff7f0e",
    "none" = "#bbbbbb",
    "HiC_baseline" = "gray50"
  )) +
  labs(title = "Overlap Type by Resolution (log scale) — nPOD2 vs HiC Baseline",
       x = "Resolution (bp)",
       y = "Log10 Count",
       fill = "Overlap Type") +
  theme_minimal()



#side plot - see if there's something distinguishing about scores of nonoverlapping 
#vs those overlapping with HiC data
library(ggplot2)
ggplot(final_final_overlap_df, aes(x = overlap_type, y = abs(rho_R2G), fill = overlap_type)) +
  geom_boxplot() +
  labs(title = "Importance R2G by Overlap Type",
       x = "Overlap Type",
       y = "Importance R2G") +
  theme_minimal()


