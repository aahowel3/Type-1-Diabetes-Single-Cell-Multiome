# Type 1 Diabetes Disease Progression Single-cell Multiomics and Gene Regulatory Elements 

## Clustering Pipeline Comparisons 

Human Pancreas Analysis Program (HPAP): 292 donors purified islets scRNA-seq 10X-Chromium-GEX-3p-v2. Pipeline 1 uses starsolo for alignment, samtools for duplicate removal, cellbender for ambient RNA removal. Secondary script `cellbender_clusters.R` used for clustering annotation. Pipeline 2 uses salmon for alignment, scanpy for clustering, trajectory inference, differential expression. 

![alt text](https://github.com/aahowel3/Type-1-Diabetes-Single-Cell-Multiome/blob/main/pipeline_comparisons.png)

## UCell Gene Enrichment Scores 

Network for Pancreatic Organ donors with Diabetes: 64 donors whole pancreas endocrine and exocrine tissues, scRNA-seq, ATAC-seq, single sample multiome, and pooled multiome data. 
