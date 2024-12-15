#demuxlet script 
#/home/.conda/envs/mamba/envs/multiome_pipeline 
export TMPDIR=Yang_demux_out
mamba activate multiome_pipeline
cd /nfs/lab/projects/
DNA=LV003

#samtools merge Yang_demux_out/$Pool.bam Yang_cellranger_atac/$DNA/outs/possorted_bam.bam Yang_cellranger_rna/$RNA/outs/possorted_genome_bam.bam
awk '{ print $1 }' good_barcodes/${DNA}_PF_cells.txt > good_barcodes/${DNA}_cells.txt
popscle demuxlet --sam cellranger.bams/dna/${DNA}.possorted_bam.bam --vcf ../genotypes/imputed/Liver_PairedTag_Pool1_All_snps.ordered.nocontig.vcf --field GT --out demux_out/$Pool --group-list good_barcodes/${DNA}_cells.txt
