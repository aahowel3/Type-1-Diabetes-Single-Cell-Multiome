#in tmux window demuxlet 
#Filter variants based on imputation quality and minor allele frequency
for chr in {1..22}; do
bcftools view chr${chr}.dose.vcf.gz -i 'R2>0.9 & MAF>0.01' -Oz -o chr${chr}.dose.r209.maf01.vcf.gz
done

# Merge chromosomes
bcftools concat chr*.dose.r209.maf01.vcf.gz -Oz -o merged.dose.r209.maf01.vcf.gz
# Change labels from rsID to position
bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' merged.dose.r209.maf01.vcf.gz -Oz -o PANC.89merged.dose.r209.maf01.snppos.vcf.gz

# Rename samples - Make rename.merged.dose.r209.maf01.vcf.gz.txt first
bcftools reheader PANC.89merged.dose.r209.maf01.snppos.vcf.gz -s rename.merged.dose.r209.maf01.vcf.gz.txt -Oz -o PANC,89merged.dose.r209.maf01.names.fixed.vcf.gz 
tabix PANC.89merged.dose.r209.maf01.names.fixed.vcf.gz

# Subset to samples in specific pool and overlapping catlas peaks
bcftools view PANC.89merged.dose.r209.maf01.names.fixed.vcf.gz -R /nfs/lab/cmcgrail/nonDR3DR4_annotations/single_cell_atlas_peaks/raw/raw_filetered_peaks/Catlas_merged_hg38.bed -Ov -o PANC_Pool1_All_snps.vcf
