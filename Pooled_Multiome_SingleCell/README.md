# Convert Genotype Array to VCF via plink
https://www.notion.so/Demuxlet-Genotype-array-to-vcfs-c16ac0bc3fd34433a81622948d89547f

# Downloading TopMed imputation results
## Unzip and check files
## Replace qvMG1IL@Gq0gj with the password in the topmed email
for file in *.zip; do unzip -P qvMG1IL@Gq0gj $file; done
md5sum *.zip

## Check proportion of 0/0, 0/1, 1/0, 1/1 calls in large Florida array v. Florida extra
bcftools query -l chr10.dose.vcf.gz > chr10_samples.txt
bcftools stats -S chr10_samples.txt chr10.dose.vcf.gz > chr10_bcftools_stats.txt

# File locations on NFS
NPOD2 processing: /nfs/lab/welison/mega_pancreas/notebooks/merge_objects/250324_nPOD_Pool_Merge.ipynb <br /> 
NPOD2 Clustering and Annotation: /nfs/lab/welison/mega_pancreas/notebooks/single_sample_processing/241125_AH_nPOD_Single_Sample_Processing_Pool.ipynb
/home/ahowell <br /> 
Spleen Clustering and Annotation: /home/ahowell/Spleen_celltype_annotation_AH_4_4_25.ipynb <br /> 
NPOD2 + Spleen + PLN Integration with Harmony: /nfs/lab/rlmelton/npod/notebooks/postdoc/pilotIntegration_Hillbloom_PLN_PancGRS.ipynb <br />

# Genotyping off ATAC
https://www.notion.so/Variant-Calling-122451be7c2d80c39987d702568d8ea9 

# Michigan HLA Imputation
https://imputationserver.sph.umich.edu/#! 
Log in to Michigan HLA and upload chr 6 file with following selections
Reference Panel: Four-digit Multi-Ethnic HLA Reference Panel (GRCh37/hg19)
Phasing: Eagle Phased
Mode: Quality Control and Imputation
Michigan will send you an email when the imputation results are ready
Follow the link in your email and find the curl command under imputation results and link for downloading all files
