# Downloading TopMed imputation results
## Unzip and check files
## Replace qvMG1IL@Gq0gj with the password in the topmed email
for file in *.zip; do unzip -P qvMG1IL@Gq0gj $file; done
md5sum *.zip

## Check proportion of 0/0, 0/1, 1/0, 1/1 calls in large Florida array v. Florida extra
bcftools query -l chr10.dose.vcf.gz > chr10_samples.txt
bcftools stats -S chr10_samples.txt chr10.dose.vcf.gz > chr10_bcftools_stats.txt
