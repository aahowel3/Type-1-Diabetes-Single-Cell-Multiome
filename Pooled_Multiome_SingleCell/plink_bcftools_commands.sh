#/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/250128_nPOD_redownload/Plink/preprocess_variants_forimputation_redo.R
#https://samtools.github.io/bcftools/howtos/plugin.fixref.html
#rscript write out rsids
#rm header selected rsids.csv
plink --bfile nPOD_UFDIchip_renamed --extract  selected_rsIDS.csv --keep-allele-order --make-bed --out nPOD_UFDIchip_renamed_filtered
#in rscript rewrite filtered bim file with old pos to new pos
#just overwrite with same name or move old file name to _old
#can use "check" variable in Rscript to lookup 2000ish variants where pos changes
#rm header nPOD_UFDIchip_renamed_filtered.bim
plink --bfile nPOD_UFDIchip_renamed_filtered --make-bed --out nPOD_UFDIchip_renamed_filtered_hg38
bcftools +fixref nPOD_UFDIchip_renamed_filtered_hg38.vcf.gz -Ob -o nPOD_UFDIchip_renamed_filtered_hg38_top.vcf.gz -- -f /tscc/nfs/home/abhowell/cistarget_genomes/hg38.fa -m top
plink --vcf nPOD_UFDIchip_renamed_filtered_hg38_top.vcf.gz --make-bed --out nPOD_UFDIchip_renamed_filtered_hg38_top

#here we are back into step 4 of emilys pipeline
#https://www.notion.so/Demuxlet-Genotype-array-to-vcfs-c16ac0bc3fd34433a81622948d89547f
plink --bfile nPOD_UFDIchip_renamed_filtered_hg38_top --geno 0.05 --maf 0.01 --hwe 1e-5 -keep-allele-order --make-bed --out nPOD_UFDIchip_renamed_filtered_hg38_top_mafgeno
#then step 3 of emilys pipeline
plink --freq --bfile nPOD_UFDIchip_renamed_filtered_hg38_top_mafgeno --keep-allele-order --out nPOD_UFDIchip_renamed_filtered_hg38_top_mafgeno
#run McCarthy tools
/tscc/nfs/home/abhowell/imputeInversion/HRC-1000G-check-bim.pl -b nPOD_UFDIchip_renamed_filtered_hg38_top_mafgeno.bim -f nPOD_UFDIchip_renamed_filtered_hg38_top_mafgeno.frq -r PASS.Variantsbravo-dbsnp-all.tab.gz -h
