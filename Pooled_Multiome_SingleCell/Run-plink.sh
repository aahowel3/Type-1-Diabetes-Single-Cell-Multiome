plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4 --exclude ./Exclude-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --out ./TEMP1
plink --bfile ./TEMP1 --update-map ./Chromosome-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --update-chr --make-bed --out ./TEMP2
plink --bfile ./TEMP2 --update-map ./Position-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --out ./TEMP3
plink --bfile ./TEMP3 --flip ./Strand-Flip-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --out ./TEMP4
plink --bfile ./TEMP4 --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --out ./nPOD_UFDIchip_renamed_filtered_step4-updated
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 1 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr1
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 2 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr2
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 3 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr3
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 4 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr4
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 5 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr5
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 6 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr6
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 7 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr7
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 8 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr8
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 9 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr9
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 10 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr10
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 11 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr11
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 12 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr12
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 13 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr13
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 14 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr14
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 15 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr15
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 16 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr16
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 17 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr17
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 18 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr18
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 19 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr19
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 20 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr20
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 21 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr21
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 22 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr22
plink --bfile ./nPOD_UFDIchip_renamed_filtered_step4-updated --reference-allele ./Force-Allele1-nPOD_UFDIchip_renamed_filtered_step4-HRC.txt --make-bed --chr 23 --out ./nPOD_UFDIchip_renamed_filtered_step4-updated-chr23
rm ./TEMP*
