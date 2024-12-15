#!/bin/bash
#has to be in conda evo get_data

while read -r i; do
pysradb gse-to-gsm $i
done < gse_list_noFluidigm_pickup.txt >> gsm_list_raw.txt

cp gsm_list_raw.txt gsm_list.txt

sed -i '/study_alias/d' gsm_list.txt
sed -i '/No results/d' gsm_list.txt

#count how many gses you started with
wc -l gse_list_noFluidigm_pickup.txt
#count how many were retruned - f1 is a mix of gsms and ssrps 
#even though you asked it for gsm - dont worry about it
cut -f1 gsm_list.txt | uniq | wc -l
#return missing gses with
grep "No" gsm_list_raw.txt
#the missing ones should all be private or no sra or unpublished - otherwise theres a problem with pysradb

for i in $(cut -f3 gsm_list.txt) 
do
pysradb srx-to-srr $i
done >> srrs_withmetadata.txt

#check which SRX's did not return results
#SRX and SRR are not 1:1 more SRR than SRX
cut -f1 srrs_withmetadata.txt | sort | uniq | wc -l
cut -f3 gsm_list.txt | uniq | wc -l
diff <(cut -f1 srrs_withmetadata.txt | sort | uniq) <(cut -f3 gsm_list.txt | sort | uniq)

#SRP270320 this srp is affiliated with SRRs ending in 3328-3318 - these entries col orders are flipped, create a seperate awk command for these 
#1218a1217
#SRX2733166 this SRX will return an SRR in the correct order, just glitched, append to srr_metadata
#1616a1616,1
#SRX8673318
#SRX8673319
#SRX8673320
#SRX8673321
#SRX8673322
#SRX8673323
#SRX8673324
#SRX8673325
#SRX8673326
#SRX8673327
#SRX8673328

#append missing entry 
cp srrs_withmetadata.txt srrs_withmetadata_added.txt
pysradb srx-to-srr SRX2733166 >> srrs_withmetadata_added.txt
awk -F $'\t' '$11 ~ "TRANSCRIPTOMIC SINGLE CELL" {print $1, $2}' srrs_withmetadata_added.txt > srrs_rnaseqonly.txt
#for atac seq datasets column 11 will not be helpful just says "GENOMIC" easiest to search for ATAC
grep "ATAC" srrs_withmetadata_added.txt | cut -f1,2 > srrs_ATAConly.txt


#grab entries that were flipped and use diff awk column selection
#not just flipped but columns only include SRP and SRX, not SRR
#actually..... when you look it up ..... there are 4,000 SRRs assoicated with the SRP here
#not adding to final call file

#add 1-length config number for fastq-dump and column headers
awk '{print $0 "\t" NR }' srrs_rnaseqonly.txt > srrs_rnaseqonly.config

#cross ref gse to srr to put SRRs in folders - unifying ID is the SRX
awk 'NR==FNR{a[$3]=$1; next}{$1=a[$1]; print}' gsm_list.txt srrs_rnaseqonly.config > srrs_rnaseqonly_crossref.config

#manually add header to config to hard to figure out - space seperated
#Var1="SRPID","SRRID","ArrayTaskID"

#adding a size check to files could be helpful 
for i in $(awk '{print $2}' srrs_rnaseqonly_crossref.config)
do
vdb-dump $i --info | grep "size"
done >> sizes.txt

#combine files together - look for a subset tester that wont take a day to run
#60,000,000,000 = 60GB
paste file1 file2 | column -s $'\t' -t

