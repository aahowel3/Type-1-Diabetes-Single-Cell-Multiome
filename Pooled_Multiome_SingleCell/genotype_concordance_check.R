library(vcfR)
#vcf_floridaog = read.vcfR("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/imputed/Florida_array_og/floridaog_concordancetester.vcf")
vcf_floridaog = read.vcfR("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/imputed/Florida_array_og/npod2.merged.dose.r209.snppos.vcf.gz")

#right do this pretop med as well - are there issues with the imputation or the liftover from the start
#vcf_floridaog = read.vcfR("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/most.raw/UFDI_array_dbsnp/merged.pretopmed.vcf.gz")

vcfr = as.data.frame(vcf_floridaog@fix)
vcfq = as.data.frame(vcf_floridaog@gt)
vcf_floridaog2 = cbind(vcfr, vcfq)

#Only need this if doing a pre-imuted vcf
#vcf_floridaog2$chrompos = paste(vcf_floridaog2$CHROM,vcf_floridaog2$POS)

#vcf_floridaog = read.vcfR("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/imputed/npod2_in_house/npod2_inhouse_concordancetester.vcf")
vcf_floridaog = read.vcfR("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/imputed/npod2_in_house/npod2.merged.dose.r209.snppos.vcf.gz")
#check pretopmed as well
#vcf_floridaog = read.vcfR("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/most.raw/PLINK_210623_0733/merged.pretopmed.vcf.gz")

vcfr = as.data.frame(vcf_floridaog@fix)
vcfq = as.data.frame(vcf_floridaog@gt)
vcf_npod2inhouse = cbind(vcfr, vcfq)

#only need to add if doing pre imputed vcf
#vcf_npod2inhouse$chrompos = paste(vcf_npod2inhouse$CHROM,vcf_npod2inhouse$POS)

#if preimputed its by chrompos not ID
df = vcf_floridaog2 %>% 
  inner_join(vcf_npod2inhouse, by="ID")

print(nrow(vcf_floridaog2))
print(nrow(vcf_npod2inhouse))
print(nrow(df))

cols <- grep("6301", names(df), value = TRUE)
df_price <- df[, c("chrompos", "REF.x", "ALT.x","REF.y","ALT.y",cols)]
#colnames(df_price) = c("INFO", "V1","V2")
colnames(df_price) = c("chrompos","REF.x", "ALT.x","REF.y","ALT.y","V1","V2")

#no need to split string unimputed - 
df_price = within(df_price, INFO<-data.frame(do.call('rbind', strsplit(as.character(INFO), ';', fixed=TRUE))))
df_price = within(df_price, V1<-data.frame(do.call('rbind', strsplit(as.character(V1), ':', fixed=TRUE))))
df_price = within(df_price, V2<-data.frame(do.call('rbind', strsplit(as.character(V2), ':', fixed=TRUE))))

table(df_price$V1 == df_price$V2)
df_price$check = (df_price$V1 == df_price$V2)

#this is not accounting for when the genotype is missing in this sample add in that clause
df_price$check = ifelse(is.na(df_price$check), "MISSING", df_price$check )

table(df_price$check)

df_price_dis=df_price[which(df_price$check == FALSE),]

table(df_price_dis$ALT.x == df_price_dis$ALT.y)
table(df_price_dis$REF.x == df_price_dis$REF.y)


table(df_price_dis$V1 == "0/1" | df_price_dis$V1 == "1/0") 

table(df_price_dis$V2 == "0/1" | df_price_dis$V2 == "1/0")


df_price$concord = df_price$V1$X1 == df_price$V2$X1
df_price_dis = df_price[df_price$concord == FALSE,]
df_price_con = df_price[df_price$concord == TRUE,]
table(df_price_dis$INFO$X1)
table(df_price_con$INFO$X1)

#so now if its not the atac peaks and its not the imputed genotypes 
#maybe it was something in the filtering with the vcf maddie used
#these are her npod1 samples
npod1_maddie = read.vcfR("/tscc/projects/ps-gaultonlab/welison/mega.panc/250613_npod1_qtls/genotypes/npod1.redo.caqtl.rename.nomiss.maf01.vcf")
vcfr = as.data.frame(npod1_maddie@fix)
vcfq = as.data.frame(npod1_maddie@gt)
npod1_maddie_vcf = cbind(vcfr, vcfq)

npod1_postimpute = read.vcfR("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/imputed/npod1_in_house/npod2.merged.dose.r209.snppos.vcf.gz")
vcfr = as.data.frame(npod1_postimpute@fix)
vcfq = as.data.frame(npod1_postimpute@gt)
npod1_postimpute_vcf = cbind(vcfr, vcfq)

df = npod1_postimpute_vcf %>% 
  inner_join(npod1_maddie_vcf, by="ID")

cols <- grep("6301", names(df), value = TRUE)
df_price <- df[, c("ID", "REF.x", "ALT.x","REF.y","ALT.y",cols)]
colnames(df_price) = c("ID","REF.x", "ALT.x","REF.y","ALT.y","V1","V2")
df_price = within(df_price, V1<-data.frame(do.call('rbind', strsplit(as.character(V1), ':', fixed=TRUE))))
df_price = within(df_price, V2<-data.frame(do.call('rbind', strsplit(as.character(V2), ':', fixed=TRUE))))

#change to unphased symbol if looking at npod1 stuff
#df_price$V2 = gsub("/","|",df_price$V2)

#tabulate non-concordance
df_price$check = (df_price$V1$X1 == df_price$V2$X1)
#this is not accounting for when the genotype is missing in this sample add in that clause
df_price$check = ifelse(is.na(df_price$check), "MISSING", df_price$check )

table(df_price$check)
df_price_dis=df_price[df_price$check == FALSE,]

#check the het swaps are comparable - ref/alt alleles should be the same in both arrays
table(df_price_dis$ALT.x == df_price_dis$ALT.y)
table(df_price_dis$REF.x == df_price_dis$REF.y)

#cehck how many sites are not het swaps
table((df_price_dis$V1$X1 == "0|1" | df_price_dis$V1$X1 == "1|0") & (df_price_dis$V2$X1 == "0|1" | df_price_dis$V2$X1 == "1|0")) 

#subtract good table (df_non) from larger table to view bad table
df_non = df_price_dis[(df_price_dis$V1$X1 == "0|1" | df_price_dis$V1$X1 == "1|0") & (df_price_dis$V2$X1 == "0|1" | df_price_dis$V2$X1 == "1|0"),]
bad = anti_join(df_price_dis, df_non, by = "ID")
