library(vroom)
library(stringr)
library(dplyr)
library(tidyr)

vcf.raw <- vroom('/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/imputed/Florida_array_extra/chr10.dose.r209.vcf.gz',delim='\t', comment='##')
#dim(vcf.raw)
#head(vcf.raw)

vcf.df <- vcf.raw
#vcf.df <- tibble::column_to_rownames(vcf.df, 'ID')
vcf.df <- select(vcf.df, -`#CHROM`, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)

extract.genotypes <- function(x) {
  base.states <- c(`0|0`=0,`0|1`=1,`1|0`=1,`1|1`=2)
  alleles <- str_extract(x, '^[01]\\|[01]')
  gene <- base.states[alleles]
  return(gene)
}

vcf.df <- mutate_all(vcf.df, extract.genotypes)

library(purrr)
vcf_table = vcf.df %>% map(table)

sink("/tscc/projects/ps-gaultonlab/projects/mega_pancreas/data/genotypes/nPOD/imputed/Florida_array_extra/chr10.r209.dose.genotypecount.txt")
print(vcf_table)
sink()
