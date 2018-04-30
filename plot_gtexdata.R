# preparing data
library(dplyr)
library(tidyr)
library(ggplot2)
samples = read.delim('GTEx_Data_V6_Annotations_SampleAttributesDS.txt', sep='\t') %>%
	select(SAMPID, primary.tissue=SMTS, tissue=SMTSD)
gtex = read.table('GTEx_Analysis_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct', skip=2, colClasses=c('character', 'character', rep('numeric', 2921)), stringsAsFactors=F, header=T)

gtex.bysample = expdata %>% 
	gather(SAMPID, expression, -Name, -Description) %>% 
	mutate(SAMPID=gsub('\\.', '-', SAMPID)) %>% 
	left_join(samples)


arhgef17 = gtex.bysample %>% filter(Description=='RUNX3')
arhgef17 %>%
	ggplot(aes(x=primary.tissue, y=expression)) + geom_boxplot()
