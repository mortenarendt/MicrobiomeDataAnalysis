# load data
library(phyloseq)
library(tidyverse)
rm(list = ls())

mapping_file <- read.delim('~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/GutMicrobiome_INULIN/mapping_INULIN+SCFAs.txt', header=TRUE)


phyX <- import_biom(BIOMfilename = '~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/GutMicrobiome_INULIN/OTU_table.gg_INULIN.biom', 
                      treefilename = '~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/GutMicrobiome_INULIN/midpoint.tre')

smp_dt <- sample_data(mapping_file)
rownames(smp_dt) <- smp_dt$X.SampleID %>% trimws()
phyX <- merge_phyloseq(phyX, smp_dt)
save(file = '~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/Rats_inulin.RData', list = c('phyX'))
