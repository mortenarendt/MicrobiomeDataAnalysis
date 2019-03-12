# load data
library(phyloseq)
rm(list = ls())

mapping_file <- read.delim('~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/mapping.txt', header=TRUE)

mapping_file$ID <- mapping_file$ID %>% trimws()
mapping_file$sex <- ifelse(mapping_file$X1male0fem==0,'Female','Male')

phyX <- import_biom(BIOMfilename = '~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/h_qual_otu_table.biom', 
                 treefilename = '~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/h_qual_rep_set.tre')

smp_dt <- sample_data(mapping_file)
rownames(smp_dt) <- smp_dt$X.SampleID %>% trimws()
phyX <- merge_phyloseq(phyX, smp_dt)
save(file = '~/Dropbox/Backup/MyDocumentsOnC/Course and teaching/Advanced_MicrobiomeAnalysis/MicrobiomeDataAnalysis/data/Mice_csec.RData', list = c('phyX'))


##
phyXlog <- transform_sample_counts(phyX , function(x) log(x + 1 / sum(x)))

GP.ord <- ordinate(phyXlog, "NMDS", "wunifrac")
p2 <-  plot_ordination(phyX, GP.ord, type="samples", color="Birth_mode")  + 
  geom_point() + 
  stat_ellipse() + 
  theme_classic()
print(p2)
p2 + geom_polygon(aes(fill=SampleType))  ggtitle("samples")
