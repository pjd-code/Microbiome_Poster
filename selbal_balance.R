rm(list = ls())

library(selbal)
library(tidyverse)
library(qiime2R)
#aged merged BAL Vsearch results as OTU table
Vsearch_A_B <- read_tsv("C:/Users/#/Desktop/selbal/ageBalTop100.tsv.txt", col_names = TRUE)
#aged merged BAL taxonomy
Tax_A_B <- read_tsv("C:/Users/#/Desktop/selbal/neat_taxonomy.txt", col_names = TRUE)
#aged merged BAL metadata
Meta_A_B <- read_tsv("C:/Users/#/Desktop/selbal/aged_merged_bal_20feb19.txt", col_names = TRUE)

###for selbal###
Meta_A_B_trans <-(Meta_A_B %>% as.data.frame() %>% column_to_rownames("SampleID"))

trans_otu_A_B <- Vsearch_A_B %>%
  gather(SampleID, valname, -'OTU ID') %>%
  spread('OTU ID', valname)

#arrange alphabetically so that vectors align
al1 <- plyr::arrange(Meta_A_B, SampleID, desc((SampleID)))
al2 <- plyr::arrange(trans_otu_A_B, SampleID, desc((SampleID)))

al3 <- as.factor(al1$Resequenced)

al477<- al2 %>% column_to_rownames(var="SampleID")
#add psuedocount
al478<- al477 + 1
al479<- al477 +.1


CV.BAL.dic <- selbal.cv(x = al477, y = al3, n.fold = 5, n.iter = 10,
                        covar = NULL, logit.acc = "AUC")



#method to join taxonomy with otu table on the NCBI ID column
otu_data <- left_join(last.ditch, taxontable2,
                      by = c("OTU ID" = "Feature ID"))



