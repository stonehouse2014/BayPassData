# clear environment
rm(list = ls())
library(tidyverse)
library(GGally)
setwd("/Volumes/GoogleDrive/My\ Drive/archive/redoBaypass/all/geneScore")
#setwd("G:/My Drive/JonSlate/PhDStudents/JoanneStonehouse/Chapters_Papers_etc/SubmittedVersion_EvolutionLetters")

#BayPassOutput <- read.table("PC1_genes_header.txt",sep="",head=T)
#BayPassOutput <- read.table("PC2_genes_header.txt",sep="",head=T)
#BayPassOutput <- read.table("PC3_genes_header.txt",sep="",head=T)
BayPassOutput <- read.table("PC4_genes_header.txt",sep="",head=T)


#Need to create header files for PC2, PC3 and PC4
head(BayPassOutput)
summary(BayPassOutput)

######################################
### Bit below in case you want    ####
### to look at the distributions  ####
### and correlations of the various ##
### metrics produced by BayPass  #####
######################################


#ggpairs(BayPassOutput[,3:10]) +
#  theme_bw()
#ggsave("BayPassOutputs_paired.png",width=8,height=8)

###################################################
##### Part 1 - read in the BayPass output, ########
##### group the results by gene, and the    #######
##### calculate the mean, max and sum of    #######
##### the statistic of choice. We also need   #####
##### to count up the number of SNPs per      #####
##### gene Here, I've used BF, but XtX is     #####
##### also justifiable.                     #######
##### Some genes have loads of SNPs, but    #######
##### that means there will few genes with  #######
##### that number of SNPs. We could group them  ###
##### e.g. all genes with 100-150 SNPs are ########
##### recoded so they have exactly 100 SNPs #######
##### Part 2 - group by the number of SNPs ########
##### per gene, and then within each group ########
##### rank the genes by max BF, in descending #####
##### order. i.e. highest BF has lowest rank ######
##### Could also rank by mean BF, but I think #####
##### max is easier to justify as it picks ########
##### up outliers, while controlling for ##########
##### the number of SNPs. Given low LD between ####
##### our SNPs, I think averaging over them #######
##### is not very sensible.  ######################

BayPassSums <- BayPassOutput %>%
  group_by (.,gene) %>%
  summarise(mean_BF = mean(BF),
            sum_BF = sum(BF),
            max_BF = max(BF),
            n_snps = n()) %>%
  mutate(n_snps = ifelse(n_snps >=100 & n_snps <150, 100, n_snps)) %>%
  mutate(n_snps = ifelse(n_snps >=150 & n_snps <200, 150, n_snps)) %>%
  mutate(n_snps = ifelse(n_snps >=200 & n_snps <300, 200, n_snps)) %>%
  mutate(n_snps = ifelse(n_snps >=300 & n_snps <500, 300, n_snps)) %>%
  mutate(n_snps = ifelse(n_snps >=500, 500, n_snps)) %>%
  group_by (.,n_snps) %>%
  mutate(rank_BF = rank(-max_BF, ties.method = "first"))


###############################################
#####  Part 3 - count up how many        ######
#####  genes have each number of SNPs      ####
#####  i.e. how many have 1 SNP, how many  ####
#####  have 2 SNPs, etc #######################
###############################################

SNPsPerGeneTally <- BayPassSums %>%
  count(.,n_snps)

###############################################
#### part 4 - combine the ranks and the #######
#### number of genes with each n_snps. ########
#### Then calculate the gene score by #########
#### dividing the rank by the number ##########
#### of genes with that many SNPs. This #######
#### gene score is like a P value, with #######
#### lower values stronger evidence for the ###
#### gene being associated with the PC ########
###############################################



BayPassSums2 <- dplyr::left_join(BayPassSums, SNPsPerGeneTally, by = "n_snps")
BayPassSums2$GeneScore <- BayPassSums2$rank_BF / BayPassSums2$n

#write.table(BayPassSums2, "BayPassSums2_PC1.txt",sep="\t",col.names=T,row.names=F,quote=F)
#write.table(BayPassSums2, "BayPassSums2_PC2.txt",sep="\t",col.names=T,row.names=F,quote=F)
#write.table(BayPassSums2, "BayPassSums2_PC3.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(BayPassSums2, "BayPassSums2_PC4.txt",sep="\t",col.names=T,row.names=F,quote=F)

