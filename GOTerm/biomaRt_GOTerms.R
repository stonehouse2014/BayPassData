# clear environment
rm(list = ls())
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt", version = "3.8")
library("biomaRt")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)

library(biomaRt)
#ensembl <- useEnsembl(biomart = "genes", dataset = "pmajor_gene_ensembl")
#ensembl <- useMart(dataset = "pmajor_gene_ensembl", biomart = "ensembl")
#ensembl <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
#ensembl <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")



#ensembl
#ensembl=useMart("ensembl",dataset="ggallus_gene_ensembl")
#Turkey
ensembl=useMart("ensembl",dataset="mgallopavo_gene_ensembl")
#Collared Flycatcher
#ensembl=useMart("ensembl",dataset="falbicollis_gene_ensembl")


filters=listFilters(ensembl)
filters[1:5,]
attributes=listAttributes(ensembl)
attributes[1:5,]

#annotate a set of gene symbols with GO annotations related to biological processes
#genesym=c("ACSL4",
#"ADCY8",
#"ADGRB1",
#"ADGRB3",
#"AKAP6",
#"ALOX5",
#"ANKRD40",
#"ARHGAP12",
#"ARHGAP15",
#"ASCC2",
#"BAZ1A",
#"BCL9",
#"BLOC1S2",
#"BMPR1A",
#"BRSK2")

geneinfo = read.table("/Users/Jo/Desktop/transfer/PC1_genes.txt", sep = " ", row.names=NULL)
#genesym=read data from column 3 of file
genesym<-geneinfo[,11]
goids=getBM(attributes=c('hgnc_symbol','go_id'),
            filters='hgnc_symbol',
            values=genesym,
            mart=ensembl)
goids

#outfile<-("/Users/Jo/Desktop/transfer/GOTerms.txt ")
#fileConn<-file(outfile, open="w")

#write.table(goids, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, fileConn)
#close.connection(fileConn)
# output 
library(dplyr)
#remove empty go_ids
report<-filter(goids, !go_id =="") %>%
  group_by(hgnc_symbol) %>%         # sort by gene symbol
  summarize(go_id=do.call(paste, c(as.list(go_id), sep=", "))) #concatenate GO id
#write.table(report,file=("/Users/Jo/Desktop/transfer/GOTerms_desc_GT.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(report,file=("/Users/Jo/Desktop/transfer/GOTerms_desc_chicken.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(report,file=("/Users/Jo/Desktop/transfer/GOTerms_desc_zebrafinch.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(report,file=("/Users/Jo/Desktop/transfer/GOTerms_desc_human.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(report,file=("/Users/Jo/Desktop/transfer/GOTerms_desc_mouse.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(report,file=("/Users/Jo/Desktop/transfer/GOTerms_desc_Turkey.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

