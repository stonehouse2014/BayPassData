# topGO Example Using Kolmogorov-Smirnov (KS) Testing

rm(list=ls())

## Load the required R packages

library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
library(AnnotationDbi)

## Initialise variables

#Initialise principal component, PC1->PC4
PC<-"PC2"
#ontology either "BP", "MF", or "CC"
#ontology<-"BP"
#ontology<-"MF"
ontology<-"CC"
nodeSize<-10

################################
## Step 4: Performing the gene enrichment tests


## Step 1: Define gene universe

BAMfilename="/Users/Jo/Desktop/transfer/GOTerms_desc_tab.txt"
geneID2GO<-readMappings(file=BAMfilename, sep='\t', IDsep=",")
str(geneID2GO, list.len=5)
geneUniverse<-names(geneID2GO)


## Step 2: Define Gene Scores
#filenameScorePC2<-file(paste0("/Volumes/GoogleDrive/My\ Drive/archive/redoBaypass/all/geneScore/BayPassSums2_",PC,".txt"))
assign(paste0("filenameScore",PC),file(paste0("/Volumes/GoogleDrive/My\ Drive/archive/redoBaypass/all/geneScore/BayPassSums2_",PC,".txt")))


#geneScorePC2<-read.table(filenameScorePC2, header=TRUE)
assign(paste0("geneScore",PC),read.table(get(eval(paste0("filenameScore",PC))),header=TRUE))

#geneListScorePC2<-geneScorePC2$GeneScore
assign(paste0("geneListScore",PC),eval(parse(text=paste0("geneScore",PC,"$GeneScore"))))

#Not dynamic variable name (need to edit each time)
names(geneListScorePC2)<-eval(parse(text=paste0("geneScore",PC,"$gene")))





## Define gene score function
topDiffGenes<-function(allScore) {
  return(allScore <0.01)
}
#x<-topDiffGenes(get(eval(paste0("geneScore",PC))))
x<-topDiffGenes(get(eval(paste0("geneListScore",PC))))
sum(x)

## Step 3: Define GOdata object
#Create topGo object
#default scoreOrder="increasing" for p-values

assign(paste0("topGoGeneScore",PC),new("topGOdata", description=paste0("ParusMajor",PC), ontology, allGenes=get(eval(paste0("geneListScore",PC))), geneSel= topDiffGenes,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = nodeSize))

assign(paste0("SVsign",PC),sigGenes(get(eval(paste0("topGoGeneScore",PC)))))


## Step 4: Perform gene enrichment tests with KS statistic
### run weight01 algorithm with KS statistic
max.print=100000

assign(paste0("resultWeight01KS",PC),runTest(get(eval(paste0("topGoGeneScore",PC))), algorithm="weight01", statistic="KS"))
#assign(paste0("resultWeight01KS",PC),runTest(get(eval(paste0("topGoGeneScore",PC))), algorithm="weight01", statistic="Fisher"))

get(eval(paste0("resultWeight01KS",PC)))
#SVweight01<-GenTable(topGoGeneScorePC1,KS = get(eval(paste0("resultWeight01KS",PC))), topNodes=20)
#SVweight01

## Step 6: Result table of default algorithm "weight01", with KS statistic
#disable scientific notation i.e. 1e-03 to use KS<=0.01 threshold later on
options(scipen = 999)
allGO=usedGO(get(eval(paste0("topGoGeneScore",PC))))
SVweight01<-GenTable(get(eval(paste0("topGoGeneScore",PC))),KS = get(eval(paste0("resultWeight01KS",PC))), orderBy=KS, topNodes=length(allGO))
#Not dynamic variable naming
SVweight01.p<-SVweight01[which(SVweight01$KS<=0.01),]
SVweight01.p
sign_weight01_KS<-nrow(SVweight01.p)

###############

### Plot enrichment Weight01 KS
#https://www.biostars.org/p/350710/
goEnrichment<-GenTable(get(eval(paste0("topGoGeneScore",PC))),KS = get(eval(paste0("resultWeight01KS",PC))), topNodes=sign_weight01_KS)
goEnrichment <- goEnrichment[goEnrichment$KS<=0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
goEnrichment$KS <- as.numeric(goEnrichment$KS)

require(ggplot2)
ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab(ontology) +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=12) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()


### Plot Go hierarchy for Weight01 KS
# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
#diagram for top nodes for topgo default, algorithm.
# print a graph (to a pdf file) with the top 'numsignif' results:
pdf(file=paste0('/Users/Jo/Google\ Drive/archive/GOtermEnrichment/network/',PC,'_',ontology,'NetworkPlot.pdf'), paper='a4r', pointsize=14)
source("/Users/Jo/Google\ Drive/archive/GOtermEnrichment/GOplotFunction.R")
par(cex=0.6)
#par(mar=c(0,0,0,0))
showSigOfNodes(get(eval(paste0("topGoGeneScore",PC))), score(get(eval(paste0("resultWeight01KS",PC)))), firstSigNodes = sign_weight01_KS, sigForAll=FALSE, useInfo="all", .NO.CHAR=20)
#.NO.CHAR=18/30/16

printGraph(get(eval(paste0("topGoGeneScore",PC))),get(eval(paste0("resultWeight01KS",PC))), firstSigNodes = sign_weight01_KS, refResult =get(eval(paste0("resultWeight01KS",PC))),  useInfo="null",
           fn.prefix =paste0(
             '/Users/Jo/Google\ Drive/archive/GOtermEnrichment/network/',
             'topGO_topNodes',PC,'_',ontology), pdfSW=TRUE)
#pdfSW=TRUE
#firstSigNodes = 18
#firstSigNodes = sign_weight01_KS
#useInfo="'all'/'null'/'pval'/'counts'/'def'/'np',"

dev.off()


## Step 7: Correcting for multiple testing for Weight01 with KS statistic

#performing BH correction on our p values
p.adj<-round(p.adjust(SVweight01$KS, method="BH"), digits=4)

# create the file with all the statistics from GO analysis
#SVweight01_final=cbind(SVweight01,p.adj)
#SVweight01_final=SVweight01_final[order(SVweight01_final$p.adj),]
#SVweight01_final[1:(sign_weight01_KS+5),]

#get list of significant GO before multiple testing correction
results.table.p2<-SVweight01[which(SVweight01$KS<=0.01),]

#get list of significant GO after multiple testing correction
results.table.bh=SVweight01[which(SVweight01$p.adj<=0.05),]






## Step 8: Get all the genes in your significant GO TERMS for weight01 KS
#print out the genes that are annotated with the significantly enriched GO terms
#R code for this sub-section edited from step 6 of web tutorial: 
#https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/


mytermsSV<-SVweight01.p$GO.ID

mygenesSV<-genesInTerm(get(eval(paste0("topGoGeneScore",PC))), mytermsSV)
var=c()

for (i in 1:length(mytermsSV))
{
  if (length(mytermsSV)==0) {
    break
  }
  mytermSV <- mytermsSV[i]
  mygenesfortermSV <- mygenesSV[mytermSV][[1]]
  mygenesfortermSV <- paste(mygenesfortermSV, collapse=',')
  #line below to print to screen
  print(paste("GOTerm",mytermSV,"genes:",mygenesfortermSV))
  
  #var[i]<-(paste("GOTerm",mytermSV,"genes:",mygenesfortermSV))
  var[i]<-(paste(mytermSV,",",mygenesfortermSV))
}
write.table(var,paste0("/Users/Jo/Desktop/SVgenetoGOmapping",PC,".txt"),sep="\t",quote=F, col.names=FALSE)


# TopGO analysis for default weight01 algorithm with Fisher statistic
#Steps 1-3 as per the previous section (KS statistic) to set up the gene list
#Fisher statistic uses rank counts instead of the gene scores.
## Step 4: run gene enrichment test with default weight01 algorithm and Fisher statistic.
#max.print=100000
assign(paste0("resultWeight01FisherGS",PC),runTest(get(eval(paste0("topGoGeneScore",PC))), algorithm="weight01", statistic="Fisher"))
get(eval(paste0("resultWeight01FisherGS",PC)))


##Step 5: comparison weight01 SV & Fisher - All annotated SNPs
#allGO=(get(eval(paste0("topGoGeneScore",PC))))
weight01comp<-GenTable(get(eval(paste0("topGoGeneScore",PC))), weight01KS = get(eval(paste0("resultWeight01KS",PC))), weight01Fisher = get(eval(paste0("resultWeight01FisherGS",PC))), orderBy = "weight01KS", ranksOf="weight01Fisher", topNodes=sign_weight01_KS)
weight01comp 


## Step 6: Result default weight01 Fisher - All annotated SNPs
allGOSV=usedGO(get(eval(paste0("topGoGeneScore",PC))))
weight01FisherGS_table<-GenTable(get(eval(paste0("topGoGeneScore",PC))), weight01FisherGS=get(eval(paste0("resultWeight01FisherGS",PC))), orderBy='weight01FisherGS', topNodes=length(allGOSV))
#weight01FisherGS_table

result.table.fis<-weight01FisherGS_table[which(weight01FisherGS_table$weight01FisherGS<=0.01),]
result.table.fis
numsigfis<-nrow(result.table.fis)


### GO hierarchy Weight01 Fisher - All annotated SNPs
#pdf(file=paste0('/Users/Jo/Fisher/','topGO_topNodesPC',PC,'_',ontology,'.pdf'), height=12, width=12, paper='special', pointsize=18)
par(cex=0.5)
showSigOfNodes(get(eval(paste0("topGoGeneScore",PC))), score(get(eval(paste0("resultWeight01FisherGS",PC)))), useInfo="all", sigForAll=FALSE, firstSigNodes = numsigfis , .NO.CHAR=50)


## Step 7: Correcting for multiple testing Weight01 Fisher - All annotated SNPs
#performing BH correction on our p values
p.adjFis<-round(p.adjust(weight01FisherGS_table$weight01FisherGS, method="BH"), digits=4)

# create the file with all the statistics from GO analysis
weight01FisherGS_all=cbind(weight01FisherGS_table,p.adjFis)
weight01FisherGS_final=weight01FisherGS_all[order(weight01FisherGS_all$p.adjFis),]
weight01FisherGS_final[1:(numsigfis+5),]

## Step 8: Get all the genes in significant GO terms for Weight01 Fisher - All annotated SNPs
mytermsFis<-result.table.fis$GO.ID

mygenesFis<-genesInTerm(get(eval(paste0("topGoGeneScore",PC))), mytermsFis)
var=c()

for (i in 1:length(mytermsFis))
{
  if (length(mytermsFis)==0) {
    break
  }
  mytermFis <- mytermsFis[i]
  mygenesfortermFis <- mygenesFis[mytermFis][[1]]
  mygenesfortermFis <- paste(mygenesfortermFis, collapse=',')
  #line below to print to screen
  print(paste("GOTerm",mytermFis,"genes:",mygenesfortermFis))
  
  #var[i]<-(paste("GOTerm",mytermSV,"genes:",mygenesfortermSV))
  var[i]<-(paste(mytermFis,",",mygenesfortermFis))
}
write.table(var,paste0("/Users/Jo/Desktop/FisgenetoGOmappingPC",PC,".txt"),sep="\t",quote=F, col.names=FALSE)

