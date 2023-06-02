#!/usr/bin/env Rscript
rm(list=ls())

library(getopt)
library(optparse)
require(data.table)
require(dplyr)
require(TissueEnrich)
require(tidyr)
library(cowplot))


cat(" ################################# \n # \n #  Read in data \n # \n ################################# \n")

inputGenes<- fread("/home1/ALC_sigGenes.csv",data.table=F)
inputGenes <- inputGenes[,1]
gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
output<-teEnrichment(inputGenes = gs, rnaSeqDataset = 1)
seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)


enrichmentOutputplot=enrichmentOutput[order(enrichmentOutput$`fold.change`,decreasing=T),]
coldata=data.frame(enrichmentOutputplot$Tissue,cols=c("#8FC1D9","#8FC18F","#F2BC78","#F09694","#C7AFD3","#FAD037",colorRampPalette(brewer.pal(12, "Paired"))(12)[seq(2,12,2)],colorRampPalette(brewer.pal(12, "Set3"))(12),colorRampPalette(brewer.pal(9, "Pastel1"))(9),colorRampPalette(brewer.pal(5, "Set2"))(2)))
coldata=coldata[order(enrichmentOutputplot$Tissue),]

p11 <- ggplot(enrichmentOutputplot,aes(x=reorder(Tissue,-fold.change),y=fold.change,label = Tissue.Specific.Genes,fill = Tissue))+
      geom_bar(stat = 'identity')+ labs(x='', y = 'Fold change')+ scale_y_continuous(breaks=seq(0,max(enrichmentOutput$`fold.change`),2))+
      theme_bw()+theme(legend.position="none")+
      theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,color = "black",size="12"),
        axis.text.y = element_text(color = "black",size="12"),
        axis.title = element_text(face="bold",size="14"),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major.y = element_line(color = "#EFEEED"),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+scale_fill_manual(values =  as.vector(coldata$cols))
p11
ggsave("/home1/Figure2C_ALC_HPA_Fold.pdf", width = 800, height = 500, scale = 1/80, dpi = 600)
fwrite(enrichmentOutputplot,"/home1/Figure2C_ALC_HPA_Fold.txt", col.names = T, row.names = F, quote = F, sep = "\t")

