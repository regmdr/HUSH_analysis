# Following instructions from Reimand 2019:
# Note, I also ran GSEa using the Java applet,
# The resulta are here, and comparable with what I get in fgsea
# For hallmark gene sets:
# //	/Users/rocio.enriquez-gasca/gsea_home/output/may29/my_analysis.GseaPreranked.1559147865379
# For reactome gene sets:
# //	/Users/rocio.enriquez-gasca/gsea_home/output/may29/my_analysis.GseaPreranked.1559148351685 



library(tidyverse)
library("ggpubr")
library("DOSE")
library("biomaRt")
library(ReactomePA)
library("pathview")
library("clusterProfiler")


comparison<-"shM-shC"

dea<-read.table(paste0(comparison, "_all_genes.txt"), sep="\t", header=T, row.names=1)

# Filter unnanotatted genes 
include <- rownames(dea)[dea$gene_name!="" & !is.na(dea$gene_name)]
annotated <- dea[include,]
annotated.filtered = annotated[!is.na(annotated$padj), ]


#calculate ranks
ranks_RNAseq = data.frame(GeneNAme=annotated.filtered$gene_name, rank=sign(annotated.filtered$log2FoldChange) * -log10(annotated.filtered$pvalue))
ranks_RNAseq <- ranks_RNAseq[order(as.numeric(ranks_RNAseq[,2]),decreasing = TRUE),]

ranks_RNAseq <- ranks_RNAseq[is.finite(ranks_RNAseq$rank),]
# had to remove Inf

# Using the Hallmark gene set from MSigDB. Hallmark gene sets summarize and represent specific well-defined biological states or processes and display coherent expression. 
# I have downloaded the Hallmark gene sets from here : http://software.broadinstitute.org/gsea/downloads.jsp
library(fgsea)
library(ggrepel)

pathways.hallmark <- gmtPathways("~/Documents/projects/MPP8_KD/gene_analysis_dose/h.all.v6.2.symbols.gmt")
ranks<-deframe(ranks_RNAseq)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


topfsgea<-fgseaResTidy[1:10,]
topfsgea$lab<-gsub("HALLMARK_", "", topfsgea$pathway)
topfsgea$lab<-gsub("_", " ", topfsgea$lab)
topfsgea$lab<-tolower(topfsgea$lab)

titl <- gsub("-", "/", comparison)

p<- ggplot(topfsgea, aes(y=NES,x=pval, label=lab)) + geom_point() + labs(x="p-value", y="Normalized Enrichment Score") +
   geom_text_repel(fontface = "bold", colour="dodgerblue4") + ggtitle(titl) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 


pdf(paste0(comparison, ".fgsea.scatter.pdf"))
print(p)
dev.off()