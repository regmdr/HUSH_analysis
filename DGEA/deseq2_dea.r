# NOTE : The pairwise comparison results are run with apeglm

library("DESeq2")
library("scales")
library("RColorBrewer")
library("gplots")
library(ggplot2)
library("calibrate")
# We assign the name of the conditions to variables. This way we can use the same script for other analyses by doing some minor changes.
options(stringsAsFactors = FALSE)

s1<-"HUSH"

# Reading the metadata file:
meta <- read.delim(file = paste("meta_", s1,".txt", sep=""), header = F)
sampleTable<-data.frame(sampleName=meta$V1, fileName=meta$V2, condition=meta$V3, experiment=meta$V4)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "../gene_counts/", design= ~ condition)
colData(ddsHTSeq)$condition <- relevel(colData(ddsHTSeq)$condition, "shC")

# To runDESeq, we just need one line:
ddsHTSeq <- DESeq(ddsHTSeq)


#############################################################################################
### Differential Gene expression analysis


ctrl<-"shC"
cnd1<-"shM"

dirname<-paste0(cnd1, "_", ctrl)

coef<-paste0("condition_",cnd1,"_vs_",ctrl)
res <- lfcShrink(ddsHTSeq, coef=coef, type="apeglm")
dir.create(file.path(getwd (), dirname), showWarnings = FALSE)

####    Getting gene names
library("biomaRt")
mart = useMart("ENSEMBL_MART_ENSEMBL", host ="www.ensembl.org")
mart=useMart(biomart = "ENSEMBL_MART_ENSEMBL", host ="www.ensembl.org", dataset="hsapiens_gene_ensembl") 
results <- getBM (attributes=c("ensembl_gene_id_version", "hgnc_symbol"), filters="ensembl_gene_id_version", values=rownames(res), mart=mart)
idx <- match(rownames(res), results$ensembl_gene_id_version)
res$gene_name <- results[idx, "hgnc_symbol"]


# We filter out rowsthat have an "NA" in the adjusted p-value column
de.filtered = res[!is.na(res$padj), ]

### Tables
# We can already print the table with the results
write.table (res, file = paste (cnd1, "-", ctrl, "_all_genes.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

# We can subset out results to include only significantly, differentially expressed genes.
significant.genes=de.filtered[(de.filtered$padj < .05),]

# Se can also use a threshold for a minimum fold change:
up<- significant.genes[(significant.genes$log2FoldChange > 0),]
down<- significant.genes[(significant.genes$log2FoldChange < 0),]

# We write the filtered results to a plain text file. These will be "tab" delimited tables:
write.table(up,file=paste("up_", cnd1, "-", ctrl, ".genes.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")
write.table(down,file=paste("down_", cnd1, "-", ctrl, ".genes.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")


### checking by log2FC
up2<- significant.genes[(significant.genes$log2FoldChange > 1),]
down2<- significant.genes[(significant.genes$log2FoldChange < (-1)),]

# We write the filtered results to a plain text file. These will be "tab" delimited tables:
write.table(up2,file=paste("log2FC1/up_", cnd1, "-", ctrl, ".genes.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")
write.table(down2,file=paste("log2FC1/down_", cnd1, "-", ctrl, ".genes.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")



### MA Plot 

# MAplot
pdf(paste0(dirname, "/MAplot.pdf", sep=""))
        t<-quantile(de.filtered$log2FoldChange[de.filtered$log2FoldChange>0],0.999)
        plot(de.filtered$baseMean[abs(de.filtered$log2FoldChange)<t[[1]]], de.filtered$log2FoldChange[abs(de.filtered$log2FoldChange)<t[[1]]], col = ifelse((de.filtered$padj[abs(de.filtered$log2FoldChange)<t[[1]]]<=0.05 & abs(de.filtered$log2FoldChange[abs(de.filtered$log2FoldChange)<t[[1]]])>=1), alpha("violetred1", 0.5), alpha("gray50", 0.5)), xlab = "mean of normalized counts", ylab=expression(log[2]~fold~change), log = "x", cex=0.45, pch=16,  ylim=c((-t[[1]]-1),(t[[1]]+1)),  main=paste("Comparison of samples ", cnd1, "vs", ctrl, sep=" "),)
        points(de.filtered$baseMean[de.filtered$log2FoldChange>=t[[1]]], rep(t[[1]]+0.5, length(de.filtered$baseMean[de.filtered$log2FoldChange>=t[[1]]])), pch=17, col = ifelse(de.filtered$padj[de.filtered$log2FoldChange>=t[[1]]]<=0.05, alpha("violetred1", 0.5), alpha("gray50", 0.5)), cex=0.7)
        points(de.filtered$baseMean[de.filtered$log2FoldChange<=(-t[[1]])], rep((-t[[1]]-0.5), length(de.filtered$baseMean[de.filtered$log2FoldChange<=(-t[[1]])])), pch=17, col = ifelse(de.filtered$padj[de.filtered$log2FoldChange>=t[[1]]]<=0.05, alpha("violetred1", 0.5), alpha("gray50", 0.5)), cex=0.7)
        abline(a=1, b=0, col="darkviolet")
        abline(a=-1, b=0, col="darkviolet")
        legend("topright", "Genes with a padj <0.05", x.intersp=0.3, pch=20, col=alpha("violetred1", 0.5), cex=0.6, horiz=TRUE)
dev.off()
