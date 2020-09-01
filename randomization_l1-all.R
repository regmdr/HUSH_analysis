library("scales")
library("RColorBrewer")
library("gplots")
library("ggpubr")
library(VennDiagram)


# kd MPP8 dea
shMPP8_file<-"shM-shC_all_genes.txt"
shMPP8 <- read.table(shMPP8_file, header=TRUE, row.names=1, sep="\t")


ens_id<-function(x){
    unlist(strsplit(x, ".", fixed = TRUE))[1]
}

mpp8_genes_ID<-as.vector(sapply(rownames(shMPP8), ens_id))
rownames(shMPP8)<-mpp8_genes_ID

nTot<-nrow(shMPP8)

#####
significant.shMPP8<-shMPP8[(shMPP8$padj < .05),]
significant.shMPP8<-significant.shMPP8[!is.na(significant.shMPP8$baseMean),]

shMPP8_up<-significant.shMPP8[(significant.shMPP8$log2FoldChange > 1),]
shMPP8_down<-significant.shMPP8[(significant.shMPP8$log2FoldChange < (-1)),]

# Choosing  upregulated
gene_subset<-shMPP8_up
suffix<-"up"

nDE<-nrow(gene_subset)

# In order to fairly assess the signficance, I need to look at the subset of ensembl IDs that were used for the L1 analysus:
 # L1 response:
l1_file<-"L1_response/ensID_htseq_output_with_stat_summary.txt"
l1_table_all<-read.table(l1_file, header=F, sep="\t")
# Will subset the significantly upregulated genes with a padj <0.05 and start with log2FC > 2
l1_filtered = l1_table_all[!is.na(l1_table_all$V13), ]


# gene_universe:
gene_universe<-l1_filtered$V14


# L1-induced ISGs
l1_file<-"L1_response/gene-symbol_to_ensembl.FC10.manual.txt"
l1_table<-read.table(l1_file, header=F, sep="\t")
l1_genes<-l1_table$V2

intr<-length(intersect(rownames(gene_subset), l1_genes))

randomOvlp_vector<-c()
nRand<-10000
for(i in 1:nRand){
	temp<-sample(gene_universe, nDE)
	temp.intr<-length(intersect(temp, l1_genes))
	randomOvlp_vector<-c(randomOvlp_vector,temp.intr)
}

rand_mean<-mean(randomOvlp_vector)
nMoreX<-length(c(1:length(randomOvlp_vector))[randomOvlp_vector>=intr])
pval1<-(nMoreX+1)/(nRand+1)


df<-data.frame(number=c(intr, rand_mean), sd=c(0, sd(randomOvlp_vector)), type=c("shM up", "randomized genes"))

library(ggplot2)
p<- ggplot(df, aes(x=type, y=number)) +  geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=number, ymax=number+sd), width=.2, position=position_dodge(.9)) +labs(title="L1-induced genes",
     x="gene subset", y = "Number of overlapping genes") + theme( plot.title = element_text(hjust = 0.5)) +
    geom_text(aes(x=2, y=3), label = paste0("p-val <", scientific(pval1)))

pdf("shM-up_l1-all.barplot.pdf")
print(p)
dev.off()

phyper(intr, nDE, (length(gene_universe)-nDE), length(l1_genes), lower.tail=F) 
# [1] 5.1203e-24


df_rand<-data.frame(number=randomOvlp_vector[1:10000], type=rep("randomized genes",10000))
p<- ggplot(df, aes(x=type, y=number)) +  geom_bar(stat="identity", fill=c("lightgrey","dodgerblue"), position=position_dodge()) +labs(title="L1-induced genes",
     x="gene subset", y = "Number of overlapping genes") + theme( plot.title = element_text(hjust = 0.5)) +
    geom_point(data=df_rand,aes(x=type, y=number),position=position_jitter(width =.15), colour="darkgrey", fill="darkgrey", alpha=0.6)+ geom_errorbar(aes(ymin=number, ymax=number+sd), width=.2, position=position_dodge(.9)) +
    geom_text(aes(x=2, y=3), label = paste0("p-val <", scientific(pval1))) +  theme_classic()



pdf("shM-up_l1-all.barplot_with_dots.pdf")
print(p)
dev.off()
