# HUSH_analysis

Repository for the code used in the analysis performed in:

	The HUSH complex is a gatekeeper of type I interferon induction through epigenetic control of LINE-1s



The main analysis described in this repo:

### Figure 2.
Differential gene expression analysis in **DGEA** folder. 

1. Alignment of reads was performed with TopHat2 for each sample _run\_tophat2.sh_. 
2. Counting of reads by gene with HTSeqCount _run\_htseq-count.sh_
3. DEA with _deseq2\_dea.r_. Gene set enrichment analysis was done with _fgsea\_dea.r_


### Figure 4. 
Differential repeat analysis in **repeat_analysis** folder. 
Counting by family using TEtranscripts with _run\_TEcount.sh_. 

##### Locus-specific differential expression analysis

* Counting of uniquely mapping reads on individual loci with _run\_htseq-count.rep-locus.sh_. 
* Differential expression analysis performed with DESeq2 for locus with at least 20 reads among all samples.
* Plots in 4d done with script _by\_fam\_log2FC.bound.py_. 
* Analysis of distance between TSS of upregulated ISGs and upregulated locus in 4f done with script _distance\_from\_loci.randomization\_median.py_ 

##### Analysis to consensus sequences. 

1. Alignment of reads to repeat consensus done with _run\_bowtie2.sh_. 
2. Read pairs were separated by strand of origin with script _split\_reads.sh_.
3. The depth at each base in the reference sequences was computed with _samtools depth_ and plotted in python. 
4. Total number of mapped reads was calculated with _samtools flagstat_. Values were normalized to account for difference in library size.

##### Assessment of bidirectional transcription. 

The steps for such analysis are described in **bidirectional\_transcription\_analysis.sh**



### Figure 5
Randmizations in figure 5d carried out with script _randomization\_l1-all.R_


### Figure 6

For 6a, TCGA data retrival and plotting was performed with _xenaPython\_GTex-TCGA\_controls.separate.py_

Plotting and statistical test in fig 6b were performed with _mpp8\_by\_immune\_subtype.py_


