# HUSH_analysis

Repository for the code used in the analysis performed in:

	The HUSH complex is a gatekeeper of type I interferon through epigenetic regulation of LINE-1s



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


##### Species specificity data.

Species-specificity of repeat families was extracted from Dfam (using the "TaxId" field, see dfam.speciesID.txt).The species relationships were obtained using:

	ete3 ncbiquery --search 1437010 314295 314293 9392 9347 32523 33554 9526 32524 314147 9989 207598 376913 9443 117571 9255 311790 9348 9605 314146 32525 9604 9606 91561 40674 9263 --tree > tree_human-repeats.TaxonID.txt


Relative proportion of families specific to each taxonomic level were calculated with _human-repeats\_by\_cladeID.py_

Depiction of such proportions into the tree were performed in python 3 using _ete3\_tree.v2.py_



##### Analysis to consensus sequences. 

1. Alignment of reads to repeat consensus done with _run\_bowtie2.sh_. 
2. Read pairs were separated by strand of origin with script _split\_reads.sh_.
3. The depth at each base in the reference sequences was computed with _samtools depth_ and plotted in python. 
4. Total number of mapped reads was calculated with _samtools flagstat_. Values were normalized to account for difference in library size.

##### Assessment of bidirectional transcription. 

The steps for such analysis are described in **bidirectional\_transcription\_analysis.sh**



### Figure 5
Randomizations in figure 5d carried out with script _randomization\_l1-all.R_


### Figure 6

For 6a, TCGA data retrieval and plotting was performed with _xenaPython\_GTex-TCGA\_controls.separate.py_

Plotting and statistical test in fig 6b were performed with _mpp8\_by\_immune\_subtype.py_


