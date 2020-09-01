"""
Availability of the following software is assumed:
	bedtools v2.27.1
	samtools 1.9
	bedGraphToBigWig (script dowloaded from: http://hgdownload.soe.ucsc.edu/admin/exe/)

The python script used for this analysis is included in the repo.
"""

### 1) Merge replicates using samtools:
# Input should be sorted and filtered for multimapping reads.
samtools merge ${samp} ${samp1} ${samp2} ${samp3}

### 2) Using merged unique alignments divide reads by strand of origin:
	# Forward strand.
	# 1. alignments of the second in pair if they map to the forward strand
	# 2. alignments of the first in pair if they map to the reverse  strand

	samtools view -b -f 128 -F 16 ${samp}.bam > ${samp}.fwd1.bam
	# 16 = reverse strand, 128 = second in pair -> alignments of the second in pair if they don't map to the reverse strand
	samtools index ${samp}.fwd1.bam

	samtools view -b -f 80 ${samp}.bam > ${samp}.fwd2.bam
	# 80 = reverse, first in pair.
	samtools index ${samp}.fwd2.bam

	# Combine alignments that originate on the forward strand.
	samtools merge -f ${samp}.fwd.bam ${samp}.fwd1.bam ${samp}.fwd2.bam
	samtools index ${samp}.fwd.bam


	# Reverse strand
	# 1. alignments of the second in pair if they map to the reverse strand
	# 2. alignments of the first in pair if they map to the forward strand

	samtools view -b -f 144 ${samp}.bam > ${samp}.rev1.bam
	# 144 = reverse strand, second in pair, 
	samtools index ${samp}.rev1.bam

	samtools view -b -f 64 -F 16 ${samp}.bam > ${samp}.rev2.bam
	# 64 = first in pair, 16 = reverse strand
	samtools index ${samp}.rev2.bam

	# Combine alignments that originate on the reverse strand.
	samtools merge -f ${samp}.rev.bam ${samp}.rev1.bam ${samp}.rev2.bam
	samtools index ${samp}.rev.bam

### 3) Calculate genome-wide coverage using bedtools. In order to make it comparable, calculate a scaling factor for the samples.
"""
Got number of aligned pairs from align_summary result from TopHat2 and used it to calculate scaling factors
# Mapped reads in shC_b2.bam 142187570 = 1
# Mapped reads in shM_b2.bam 137846343 = 1.031493233
"""
	LC_COLLATE=C
	# for reverse pairs:
	bedtools genomecov -bg -split -scale $scale -ibam ${samp}.rev.bam -g ../hg38.genome > ${samp}.rev.bg
	sort -k1,1 -k2,2n ${samp}.rev.bg > ${samp}.rev.sorted.bdg
	bedGraphToBigWig ${samp}.rev.sorted.bdg ../hg38.genome ${samp}.rev.bw

	# for forward pairs:
	bedtools genomecov -bg -split -scale $scale -ibam ${samp}.fwd.bam -g ../hg38.genome > ${samp}.fwd.bg
	sort -k1,1 -k2,2n ${samp}.fwd.bg > ${samp}.fwd.sorted.bdg
	bedGraphToBigWig ${samp}.fwd.sorted.bdg ../hg38.genome ${samp}.fwd.bw


### 4) Run script to assess if a window is expressed from both strands:
	python bin_expression_by_strand.py ${samp}.fwd.bw ${samp}.rev.bw shM_b2

### 5) Get the subset of L1s in the human genome that overlap bidirectional locus
bedtools intersect -u -a hg38.L1s.chrFilt.gtf  -b shM_b2_expressed_both_strands85_non-zero.500.values.txt > bidirectional_L1s.shM_b2.85_500.gtf

### 6) Which repeat families are most represented?
cat bidirectional_L1s.shM_b2.85_500.gtf | awk '{split($10, VNAME, "_"); print VNAME[1]}' | sed 's/\"//' | sort | uniq -c | sort -nr
 # 199 L1HS
 #  74 L1PA2
 #  42 L1PA3
 #  26 L1MC4