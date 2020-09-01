#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=10:00:0
#$ -l h_vmem=7G
#$ -N htseqCount
#$ -cwd
#$ -j y

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "=========================================================="


module load python
module load samtools

samp=$1
ref_file=~/lab/reference_files/Homo_sapiens/UCSC/transcriptome_data/gencode.v30.annotation.gff3

samtools view ${samp}.sorted.bam | htseq-count -s reverse --type exon --idattr gene_id --mode intersection-nonempty - ${ref_file} > ${samp}.htseq-count.genes.out

echo $samp
echo "*******************"
echo "Finished on : $(date)"
echo "*******************"
