#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=10:00:0
#$ -l mem=7G
#$ -l tmpfs=10G
#$ -N htseqCount
#$ -cwd
#$ -j y

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "=========================================================="



module purge
module load gcc-libs/4.9.2
module load python2/recommended
module load samtools


samp=$1

ref_file=~/by_instance_counts_human/hg38.by_instance.no_chrUn.gtf 

samtools sort -m 5G -n -o ${samp}.sorted.bam -T ${samp} ${alignment_path}/${samp}/accepted_hits.bam
samtools view ${samp}.sorted.bam | htseq-count --order name -s no --type repeat --idattr repeat_ID --mode intersection-nonempty - ${ref_file} > ${samp}.htseq-count.repeats.out


echo $samp
echo "*******************"
echo "Finished on : $(date)"
echo "*******************"
