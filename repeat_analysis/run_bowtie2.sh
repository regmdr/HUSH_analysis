#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=3:00:0
#$ -l h_vmem=15G
#$ -N bowtie
#$ -cwd
#$ -j y

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "=========================================================="


module load samtools
module load bowtie2

samp=$1
ref=$2


read_pth=/data/Blizard-RoweLab/HUSH_KD/reads/trimmed_reads

bowtie2 --no-unal -x alignment_to_rep-consensus/${ref} -1 ${read_pth}/${samp}_R1_001_val_1.fq.gz -2 ${read_pth}/${samp}_R2_001_val_2.fq.gz | samtools view -bS - > ${samp}.${ref}.bam


samtools sort -m 5G -o ${samp}.${ref}.sorted.bam -T ${samp} ${samp}.${ref}.bam
samtools index ${samp}.${ref}.sorted.bam

rm ${samp}.${ref}.bam


echo $samp ${ref}

echo "*******************"
echo "Finished on : $(date)"
echo "*******************"
