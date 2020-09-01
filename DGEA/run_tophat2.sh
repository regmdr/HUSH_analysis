#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=20:00:0
#$ -l mem=35G
#$ -pe smp 4
#$ -N tophat2
#$ -cwd
#$ -j y

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "=========================================================="

module load gcc-libs/4.9.2
module load python/2.7.9
module load bowtie/1.1.2
module load bowtie2/2.2.5
module load tophat

samp=$1


tophat2 -o $samp -g 100  -p 4 \
                --transcriptome-index=/home/rekgmdr/Scratch/reference_files/Homo_sapiens/UCSC/transcriptome_data/hg38_transcriptome/gencode.v30.annotation \
        -G /home/rekgmdr/Scratch/reference_files/Homo_sapiens/UCSC/transcriptome_data/gencode.v30.annotation.gff3 \
        /home/rekgmdr/Scratch/reference_files/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
        ../reads/trimmed_reads/${samp}_R1_001_val_1.fq.gz ../reads/trimmed_reads/${samp}_R2_001_val_2.fq.gz


echo $samp
echo "*******************"
echo `qstat -j $JOB_ID | grep usage`
echo "*******************"
