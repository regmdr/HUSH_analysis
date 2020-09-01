#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=3:00:0
#$ -l mem=30G
#$ -N TEcount
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


TEcount --sortByPos --stranded reverse --format BAM --mode multi --GTF ~/reference_files/Homo_sapiens/UCSC/transcriptome_data/gencode.v30.annotation.gff3 --TE ~/hg38_rmsk_TE.gtf -b alignment/${samp}/accepted_hits.bam --project ${samp}


echo $samp
echo "*******************"
echo `qstat -j $JOB_ID | grep usage`
echo "*******************"
