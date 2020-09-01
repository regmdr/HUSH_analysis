#!/bin/bash


samp=$1
ref=$2

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand

samtools view -b -f 128 -F 16 ${samp}.${ref}.sorted.bam > ${ref}/${samp}.fwd1.bam
# 16 = reverse strand, 128 = second in pair -> alignments of the second in pair if they don't map to the reverse strand
samtools index ${ref}/${samp}.fwd1.bam

samtools view -b -f 80 ${samp}.${ref}.sorted.bam > ${ref}/${samp}.fwd2.bam
# 80 = reverse, first in pair.
samtools index ${ref}/${samp}.fwd2.bam

# Combine alignments that originate on the forward strand.
#
samtools merge -f ${ref}/${samp}.fwd.bam ${ref}/${samp}.fwd1.bam ${ref}/${samp}.fwd2.bam
samtools index ${ref}/${samp}.fwd.bam

rm  ${ref}/${samp}.fwd1.bam ${ref}/${samp}.fwd2.bam 
rm  ${ref}/${samp}.fwd1.bam.bai ${ref}/${samp}.fwd2.bam.bai

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 ${samp}.${ref}.sorted.bam > ${ref}/${samp}.rev1.bam
# 144 = reverse strand, second in pair, 
samtools index ${ref}/${samp}.rev1.bam

samtools view -b -f 64 -F 16 ${samp}.${ref}.sorted.bam > ${ref}/${samp}.rev2.bam
# 64 = first in pair, 16 = reverse strand
samtools index ${ref}/${samp}.rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${ref}/${samp}.rev.bam ${ref}/${samp}.rev1.bam ${ref}/${samp}.rev2.bam
samtools index ${ref}/${samp}.rev.bam

rm ${ref}/${samp}.rev1.bam ${ref}/${samp}.rev2.bam
rm ${ref}/${samp}.rev1.bam.bai ${ref}/${samp}.rev2.bam.bai


echo $samp $ref

echo "*******************"
