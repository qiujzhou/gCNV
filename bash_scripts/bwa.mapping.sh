#!/bin/bash
#SBATCH -J BWA
#SBATCH -A uppmax2021-2-19
#SBATCH -p core -n 10
#SBATCH -t 0-20:00:00
#SBATCH -M snowy

ulimit -c unlimited

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.8
module load picard/2.23.4

cd /crex/proj/snic2020-6-184/private/PM_LL/Jay/drones

genome=/crex/proj/snic2020-6-185/Jay/drones/ncbi-genomes-2021-10-18/GCF_003254395.2_Amel_HAv3.1_genomic.fna

sample_id=$1

#bwa index -a bwtsw -b 5000000 $genome
bwa mem -t 10 -M $genome $sample_id\_1.fastq.gz $sample_id\_2.fastq.gz > $sample_id.sam
samtools view -bS -q 10 $sample_id.sam > $sample_id.q10.bam
rm $sample_id.sam

ls $SNIC_TMP

java -Xmx40g  -jar $PICARD_ROOT/picard.jar FixMateInformation \
    INPUT=$sample_id.q10.bam \
    OUTPUT=$sample_id.FM.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=STRICT \
    ADD_MATE_CIGAR=True \
    ASSUME_SORTED=TRUE
rm $sample_id.q10.bam

echo "Start program CleanSam"

java -Xmx40g  -jar $PICARD_ROOT/picard.jar CleanSam \
INPUT=$sample_id.FM.bam \
OUTPUT=$sample_id.clean.bam \
VALIDATION_STRINGENCY=STRICT
rm $sample_id.FM.bam


echo "done"
