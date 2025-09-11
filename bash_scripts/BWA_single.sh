#!/bin/bash
#SBATCH -J BWA
#SBATCH -A snic2020-5-508
#SBATCH -p core -n 6
#SBATCH -t 0-06:00:00



ulimit -c unlimited

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.8
module load picard/2.23.4

ref=/proj/snic2020-6-185/Jay/Yuxuan/ref/ #<reference.fa>
SAMPLE=$1
grp_head=/crex/proj/snic2020-6-184/private/mtiret/input/fastq/grep_header.txt
outdir=/crex/proj/snic2020-6-185/Jay/proj5000/cp_dataset/reads_from_sra/bam


#cd /crex/proj/snic2020-6-185/Jay/rawReads/combined
cd /crex/proj/snic2020-6-185/Jay/proj5000/cp_dataset/reads_from_sra
#gunzip *$SAMPLE*

#cat *$SAMPLE*R1* >$SAMPLE.combined.R1.fastq
#cat *$SAMPLE*R2* >$SAMPLE.combined.R2.fastq 

#gzip *$SAMPLE*


#mkdir $outdir/bam
#mkdir $outdir/sam
#mkdir $outdir/clean

reads1=$SAMPLE\_1.fastq.gz
reads2=`echo $reads1 | sed s/_1/_2/g`;

echo "Sync files"
rsync -ta $reads1 $SNIC_TMP
rsync -ta $reads2 $SNIC_TMP
rsync -ta $ref/picea* $SNIC_TMP
rsync -ta $grp_head $SNIC_TMP
#rsync -ta ./picard.jar $SNIC_TMP



echo "Start program"

cd $SNIC_TMP
echo
ls $SNIC_TMP
echo
#R1=`ls *R1*fastq.gz`
echo $R1
#R2=`ls *R2*fastq.gz`
echo $R2
genome=`ls *.fa`
echo $genome
echo
echo
sample_id=$SAMPLE


echo $sample_id
echo
header=`grep "$sample_id" $grp_head `
echo $header
echo
echo "Start BWA"
echo
#bwa mem -t 6 -M -R "$header" $genome $reads1 $reads2 > $sample_id.sam

bwa mem -t 6 -M $genome $reads1 $reads2 > $sample_id.sam

ls $SNIC_TMP

echo "Start samtools"

samfile=`ls *.sam`

samtools view -bS -q 10 $samfile > $sample_id.q10.bam
ls $SNIC_TMP

echo "Start program FixMate"
ls $SNIC_TMP

java -Xmx40g  -jar $PICARD_ROOT/picard.jar FixMateInformation \
    INPUT=$sample_id.q10.bam \
    OUTPUT=$sample_id.FM.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=STRICT \
    ADD_MATE_CIGAR=True \
    ASSUME_SORTED=TRUE

echo "Start program CleanSam"

java -Xmx40g  -jar $PICARD_ROOT/picard.jar CleanSam \
INPUT=$sample_id.FM.bam \
OUTPUT=$sample_id.clean.bam \
VALIDATION_STRINGENCY=STRICT

cd -

#rsync $SNIC_TMP/*.FM.bam $outdir/FM
rsync $SNIC_TMP/*clean.bam $outdir/clean/

echo "Done"
