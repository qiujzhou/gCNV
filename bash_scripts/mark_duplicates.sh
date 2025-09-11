#!/bin/bash
#SBATCH -J Mark_duplicates
#SBATCH -A snic2020-5-508
#SBATCH -p core -n 6
#SBATCH -t 0-2:00:00

ulimit -c unlimited

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.8
module load picard/2.23.4

DATA=$1  #SAMPLE NAME
outdir=/proj/snic2020-6-185/Jay/drones  #output directory, created before running



echo "Sync files"
rsync -ta /proj/snic2020-6-185/Jay/drones/$DATA.fix.bam $SNIC_TMP
#rsync -ta ./picard.jar $SNIC_TMP

echo

cd $SNIC_TMP
ls $SNIC_TMP
ind=$DATA.fix.bam

echo "Start program"

sample_id=$DATA
    echo $sample_id
    echo
    echo 'Markduplicates'
    echo

java -Xmx40g  -jar $PICARD_ROOT/picard.jar MarkDuplicatesWithMateCigar \
    INPUT=$ind\
    OUTPUT=$sample_id.dedup.bam \
    REMOVE_DUPLICATES=True \
    CREATE_INDEX=True \
    VALIDATION_STRINGENCY=LENIENT \
    METRICS_FILE=$sample_id.duplicates_metrics.txt \
    MINIMUM_DISTANCE=400 \
    BLOCK_SIZE=10000000
    echo
    echo "Done"
    echo
cd -

rsync $SNIC_TMP/*dedup.bam $outdir
rsync $SNIC_TMP/*dedup.bai $outdir
rsync $SNIC_TMP/*.txt $outdir
