#!/bin/bash
#SBATCH -J HaplotypeCaller
#SBATCH -A snic2021-5-540
#SBATCH -p core -n 12
#SBATCH -t 0-08:00:00
#SBATCH -M snowy

ulimit -c unlimited

module load bioinfo-tools
module load GATK/4.1.4.1
#module load bwa
#module load samtools
#module load R

DATA=$1 #list of ".sam" file
ref=/proj/snic2020-6-185/Jay/Yuxuan/ref/
target=/crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/depth__filter_all.5n.bed
outdir=/crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/gvcf

#mkdir $outdir

cd /crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/mark_dup

echo "Sync files"
rsync -ta $DATA.fix.bam $SNIC_TMP
rsync -ta $DATA.*.bai $SNIC_TMP
rsync -ta $target $SNIC_TMP
rsync -ta $ref/picea* $SNIC_TMP
echo

cd $SNIC_TMP

ref2=`ls *.fa`
#target2=`ls *.bed`
echo "Start program"


    echo
    echo 'HaplotyCaller'
    echo
    gatk --java-options "-Xmx100G" HaplotypeCaller -R $ref2 -I $DATA.fix.bam -L $target --native-pair-hmm-threads 12 -ERC GVCF -O $DATA.exact_target.g.vcf.gz
    #gatk --java-options "-Xmx64G" HaplotypeCaller -R $ref2 -I $DATA.dedup.bam --emit-ref-confidence GVCF -O $DATA.exact_target.g.vcf.gz

	echo
    echo

cd -

rsync $SNIC_TMP/*.vcf.gz* $outdir

