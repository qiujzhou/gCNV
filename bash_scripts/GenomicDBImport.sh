#!/bin/bash
#SBATCH -J genomicsDB
#SBATCH -A snic2021-5-540
#SBATCH -p core -n 16
#SBATCH -t 0-08:00:00



ulimit -c unlimited

module load bioinfo-tools
module load bwa
module load samtools
module load GATK/4.1.4.1


IN=/crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/fixheader ## without the extenstion (.bam)
inter=$1
outdir=/crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/GenomicDB/

#INbai=`echo $(basename $IN) | sed 's/\.bam/\.bai/'`

#mkdir $outdir

#echo "Sync files"
#rsync -ta ./gatk-4.0.10.0 $SNIC_TMP
#rsync -ta $IN/*.gz* $SNIC_TMP
#rsync -ta $cont $SNIC_TMP


echo

#cd $SNIC_TMP

#ls picea*


#ls


#echo $inter

out=`echo $inter | sed 's/\_cont.*//'`
echo $out


cd $IN



ref=/proj/snic2020-6-185/Jay/Yuxuan/ref/picea_abies.master-rna-scaff.nov2012_sorted.fa
j=`ls *.gz | sed 's/^/\-V /g' | tr '\n' ' '`

#echo $j
echo

echo

echo "Start program HaplotyCaller"

gatk --java-options "-Xmx120G" GenomicsDBImport \
    --genomicsdb-workspace-path $out \
    $j \
    --merge-input-intervals true \
    -L $inter  \
    

#gatk --java-options "-Xmx120G" GenotypeGVCFs \
#    -R $ref \
#    -V gendb://$out \
#    -O $out.vcf.gz \
#    -L $inter \
#    -new-qual
#

cd -

#rsync -rta $SNIC_TMP/$out $outdir

