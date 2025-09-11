#!/bin/bash
#SBATCH -J genomicsDB
#SBATCH -A snic2021-5-540
#SBATCH -p core -n 16
#SBATCH -t 1-02:00:00
#SBATCH -M snowy

ulimit -c unlimited

module load bioinfo-tools
module load bwa
module load samtools
module load GATK/4.1.4.1


#IN=/crex/proj/snic2020-6-185/Jay/proj5000/fixheader ## without the extenstion (.bam)
#cont=$1
#outdir=/crex/proj/snic2020-6-185/Jay/proj5000/GenomicDB/

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

inter=$1
#echo $inter

out=`echo $inter | sed 's/\_cont.*//'`
echo $out


cd /crex/proj/snic2020-6-185/Jay/proj5000/fixheader/$out



ref=/proj/snic2020-6-185/Jay/Yuxuan/ref/picea_abies.master-rna-scaff.nov2012_sorted.fa
#j=`ls *.gz | sed 's/^/\-V /g' | tr '\n' ' '`
ls *.gz >input.list

#echo $j
echo

echo

echo "Start program HaplotyCaller"

gatk --java-options "-Xmx110G" CombineGVCFs -R $ref \
	--variant input.list \
	-G StandardAnnotation \
	-O $out.combined.g.vcf.gz

#gatk --java-options "-Xmx90G" GenomicsDBImport \
#    --genomicsdb-workspace-path $out \
#    $j \
#    --merge-input-intervals true \
#    -L $inter  \
    

gatk --java-options "-Xmx90G" GenotypeGVCFs \
    -R $ref \
    -V $out.combined.g.vcf.gz \
    -O $out.vcf.gz \
    -L ../$inter \
    -new-qual


cd -

#rsync -rta $SNIC_TMP/$out $outdir
 
