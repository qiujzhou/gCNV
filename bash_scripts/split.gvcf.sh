#!/bin/bash
#SBATCH -J genomicsDB
#SBATCH -A snic2021-5-540
#SBATCH -p core -n 16
#SBATCH -t 1-20:00:00
#SBATCH -M snowy

ulimit -c unlimited

module load bioinfo-tools
#module load bwa
#module load samtools
module load GATK/4.1.4.1


batch=$1
batch_sample=$2

cd /crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/batch.list/


inter=`ls $batch | sed "s/\_contig.list//g" `

mkdir $inter
cp $batch $inter
cd $inter

grep -vwFf $batch ../../keep.contig.list > excluded.contigs.list

for te in $(cat /crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/fixheader/$batch_sample);
  do zgrep -vwFf excluded.contigs.list ../../fixheader/complete_gvcf/$te\_reduced.g.vcf.gz | bgzip -c > $te.$inter.vcf.gz && tabix -p vcf $te.$inter.vcf.gz;
  done
  
  
#cd /crex/proj/snic2020-6-185/Jay/proj5000/fixheader/$inter


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
	-O $inter.combined.g.vcf.gz

#gatk --java-options "-Xmx90G" GenomicsDBImport \
#    --genomicsdb-workspace-path $out \
#    $j \
#    --merge-input-intervals true \
#    -L $inter  \
    

gatk --java-options "-Xmx110G" GenotypeGVCFs \
    -R $ref \
    -V $inter.combined.g.vcf.gz \
    -O $inter.vcf.gz \
    -L ../$batch \
    -new-qual


cd -

#rsync -rta $SNIC_TMP/$out $outdir
 
