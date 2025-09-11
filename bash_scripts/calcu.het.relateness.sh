#!/bin/bash
#SBATCH -J cal_het_relatedness
#SBATCH -A snic2020-5-508
#SBATCH -p core -n 1
#SBATCH -t 2-20:00:00


module load bioinfo-tools
module load vcftools/0.1.16
#module load samtools
#module load GATK/4.1.4.1

cd /proj/snic2020-6-185/Jay/proj5000/vcf

#batch=$1
vcf=proj5000.passfilter.biallelic.mind0.3.geno0.3.vcf
#vcf_prefix=`echo $batch | sed "s/\_SNPs.list//g"`

#gunzip raw_SNPs.hardfilter.vcf.gz
#vcf=proj5000.passfilter.biallelic2
#vcf_prefix=proj5000.passfilter.biallelic2
#echo $vcf_prefix
#zcat raw_SNPs.vcf.gz | grep -v "#" | awk -F "\t" '{print $1"\t"$2}' > raw.snp.position
#vcftools --vcf proj5000.passfilter.withfixedHeader.vcf --min-alleles 2 --max-alleles 2 --recode --out proj5000.passfilter.biallelic2
#grep -v "#" | awk -F "\t" '{print $1"\t"$2}' proj5000.passfilter.biallelic*vcf > passfilter.snp.position

echo "calculate hets"

vcftools --gzvcf proj5000.passfilter.biallelic.mind0.3.geno0.3.maf0.01.fixSNPname.pruned.random0.25.vcf.gz --het --out proj5000.passfilter.biallelic.mind0.3.geno0.3.maf0.01.random0.25

echo "calculate relatedness"

vcftools --gzvcf proj5000.passfilter.biallelic.mind0.3.geno0.3.maf0.01.fixSNPname.pruned.random0.25.vcf.gz --relatedness --out proj5000.passfilter.biallelic.mind0.3.geno0.3.maf0.01.ramdom0.25
#grep -wFf $batch ../$vcf > $vcf_prefix.SNPs

#echo "start merge Header"

#cat ../header.passfilter.biallelic.mind0.3.geno0.3.txt $vcf_prefix.SNPs > $vcf_prefix.vcf

gunzip proj5000.passfilter.withfixedHeader.vcf.gz
gatk --java-options "-Xmx35G" IndexFeatureFile -I proj5000.passfilter.withfixedHeader.vcf
gatk SelectVariants -restrict-alleles-to BIALLELIC -V proj5000.passfilter.withfixedHeader.vcf -O proj5000.passfilter.biallelic.vcf.gz
gatk SelectVariants -V proj5000.passfilter.biallelic.geno0.3.vcf --exclude-sample-name samples.mind0.3.list.args -O proj5000.passfilter.biallelic.mind0.3.geno0.3.vcf
gatk --java-options "-Xmx35G" VariantsToTable -V proj5000.passfilter.withfixedHeader.vcf -F QD -F MQ -F SOR -F QUAL -F MQRankSum -F ReadPosRankSum -F FS -O proj5000.passfilter.withfixedHeader.table

echo "done"
