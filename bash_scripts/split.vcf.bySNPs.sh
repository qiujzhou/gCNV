#!/bin/bash
#SBATCH -J SNP_split
#SBATCH -A g2021010
#SBATCH -M snowy
#SBATCH -p core -n 1
#SBATCH -t 0-20:00:00


module load bioinfo-tools
module load vcftools/0.1.16
#module load samtools
module load GATK/4.1.4.1

cd /proj/snic2020-6-185/Jay/proj5000/vcf/splitted

batch=$1
vcf=proj5000.passfilter.biallelic.mind0.3.geno0.3.vcf
vcf_prefix=`echo $batch | sed "s/\_SNPs.list//g"`

#gunzip raw_SNPs.hardfilter.vcf.gz
#vcf=proj5000.passfilter.biallelic2
#vcf_prefix=proj5000.passfilter.biallelic2
#echo $vcf_prefix
#zcat raw_SNPs.vcf.gz | grep -v "#" | awk -F "\t" '{print $1"\t"$2}' > raw.snp.position
#vcftools --vcf proj5000.passfilter.withfixedHeader.vcf --min-alleles 2 --max-alleles 2 --recode --out proj5000.passfilter.biallelic2
#grep -v "#" | awk -F "\t" '{print $1"\t"$2}' proj5000.passfilter.biallelic*vcf > passfilter.snp.position

echo "start grep"

grep -wFf $batch ../$vcf > $vcf_prefix.SNPs

echo "start merge Header"

cat ../header.passfilter.biallelic.mind0.3.geno0.3.txt $vcf_prefix.SNPs > $vcf_prefix.vcf

#gunzip proj5000.passfilter.withfixedHeader.vcf.gz
#gatk --java-options "-Xmx35G" IndexFeatureFile -I proj5000.passfilter.withfixedHeader.vcf
#gatk SelectVariants -restrict-alleles-to BIALLELIC -V proj5000.passfilter.withfixedHeader.vcf -O proj5000.passfilter.biallelic.vcf.gz
#gatk SelectVariants -V proj5000.passfilter.biallelic.geno0.3.vcf --exclude-sample-name samples.mind0.3.list.args -O proj5000.passfilter.biallelic.mind0.3.geno0.3.vcf
#gatk --java-options "-Xmx35G" VariantsToTable -V proj5000.passfilter.withfixedHeader.vcf -F QD -F MQ -F SOR -F QUAL -F MQRankSum -F ReadPosRankSum -F FS -O proj5000.passfilter.withfixedHeader.table

echo "done"
