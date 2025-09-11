#!/bin/bash
#SBATCH -J gatk_select
#SBATCH -A snic2021-5-540
#SBATCH -p core -n 16
#SBATCH -t 1-06:00:00
#SBATCH -M snowy


ulimit -c unlimited

module load bioinfo-tools
#module load bwa
#module load samtools
module load GATK/4.1.4.1
module load bcftools/1.12

#batch=$1
#batch_sample=$2

cd /crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/raw/

bcftools concat --file-list vcf.list --output-type z --output variants.vcf.gz

gatk --java-options "-Xmx30G" IndexFeatureFile -I variants.vcf.gz
  
gatk --java-options "-Xmx30G" SelectVariants \
-R /proj/snic2020-6-185/Jay/Yuxuan/ref/picea_abies.master-rna-scaff.nov2012_sorted.fa \
-V variants.vcf.gz \
--select-type-to-include SNP \
-O raw_SNPs.vcf.gz

gatk --java-options "-Xmx30G" VariantFiltration \
        -V raw_SNPs.vcf.gz \
        --verbosity ERROR \
        --filter-name "QDfilter" \
        --filter-expression "QD < 2.0 || MQ < 40.0 || SOR > 3.0 || QUAL < 20.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" \
        -O raw_SNPs.hardfilter.vcf
gatk --java-options "-Xmx100G" SelectVariants \
        -V raw_SNPs.hardfilter.vcf \
        --verbosity ERROR \
        --exclude-filtered true \
        --exclude-non-variants true \
        --restrict-alleles-to BIALLELIC \
        -O cline.passfilter.vcf

#gatk SelectVariants -V proj5000.passfilter.biallelic.vcf.gz --exclude-sample-name samples.mind0.3.list.args -L mind0.15.geno0.3.maf0.01.SNPs.list.bed -O proj5000.passfilter.biallelic.mind0.3.geno0.3.maf0.01.vcf