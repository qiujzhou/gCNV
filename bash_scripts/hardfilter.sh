#!/bin/bash
#SBATCH -A snic2020-5-508
#SBATCH -p core -n 5
#SBATCH -t 08:00:00
#SBATCH -J gatKfilter


ulimit -c unlimited
module load bioinfo-tools
module load gcc/10.1.0
module load GATK/4.0.8.0


#vcf=$data
#out=$out
#sp=$sp
#rsync -ta $data $SNIC_TMP
#cd $S
cd /crex/proj/snic2020-6-185/Jay/drones/vcf

gatk --java-options "-Xmx30G" VariantFiltration \
        -V drones.raw.vcf.gz \
        --verbosity ERROR \
        --filter-name "QDfilter" \
        --filter-expression "QD < 2.0 || MQ < 40.0 || SOR > 3.0 || QUAL < 20.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0" \
        -O raw_SNPs.hardfilter.vcf
gatk --java-options "-Xmx30G" SelectVariants \
        -V raw_SNPs.hardfilter.vcf \
        --verbosity ERROR \
        --exclude-filtered true \
        --exclude-non-variants true \
        --restrict-alleles-to BIALLELIC \
        -O drones.passfilter.vcf
#rsync -ta ${sp}_f2.vcf $out
