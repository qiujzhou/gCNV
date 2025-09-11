#!/bin/bash
#SBATCH -J merge
#SBATCH -A g2021010
#SBATCH -p core -n 8
#SBATCH -t 0-06:00:00
#SBATCH -M snowy


ulimit -c unlimited

module load bioinfo-tools
#module load bwa/0.7.17
module load samtools/1.8
#module load picard/2.23.4

batch=$1

cd /crex/proj/snic2020-6-185/Jay/proj5000/mark_dup


mkdir ../temp.$batch
sed -n "$(($batch*100-100+1)),$(($batch*100))p" ../list.4795.txt >../temp.$batch/batch.$batch.list

for te in $(cat ../temp.$batch/batch.$batch.list); do cp $te.*bam  ../temp.$batch/ ;done


cd ../temp.$batch 

echo "Merge"
samtools merge \
--threads 8 \
merged_$batch.samples.bam \
*.bam > samtools_merging.log 2>&1

samtools depth -aa  merged_$batch.samples.bam | gzip -c - > samtools_depth_$batch.samples.tsv.gz
mv merged_all.samples.bam ../temp

cd ../
#rm -rf ./temp.$batch
