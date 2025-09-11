#!/bin/bash
#SBATCH -J depth
#SBATCH -A snic2021-5-540
#SBATCH -p core -n 1
#SBATCH -t 0-06:00:00
#SBATCH -M snowy


ulimit -c unlimited

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.8
module load picard/2.23.4

cd /proj/snic2020-6-184/private/PM_LL/Jay/backup_proj5000/mark_dups

batch=$1

samtools depth -aa -b ../Pabies.probe60bp.bed $batch.dedup.bam | gzip -c - > $batch.samtools_depth_probe60bp.tsv.gz 
