#!/bin/bash
#SBATCH -J BWA
#SBATCH -A g2021010
#SBATCH -p core -n 3
#SBATCH -t 0-18:00:00
#SBATCH -M snowy

ulimit -c unlimited

module load bioinfo-tools
module load GATK/4.1.4.1
module load bwa
module load samtools
module load BEDTools/2.29.2


#batch=$1
#DB=$2
#outdir=$3

#mkdir $outdir

#echo "Sync files"
#rsync -ta ./gatk-4.0.10.0 $SNIC_TMP
#rsync -ta $DB $SNIC_TMP
#rsync -ta $cont $SNIC_TMP
#rsync -ta ./sub* $SNIC_TMP

echo

cd /crex/proj/snic2020-6-185/Jay/proj5000/temp

#zcat samtools_depth_$batch.samples.tsv.gz | cut -f 3 > $batch.depth

# ls
gunzip samtools_depth_1.samples.tsv.gz

paste samtools_depth_1.samples.tsv *depth > samtools_depth_all.samples.tsv


#gunzip samtools_depth_$batch.samples.tsv.gz
grep  -v 0$'\t'0$'\t'0$'\t'0$'\t'0$'\t' samtools_depth_all.samples.tsv | awk -v val=24000 'BEGIN{OFS="\t"; FS="\t"} {if($3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30+$31+$32+$33+$34+$35+$36+$37+$38+$39+$40+$41+$42+$43+$44+$45+$46+$47+$48 >= val){print $1, ($2 -1), $2}}' | bedtools merge -i - > depth__filter_all.5n.bed
#awk -v val=4800 'BEGIN{OFS="\t"; FS="\t"} {if($3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14 >= val){print $1, ($2 - 1), $2}}' samtools_depth_all.samples.rm0.tsv > depth_summary_filter_1n.bed


#rsync -rta $SNIC_TMP/$out*vcf* $outdir

