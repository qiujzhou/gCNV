#!/bin/bash
#SBATCH -J fix_header_gentree
#SBATCH -A snic2021-5-540
#SBATCH -p core -n 2
#SBATCH -t 0-01:00:00
#SBATCH --mail-user qiujie.zhou@ebc.uu.se
#SBATCH --mail-type=FAIL


ulimit -c unlimited

module load bioinfo-tools
module load tabix


IN=$1 ## without the extenstion (g.vcf.gz)
outdir=/crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/fixheader

#mkdir $outdir

cd /crex/proj/snic2020-6-184/private/PM_LL/Jay/analysis_add_russia/gvcf

echo "Sync files"
#rsync -ta ./gatk-4.0.10.0 $SNIC_TMP
rsync -ta $IN.exact_target.g.vcf.gz $SNIC_TMP
#rsync -ta $ref/sub* $SNIC_TMP
rsync -ta exclude.contig.list $SNIC_TMP

echo

cd $SNIC_TMP

#ls picea*
ls

ind=`ls *vcf.gz`

echo
echo $ind
echo
sample_id=`echo $ind | sed 's/\..*//g'`
echo $sample_id
echo
outfile=$sample_id\_reduced.g.vcf.gz
echo $outfile
echo
echo "Start program fixheader"

zgrep -vwFf exclude.contig.list $ind | bgzip -c > $outfile && tabix -p vcf $outfile

cd -

rsync $SNIC_TMP/*reduced* $outdir
