#!/bin/bash
#SBATCH -J GenotypeGvcfs
#SBATCH -A snic2017-7-328
#SBATCH -p core -n 2
#SBATCH -t 0-03:00:00
#SBATCH --mail-user pascal.milesi@ebc.uu.se
#SBATCH --mail-type=ALL
#SBATCH -M snowy

ulimit -c unlimited


cont=$1
DB=$2
outdir=$3

mkdir $outdir

echo "Sync files"
rsync -ta ./gatk-4.0.10.0 $SNIC_TMP
rsync -ta $DB $SNIC_TMP
rsync -ta $cont $SNIC_TMP
rsync -ta ./sub* $SNIC_TMP

echo

cd $SNIC_TMP

ls

echo
ref=`ls sub*.fa`
echo $ref
echo
inter=`ls *.list`
echo $inter

out=`echo $inter | sed 's/\_list.*//'`
echo $out

out2=`echo $inter | sed 's/\_cont.*//'`
echo
echo $out2
echo
echo "Start program GenptypeGVFs"

gatk --java-options "-Xmx12G" GenotypeGVCFs \
    -R $ref \
    -V gendb://$out2 \
    -O $out.vcf.gz \
    -L $inter \
    -new-qual

cd -

rsync -rta $SNIC_TMP/$out*vcf* $outdir
