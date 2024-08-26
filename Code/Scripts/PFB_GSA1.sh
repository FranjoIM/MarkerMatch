#!/bin/bash
#SBATCH --job-name=GCA_GSA1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fivankovic@broadinstitute.org
#SBATCH --ntasks=1
#SBATCH --mem=7gb
#SBATCH --qos=carolmathews-b
#SBATCH --time=4-00:00:00
#SBATCH --output=/Logs/LOG_%j-%x.log

date; hostname; pwd

module load penncnv

WKD="..."
PCN="/apps/penncnv/1.0.5"
GSA1SNPs="$WKD/SupportingFiles/GSA1_SNPPOS.txt"
GSA1samples="$WKD/SupportingFiles/GSA1_Samples.List"

cat ${GSA1samples} | while read item || [[ -n $line ]]
do
    $PCN/genomic_wave.pl \
        --prefix $WKD/Adjusted/TS/ \
        -adjust \
        -gcmodel $WKD/SupportingFiles/GSA1.gc \
        $item
done
