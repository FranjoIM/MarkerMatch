#!/bin/bash
#SBATCH --job-name=GCA_OEE
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
OEESNPs="$WKD/SupportingFiles/OEE_SNPPOS.txt"
OEEsamples="$WKD/SupportingFiles/OEE_Samples.List"

cat ${OEEsamples} | while read item || [[ -n $line ]]
do
    $PCN/genomic_wave.pl \
        --prefix $WKD/Adjusted/OEE/ \
        -adjust \
        -gcmodel $WKD/SupportingFiles/OEE.gc \
        $item
done
