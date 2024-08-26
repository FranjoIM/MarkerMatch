#!/bin/bash
#SBATCH --job-name=GCA_OMNI
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fivankovic@broadinstitute.org
#SBATCH --ntasks=1
#SBATCH --mem=7gb
#SBATCH --time=7-00:00:00
#SBATCH --output=/Logs/LOG_%j-%x.log

date; hostname; pwd

module load penncnv

WKD="..."
PCN="/apps/penncnv/1.0.5"
OMNISNPs="$WKD/SupportingFiles/ASD_SNPPOS.txt"
OMNIsamples="$WKD/SupportingFiles/OMNI_Samples.List"

cat ${OMNIsamples} | while read item || [[ -n $line ]]
do
    $PCN/genomic_wave.pl \
        --prefix $WKD/Adjusted/ASD/ \
        -adjust \
        -gcmodel $WKD/SupportingFiles/OMNI.gc \
        $item
done
