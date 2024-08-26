#!/bin/bash
#SBATCH --job-name=PFB_OEE
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
OEESNPs="$WKD/SupportingFiles/OEE_SNPPOS.txt"

$PCN/compile_pfb.pl \
    --output $WKD/SupportingFiles/OEE_FullSet.pfb \
    --snpposfile $OEESNPs \
    --listfile $WKD/SupportingFiles/OEE_1000_SamplesList.txt
