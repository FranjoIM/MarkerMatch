#!/bin/bash
#SBATCH --job-name=PFB_GSA1
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
GSA1SNPs="$WKD/SupportingFiles/GSA1_SNPPOS.txt"

$PCN/compile_pfb.pl \
    --output $WKD/SupportingFiles/GSA1_FullSet.pfb \
    --snpposfile $GSA1SNPs \
    --listfile $WKD/SupportingFiles/GSA1_1000_ParentsList.txt
