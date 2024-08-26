#!/bin/bash
#SBATCH --job-name=PFB_OMNI
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fivankovic@broadinstitute.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7gb
#SBATCH --time=10-00:00:00
#SBATCH --output=/Logs/LOG_%j-%x.log

date; hostname; pwd

module load penncnv

WKD="..."
PCN="/apps/penncnv/1.0.5"
OMNISNPs="$WKD/SupportingFiles/ASD_SNPPOS.txt"
OUTPFB="$WKD/SupportingFiles/OMNI_FullSet.pfb"

$PCN/compile_pfb.pl \
    --output $OUTPFB \
    --snpposfile $OMNISNPs \
    --listfile $WKD/SupportingFiles/OMNI_1000_ParentsList.txt
