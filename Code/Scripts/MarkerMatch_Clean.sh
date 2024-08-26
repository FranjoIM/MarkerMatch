#!/bin/bash

# Pull environment variables
Var_Sets=$1

# Define directories and files, make directories
WKD="..."
PCN=/apps/penncnv/1.0.5
HMM=$PCN/lib/hhall.hmm
ADJ=$WKD/Adjusted/ASD
SNPPOS=$WKD/SupportingFiles/ASD_SNPPOS.txt
GCM=$WKD/SupportingFiles/OMNI.gc
F_TEL=$WKD/SupportingFiles/TelomereFilter.txt
F_CEN=$WKD/SupportingFiles/CentromereFilter.txt
F_IMU=$WKD/SupportingFiles/ImmuneFilter.txt
F_SED=$WKD/SupportingFiles/SegDupFilter.txt
RAW=$WKD/Validation/SSC/Raw/$Var_Sets.Raw
LOG=$WKD/Validation/SSC/Log/$Var_Sets.Log
MER=$WKD/Validation/SSC/Mer/$Var_Sets.Mer
O_TEL=$WKD/Validation/SSC/Extra/$Var_Sets.Telomere
O_CEN=$WKD/Validation/SSC/Extra/$Var_Sets.Centromere
O_IMU=$WKD/Validation/SSC/Extra/$Var_Sets.Imumune
O_SED=$WKD/Validation/SSC/Extra/$Var_Sets.SegDup
SLG=$WKD/Logs/Validation/PennCNV_Clean/

mkdir -p $WKD/Validation/SSC/Raw $WKD/Validation/SSC/Log $WKD/Validation/SSC/Mer $WKD/Validation/SSC/Extra $SLG

sbatch<<EOT
#!/bin/bash
#SBATCH --job-name=PennCNV-Clean_Set-${Var_Sets}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=fivankovic@broadinstitute.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=carolmathews
#SBATCH --mem=7gb
#SBATCH --time=4-00:00:00
#SBATCH --output=${SLG}/$Var_Sets

date; hostname; pwd

module load penncnv

# MERGE CNVS
time($PCN/clean_cnv.pl \
    --output $MER \
    --signalfile $SNPPOS \
    --bp \
    combineseg \
    $RAW)

# ANNOTATE TELOMERES
time($PCN/scan_region.pl \
    --verbose \
    --minqueryfrac 0.5 \
    --append \
    $MER \
    $F_TEL > \
    $O_TEL)

# ANNOTATE CENTROMERES
time($PCN/scan_region.pl \
    --verbose \
    --minqueryfrac 0.5 \
    --append \
    $MER \
    $F_CEN > \
    $O_CEN)

# ANNOTATE IMMUNE
time($PCN/scan_region.pl \
    --verbose \
    --minqueryfrac 0.5 \
    --append\
    $MER \
    $F_IMU > \
    $O_IMU)

# ANNOTATE SEGDUPS
time($PCN/scan_region.pl \
    --verbose \
    --minqueryfrac 0.5 \
    --append \
    $MER \
    $F_SED > \
    $O_SED)

EOT
