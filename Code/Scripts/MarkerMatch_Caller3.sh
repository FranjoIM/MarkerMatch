#!/bin/bash

# Pull environment variables
Var_Sets=$1
Var_Sams=$2

# Define directories and files, make directories
WKD="..."
PCN=/apps/penncnv/1.0.5
HMM=$PCN/lib/hhall.hmm
ADJ=$WKD/Adjusted/ASD
PFB=$WKD/SupportingFiles/MarkerMatch/$Var_Sets.pfb
SAM=$WKD/SupportingFiles/OMNI_Samples.List$Var_Sams
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
SLG=$WKD/Logs/Validation/PennCNV_Caller3/$Var_Sets

mkdir -p $WKD/Validation/SSC/Raw $WKD/Validation/SSC/Log $WKD/Validation/SSC/Mer $WKD/Validation/SSC/Extra $SLG

sbatch<<EOT
#!/bin/bash
#SBATCH --job-name=PennCNV-Call_Set-${Var_Sets}_${Var_Sams}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=fivankovic@broadinstitute.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=carolmathews-b
#SBATCH --mem=7gb
#SBATCH --time=4-00:00:00
#SBATCH --output=${SLG}/SampleSet_${Var_Sams}

date; hostname; pwd

module load penncnv

# CALL CNVS
time($PCN/detect_cnv.pl \
    -test \
    -hmm $HMM \
    -pfb $PFB \
    --listfile $SAM \
    --gcmodelfile $GCM \
    -log ${LOG}_${Var_Sams} \
    --confidence \
    --output ${RAW}_${Var_Sams})
EOT
