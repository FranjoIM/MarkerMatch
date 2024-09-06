### Environments
#### Local
Local analyses were ran on a Windows machine. Specifications:
|Component| Value|
|---|---|
|Processor|	12th Gen Intel(R) Core(TM) i9-12900KS   3.40 GHz|
|Installed RAM|	128 GB (128 GB usable)|
|System type|	64-bit operating system, x64-based processor|
|Edition|	Windows 11 Pro|
|Version|	23H2|
|OS build|	22631.4037|

#### Cloud
Cloud analyses were ran oh HiPerGator, a Linux cloud computing service, using SLURM scheduler to schedule/run jobs.

### Prepare Data

Ran locally. Using Illumina GenomeStudio, prepare data by clustering microarray genotyping samples _de novo_. Once done, exclude non-autosomes using the filter 
function, and export following files:
1. Illumina Final Report (Columns: `SNP Name`, `Sample ID`, `Chr`, `Position`, `B Allele Freq`, `Log R Ratio`),
2. PLINK Input Report (Leaving default options, except `UseForwardStrand` should be set to `True`),
3. Copy-numer metrics (CNM file, export `Name`, `Chr`, `Position` and `Log R Ratio` for all samples columns from the `Full Data Table` tab as tab-delimeted file),
4. Manifest file (MAN file, export `Name`, `Chr`, `Position` and `Minor Freq` columns from the `SNP Table` tab as a tab-delimeted file).
5. SNP Pos file (export `Name`, `Chr`, and `Position` columns from the `SNP Table` tab as a tab-delimeted file)

### Complete the Manifest Files

In R, run on a local machine, use thr [ManifSum.R](Scripts/ManifSum.R) function to suummarize the copy-number metrics (CNM) export into copu-number summary (CNS) file.

```R
# Set working directory
W_DIR <- "..."
setwd(W_DIR)

# Load packages
library(tidyverse)
library(matrixStats)

# Source function
source("./ManifSum.R")

# Run functions
ManifSum(fPath="DATA/OMNI_CNM.txt", oPath="DATA/OMNI_CNS.csv", nSam=4239)
ManifSum(fPath="DATA/GSA_CNM.txt", oPath="DATA/GSA_CNS.csv", nSam=4389)
ManifSum(fPath="DATA/OEE_CNM.txt", oPath="DATA/OEE_CNS.csv", nSam=3920)

# Import LRR summaries and manifests
OMNI_LRR <- read_csv("DATA/OMNI_CNS.csv")
OMNI_MAN <- read_delim("DATA/OMNI_MAN.txt", delim="\t")
GSA_LRR <- read_csv("DATA/GSA_CNS.csv")
GSA_MAN <- read_delim("DATA/GSA_MAN.txt", delim="\t")
OEE_LRR <- read_csv("DATA/OEE_CNS.csv")
OEE_MAN <- read_delim("DATA/OEE_MAN.txt", delim="\t")

# Remove irrelevant chromosomes and ensure consistend column names
OMNI_LRR <- OMNI_LRR %>%
  filter(!chr %in% c("0", "MT", "X", "Y", "XY")) %>%
  rename(Name=snpID, Chr=chr, Position=pos, LRR_mean=mean, LRR_sd=sd) %>%
  mutate(Chr=as.character(Chr))
OMNI_MAN <- OMNI_MAN %>%
  filter(!Chr %in% c("0", "MT", "X", "Y", "XY")) %>%
  rename(BAF=`Minor Freq`)
GSA_LRR <- GSA_LRR %>%
  filter(!chr %in% c("0", "MT", "X", "Y", "XY")) %>%
  rename(Name=snpID, Chr=chr, Position=pos, LRR_mean=mean, LRR_sd=sd) %>%
  mutate(Chr=as.character(Chr))
GSA_MAN <- GSA_MAN %>%
  filter(!Chr %in% c("0", "MT", "X", "Y", "XY"))  %>%
  rename(BAF=`Minor Freq`)
OEE_LRR <- OEE_LRR %>%
  filter(!Chr %in% c("0", "MT", "X", "Y", "XY")) %>%
  mutate(Chr=as.character(Chr))
OEE_MAN <- OEE_MAN %>%
  filter(!Chr %in% c("0", "MT", "X", "Y", "XY"))  %>%
  rename(BAF=`Minor Freq`)

# Populate manifest data with LRR data
OMNI_MAN <- OMNI_MAN %>%
  left_join(OMNI_LRR, by=c("Name", "Chr", "Position"))
GSA_MAN <- GSA_MAN %>%
  left_join(GSA_LRR, by=c("Name", "Chr", "Position"))
OEE_MAN <- OEE_MAN %>%
  left_join(OEE_LRR, by=c("Name", "Chr", "Position"))

# Export the complete manifests
OMNI_MAN %>% write_delim("DATA/OMNI_MAN.txt", delim="\t", col_names=T)
GSA_MAN %>% write_delim("DATA/GSA_MAN.txt", delim="\t", col_names=T)
OEE_MAN %>% write_delim("DATA/OEE_MAN.txt", delim="\t", col_names=T)
```

### Run MarkerMatch Algorithm
In R, run on a local machine, use thr [MarkerMatch.R](MarkerMatch.R) function to run the MarkerMatch algorithm. The function will print out time elapsed for each run.

```R
# Set working directory
W_DIR <- "..."
setwd(W_DIR)

# Load packages
library(tidyverse)

# Source function
source("./MarkerMatch.R")

# Import array manifests
OMNI_MAN <- read_delim("DATA/OMNI_MAN.txt", delim="\t")
GSA_MAN <- read_delim("DATA/GSA_MAN.txt", delim="\t")
OEE_MAN <- read_delim("DATA/OEE_MAN.txt", delim="\t")

# Calculate perfect overlaps
OMNI_MAN %>%
  filter(Name %in% GSA_MAN$Name) %>%
  select(Name) %>%
  write_delim("DATA/PerfectMatch_OMNI_MatSelect", col_names=FALSE)

OEE_MAN %>%
  filter(Name %in% GSA_MAN$Name) %>%
  select(Name) %>%
  write_delim("DATA/PerfectMatch_OEE_MatSelect", col_names=FALSE)

GSA_MAN %>%
  filter(Name %in% OEE_MAN$Name) %>%
  select(Name) %>%
  write_delim("DATA/PerfectMatch_GSA_MatSelect", col_names=FALSE)

# Run MarkerMatch
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10, Factor="Distance", OutPath="DATA/Pos_10")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10, Factor="BAF.delta", OutPath="DATA/BAF_10")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_10")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_10")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50, Factor="Distance", OutPath="DATA/Pos_50")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50, Factor="BAF.delta", OutPath="DATA/BAF_50")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_50")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_50")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100, Factor="Distance", OutPath="DATA/Pos_100")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100, Factor="BAF.delta", OutPath="DATA/BAF_100")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_100")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_100")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500, Factor="Distance", OutPath="DATA/Pos_500")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500, Factor="BAF.delta", OutPath="DATA/BAF_500")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_500")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_500")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000, Factor="Distance", OutPath="DATA/Pos_1000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000, Factor="BAF.delta", OutPath="DATA/BAF_1000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_1000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_1000")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000, Factor="Distance", OutPath="DATA/Pos_5000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000, Factor="BAF.delta", OutPath="DATA/BAF_5000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_5000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_5000")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10000, Factor="Distance", OutPath="DATA/Pos_10000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10000, Factor="BAF.delta", OutPath="DATA/BAF_10000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_10000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=10000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_10000")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50000, Factor="Distance", OutPath="DATA/Pos_50000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50000, Factor="BAF.delta", OutPath="DATA/BAF_50000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_50000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=50000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_50000")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100000, Factor="Distance", OutPath="DATA/Pos_100000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100000, Factor="BAF.delta", OutPath="DATA/BAF_100000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_100000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=100000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_100000")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500000, Factor="Distance", OutPath="DATA/Pos_500000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500000, Factor="BAF.delta", OutPath="DATA/BAF_500000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_500000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=500000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_500000")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000000, Factor="Distance", OutPath="DATA/Pos_1000000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000000, Factor="BAF.delta", OutPath="DATA/BAF_1000000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_1000000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=1000000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_1000000")

MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000000, Factor="Distance", OutPath="DATA/Pos_5000000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000000, Factor="BAF.delta", OutPath="DATA/BAF_5000000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000000, Factor="LRR_sd.delta", OutPath="DATA/LRRsd_5000000")
MarkerMatch(Reference=GSA_MAN, Matching=OMNI_MAN, MaxD=5000000, Factor="LRR_mean.delta", OutPath="DATA/LRRmean_5000000")

MarkerMatch(Reference=GSA_MAN, Matching=OEE_MAN, MaxD=10000, Factor="Distance", OutPath="DATA/OEE_Pos_10000")
MarkerMatch(Reference=GSA_MAN, Matching=OEE_MAN, MaxD=10000, Factor="BAF.delta", OutPath="DATA/OEE_BAF_10000")
MarkerMatch(Reference=GSA_MAN, Matching=OEE_MAN, MaxD=10000, Factor="LRR_sd.delta", OutPath="DATA/OEE_LRRsd_10000")
MarkerMatch(Reference=GSA_MAN, Matching=OEE_MAN, MaxD=10000, Factor="LRR_mean.delta", OutPath="DATA/OEE_LRRmean_10000")
```

### Graph Runtimes
In R, ran locally. Code for creating runtime graphs for MarkerMatch algorithms, manuscript **Figure 2**.
```R
# Set working directory
W_DIR <- "..."
setwd(W_DIR)

# Load packages
library(tidyverse)
library(ggplot)

# Import array manifests
OMNI_MAN <- read_delim("DATA/OMNI_MAN.txt", delim="\t")
GSA_MAN <- read_delim("DATA/GSA_MAN.txt", delim="\t")
OEE_MAN <- read_delim("DATA/OEE_MAN.txt", delim="\t")
```

### Graph Marker Density Plots
In R, ran locally. Code for creating marker density plots, manuscript **Figure 3**. First, the data cen be prepared using [UsedSNPs.R](Scripts/UsedSNPs.R) script, and then plotted using [Plot3.R](Scripts/Plot3.R).

### Graph Gaps Plots
In R, ran locally. Code for creating marker gaps plots, manuscript **Figure 4**. First, the data cen be plotted using [Plot4.R](Scripts/Plot4.R).

### Graph BAF Plots
In R, ran locally. Code for creating BAF plots, manuscript **Figure 5**. First, the data cen be plotted using [Plot5.R](Scripts/Plot5.R).

### Graph LRR-SD Plots
In R, ran locally. Code for creating LRR-sd plots, manuscript **Figure 6**. First, the data cen be plotted using [Plot6.R](Scripts/Plot6.R).

### Graph LRR-Mean Plots
In R, ran locally. Code for creating LRR-mean plots, manuscript **Figure 7**. First, the data cen be plotted using [Plot7.R](Scripts/Plot7.R).

### Split Final Reports
On the cloud, after Illumina final reports have been uploaded, split the final reports. `split_illumina_report.pl` is a part of PennCNV software.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# CREATE NEW DIRECTORIES
mkdir -p ${WKD}/Intensities/TS ${WKD}/Intensities/ASD ${WKD}/Intensities/OEE

# SPLIT FINAL REPORTS
$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/TS/ \
    $OLD/MarkerMatch/Experiment1/ChipData/TS_Wave1_FR.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/TS/ \
    $OLD/MarkerMatch/Experiment1/ChipData/TS_Wave2_FR.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/TS/ \
    $OLD/MarkerMatch/Experiment1/ChipData/TS_Wave3_FR.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/ASD/ \
    $OLD/MarkerMatch/SSC_Full_Autosomes.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/OEE/ \
    $NEW/Intensities/OEE/OEE_033s.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/OEE/ \
    $NEW/Intensities/OEE/OEE_328s1.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/OEE/ \
    $NEW/Intensities/OEE/OEE_328s2.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/OEE/ \
    $NEW/Intensities/OEE/OEE_328s3.txt

$PCN/split_illumina_report.pl \
    -p $WKD/Intensities/OEE/ \
    $NEW/Intensities/OEE/OEE_328s4.txt
```

### Make GC Model Files for Datasets
On the cloud, create GC model files for the datasets. Also upload SNP Pos files for this purpose. Scripts for [cal_gc_snp.pl](Scripts/cal_gc_snp.pl) and [gc_wrangle.r](Scripts/gc_wrangle.r) are provided.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# SPECIFY SNP POS FILE LOCATIONS
GSA1SNPs="$WKD/SupportingFiles/GSA1_SNPPOS.txt"
OMNISNPs="$WKD/SupportingFiles/ASD_SNPPOS.txt"
OEESNPs="$WKD/SupportingFiles/OEE_SNPPOS.txt"

# SORT GC MODEL FILE
zcat $PCN/gc_file/hg19.gc5Base.txt.gz |\
  sort -k2,2 -k3,3n \
  > $WKD/SupportingFiles/hg19.gc5Base_sorted.txt

# MAKE GC MODEL FILES FOR DATASETS
# cal_gc_snp.pl was adjusted slightly to resolve a bug that affects
# delimation of the GC percetage and column positioning
$WKD/Scripts/cal_gc_snp.pl \
    $WKD/SupportingFiles/hg19.gc5Base_sorted.txt \
    $GSA1SNPs > \
    $WKD/SupportingFiles/GSA1.hg19.gcmodel

$WKD/Scripts/cal_gc_snp.pl \
    $WKD/SupportingFiles/hg19.gc5Base_sorted.txt \
    $OMNISNPs > \
    $WKD/SupportingFiles/OMNI.hg19.gcmodel

$WKD/Scripts/cal_gc_snp.pl \
    $WKD/SupportingFiles/hg19.gc5Base_sorted.txt \
    $OEESNPs > \
    $WKD/SupportingFiles/OEE.hg19.gcmodel

# FIX GC MODEL FILE OUTPUT FOR THE PIPELINE
module load R
Rscript --vanilla $WKD/Scripts/gc_wrangle.r $GSA1SNPs $WKD/SupportingFiles/GSA1.hg19.gcmodel $WKD/SupportingFiles/GSA1.gc
Rscript --vanilla $WKD/Scripts/gc_wrangle.r $OMNISNPs $WKD/SupportingFiles/OMNI.hg19.gcmodel $WKD/SupportingFiles/OMNI.gc
Rscript --vanilla $WKD/Scripts/gc_wrangle.r $OEESNPs $WKD/SupportingFiles/OEE.hg19.gcmodel $WKD/SupportingFiles/OEE.gc
```

### Generate PFB Files for All Datasets
On the cloud, generate PFB files for datasets, used for CNV calling. Scripts for [PFB_GSA1.sh](Scripts/PFB_GSA1.sh), [PFB_OMNI.sh](Scripts/PFB_OMNI.sh), and [PFB_OEE.sh](Scripts/PFB_OEE.sh) are provided.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# GENERATE A LIST OF 1000 PARENTS FROM OMNI DATA
ls -d \
    "$WKD"/Intensities/ASD/* > \
    $WKD/SupportingFiles/OMNI_Samples.List

printf ".fa_\n.mo_" > \
    $WKD/SupportingFiles/Parent.Pattern

grep -Ff \
    $WKD/SupportingFiles/Parent.Pattern \
    $WKD/SupportingFiles/OMNI_Samples.List > \
    $WKD/SupportingFiles/OMNI_ParentsList.txt

shuf -n 1000 \
    $WKD/SupportingFiles/OMNI_ParentsList.txt > \
    $WKD/SupportingFiles/OMNI_1000_ParentsList.txt

# GENERATE A LIST OF 1000 PARENTS FROM GSA DATA (PARENTS LIST GENERATED FROM EXCEL SAMPLE SHEET)
ls -d \
    "$WKD"/Intensities/TS/* > \
    $WKD/SupportingFiles/GSA1_Samples.List

sort $WKD/SupportingFiles/GSA1_ParentsList.txt | \
    uniq -u | \
    shuf -n 1000 > \
    $WKD/SupportingFiles/GSA1_1000_ParentsList.txt

sed -i -e 's/^/\/Intensities\/TS\//' $WKD/SupportingFiles/GSA1_1000_ParentsList.txt

# GENERATE A LIST OF 1000 SAMPLES FROM OEE DATA
ls -d \
    "$WKD"/Intensities/OEE/* > \
    $WKD/SupportingFiles/OEE_Samples.List

shuf -n 1000 \
    $WKD/SupportingFiles/OEE_Samples.List > \
    $WKD/SupportingFiles/OEE_1000_SamplesList.txt

# CREATE A FULL SET PFB FILE (MARKER MATCH WILL SUBSET IT)
sbatch $WKD/Scripts/PFB_GSA1.sh
sbatch $WKD/Scripts/PFB_OMNI.sh
sbatch $WKD/Scripts/PFB_OEE.sh
```

### Adjust Intensities for Genomic Wave
On the cloud, GC-adjust intensity files, used for CNV calling. Scripts for [GCA_GSA1.sh](Scripts/GCA_GSA1.sh), [GCA_OMNI.sh](Scripts/GCA_OMNI.sh), and [GCA_OEE.sh](Scripts/GCA_OEE.sh) are provided.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# CREATE HOLDING DIRECTORIES
mkdir -p $WKD/Adjusted/TS $WKD/Adjusted/ASD $WKD/Adjusted/OEE

# RUN GENOMIC WAVE ADJUSTMENT
sbatch $WKD/Scripts/GCA_GSA1.sh
sbatch $WKD/Scripts/GCA_OMNI.sh
sbatch $WKD/Scripts/GCA_OEE.sh

# GENERATE ADJUSTED LIST OMNI SAMPLES
ls -d \
    "$WKD"/Adjusted/ASD/* > \
    $WKD/SupportingFiles/OMNI_Samples.List

# GENERATE ADJUSTED LIST OF GSA SAMPLES
ls -d \
    "$WKD"/Adjusted/TS/* > \
    $WKD/SupportingFiles/GSA1_Samples.List

# GENERATE ADJUSTED LIST OF OEE SAMPLES
ls -d \
    "$WKD"/Adjusted/OEE/* > \
    $WKD/SupportingFiles/OEE_Samples.List
```

### Generate Marker Match Specific Supporting Files
On the cloud, use [MarkerMatch_Extract_SSC.r](Scripts/MarkerMatch_Extract_SSC.r) and [MarkerMatch_Extract_OEE.r](Scripts/MarkerMatch_Extract_OEE.r) to generate supporting files for each MarkerMatch configuration used in this validation study.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# RUN R SCRIPT TO EXTRACT MARKERMATCHED DATA
module load R
Rscript --vanilla $WKD/Scripts/MarkerMatch_Extract_SSC.r
Rscript --vanilla $WKD/Scripts/MarkerMatch_Extract_OEE.r
```

### Call Validation CNVs
On the cloud, run jobs to create MarkerMatch CNV callsets. Scripts for [MarkerMatch_Caller5.sh](Scripts/MarkerMatch_Caller3.sh), [MarkerMatch_Caller5.sh](Scripts/MarkerMatch_Caller4.sh), and [MarkerMatch_Caller5.sh](Scripts/MarkerMatch_Caller5.sh) are provided.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# FRAGMENT THE SAMPLE FILES
split -l 710 --numeric-suffixes $WKD/SupportingFiles/OMNI_Samples.List $WKD/SupportingFiles/OMNI_Samples.List
split -l 654 --numeric-suffixes $WKD/SupportingFiles/OEE_Samples.List $WKD/SupportingFiles/OEE_Samples.List
split -l 732 --numeric-suffixes $WKD/SupportingFiles/GSA1_Samples.List $WKD/SupportingFiles/GSA1_Samples.List

# RUN FULL SET AND PERFECT MATCH
for j in {0..5}; do bash Scripts/MarkerMatch_Caller3.sh FullSet 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller3.sh PerfectMatch 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller4.sh OEE_FullSet 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller4.sh OEE_PerfectMatch_0 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller5.sh GSA_FullSet 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller5.sh GSA_PerfectMatch_0 0$j; done;

# RUN THE MARKER MATCHED SAMPLES
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do for j in {0..5}; do bash Scripts/MarkerMatch_Caller3.sh BAF_$i 0$j; done; done;
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do for j in {0..5}; do bash Scripts/MarkerMatch_Caller3.sh LRRmean_$i 0$j; done; done;
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do for j in {0..5}; do bash Scripts/MarkerMatch_Caller3.sh LRRsd_$i 0$j; done; done;
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do for j in {0..5}; do bash Scripts/MarkerMatch_Caller3.sh Pos_$i 0$j; done; done;

# RUN THE MARKER MATCHED VALIDATION 2 SAMPLES OEE
for j in {0..5}; do bash Scripts/MarkerMatch_Caller4.sh OEE_Pos_10000 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller4.sh OEE_BAF_10000 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller4.sh OEE_LRRsd_10000 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller4.sh OEE_LRRmean_10000 0$j; done;

# RUN THE MARKER MATCHED VALIDATION 2 SAMPLES GSA
for j in {0..5}; do bash Scripts/MarkerMatch_Caller5.sh GSA_Pos_10000 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller5.sh GSA_BAF_10000 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller5.sh GSA_LRRsd_10000 0$j; done;
for j in {0..5}; do bash Scripts/MarkerMatch_Caller5.sh GSA_LRRmean_10000 0$j; done;
```

### Preprocess Validation CNVs
On the cloud, merge fragmented CNV callsets.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# DEFINE FUNCTION
CNV_merger () {
	cat Validation/$1/Log/$2.Log_0* > Validation/$1/Log/$2.Log
    cat Validation/$1/Raw/$2.Raw_0* > Validation/$1/Raw/$2.Raw
}

# RUN FUNCTIONS
CNV_merger SSC FullSet
CNV_merger SSC PerfectMatch

CNV_merger SSC BAF_10       # repeat for all BAF configurations
CNV_merger SSC LRRmean_10   # repeat for all LRRmean configurations
CNV_merger SSC LRRsd_10     # repeat for all LRRsd configurations
CNV_merger SSC Pos_10       # repeat for all Pos configurations

CNV_merger GSA GSA_FullSet
CNV_merger GSA GSA_PerfectMatch_0
CNV_merger GSA GSA_BAF_10000       
CNV_merger GSA GSA_LRRmean_10000   
CNV_merger GSA GSA_LRRsd_10000    
CNV_merger GSA GSA_Pos_10000

CNV_merger OEE OEE_FullSet
CNV_merger OEE OEE_PerfectMatch_0
CNV_merger OEE OEE_BAF_10000       
CNV_merger OEE OEE_LRRmean_10000   
CNV_merger OEE OEE_LRRsd_10000    
CNV_merger OEE OEE_Pos_10000
```

### Clean Validation CNVs
On the cloud, run scripts to clean up CNV callsets. Scripts for [MarkerMatch_Clean.sh](Scripts/MarkerMatch_Clean.sh), [MarkerMatch_Clean2.sh](Scripts/MarkerMatch_Clean2.sh), and [MarkerMatch_Clean3.sh](Scripts/MarkerMatch_Clean3.sh) are provided.

```bash
# ASSIGN PATHS TO VARIABLE NAMES
WKD="..."
PCN="/apps/penncnv/1.0.5"

# RUN FULL SET AND PERFECT MATCH
bash Scripts/MarkerMatch_Clean.sh FullSet
bash Scripts/MarkerMatch_Clean.sh PerfectMatch
bash Scripts/MarkerMatch_Clean2.sh OEE_FullSet
bash Scripts/MarkerMatch_Clean2.sh OEE_PerfectMatch_0
bash Scripts/MarkerMatch_Clean3.sh GSA_FullSet
bash Scripts/MarkerMatch_Clean3.sh GSA_PerfectMatch_0

# RUN THE MARKER MATCHED SAMPLES
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do bash Scripts/MarkerMatch_Clean.sh BAF_$i; done;
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do bash Scripts/MarkerMatch_Clean.sh LRRmean_$i; done;
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do bash Scripts/MarkerMatch_Clean.sh LRRsd_$i; done;
for i in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do bash Scripts/MarkerMatch_Clean.sh Pos_$i; done;

bash Scripts/MarkerMatch_Clean2.sh OEE_BAF_10000
bash Scripts/MarkerMatch_Clean2.sh OEE_LRRmean_10000
bash Scripts/MarkerMatch_Clean2.sh OEE_LRRsd_10000
bash Scripts/MarkerMatch_Clean2.sh OEE_Pos_10000

bash Scripts/MarkerMatch_Clean3.sh GSA_BAF_10000
bash Scripts/MarkerMatch_Clean3.sh GSA_LRRmean_10000
bash Scripts/MarkerMatch_Clean3.sh GSA_LRRsd_10000
bash Scripts/MarkerMatch_Clean3.sh GSA_Pos_10000
```

### Prepare Callsets for Analysis
On the cloud, run the following R code to prepare callsets for analysis.
```bash
# BOOT UP AN INTERACTIVE SESSION
srun --mem=32gb --time=08:00:00 --pty bash -i

# BOOT UP R
module load R
R
```
```R
# SET UP WORKING DIRECTORIES
WORK_DIR="/Validation"
dir.create(paste0(WORK_DIR,"/Analysis"))
setwd("Validation")

# LOAD LIBRARIES
library(tidyverse)

# INITIALIZE HOLDING DATAFRAMES
DataQC <- data.frame(
    lines=as.character(NULL),
    stringsAsFactors=FALSE)

DataCNV <- data.frame(
    lines=as.character(NULL),
    stringsAsFactors=FALSE)

# CREATE DATAFRAME OF INPUT LISTS
DataFileNames <- data.frame(
    Subdirectory=c(rep("SSC", 50), rep("GSA", 6), rep("OEE", 6)),
    Factor=c("FullSet", "PerfectMatch", rep("BAF", 12), rep("LRRmean", 12), rep("LRRsd", 12), rep("Pos", 12), "GSA_FullSet", "GSA_PerfectMatch", "GSA_BAF", "GSA_LRRmean", "GSA_LRRsd", "GSA_Pos", "OEE_FullSet", "OEE_PerfectMatch", "OEE_BAF", "OEE_LRRmean", "OEE_LRRsd", "OEE_Pos"),
    MaxD=c(NA, NA, rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4), NA, "0", rep("10000", 4),  NA, "0", rep("10000", 4)),
    FactorLab=c("FullSet", "PerfectMatch", rep("BAF", 12), rep("LRRmean", 12), rep("LRRsd", 12), rep("Pos", 12), rep(c("FullSet", "PerfectMatch", "BAF", "LRRmean", "LRRsd", "Pos"), 2)),
    MaxDLab=c("0", "0", rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4), "0", "0", rep("10000", 4),  "0", "0", rep("10000", 4)),
    stringsAsFactors=FALSE)

# LOAD CNV AND QC DATA INTO A LIST
DATA <- list()

# IMPORT DATA
WORK_DIR <- getwd()
for(i in 1:nrow(DataFileNames)){
    FIL_NAME <- ifelse(is.na(DataFileNames$MaxD[i]), DataFileNames$Factor[i], paste0(DataFileNames$Factor[i], "_", DataFileNames$MaxD[i]))

    con <- file(paste0(WORK_DIR, "/", DataFileNames$Subdirectory[i], "/Log/", FIL_NAME, ".Log"), open="r")
    DATA[["Raw"]][[DataFileNames$Subdirectory[i]]][[DataFileNames$FactorLab[i]]][[DataFileNames$MaxDLab[i]]][["QC"]] <- data.frame(lines=readLines(con))
    close(con)

    DATA[["Raw"]][[DataFileNames$Subdirectory[i]]][[DataFileNames$FactorLab[i]]][[DataFileNames$MaxDLab[i]]][["CNV"]] <- read.table(paste0(WORK_DIR, "/", DataFileNames$Subdirectory[i], "/Mer/", FIL_NAME, ".Mer"), sep="", col.names=c("Position", "NumSNP", "Length", "StateCopyNumber", "ID", "StartSNP", "EndSNP", "Confidence"), skip=0, blank.lines.skip=TRUE, header=FALSE, stringsAsFactors=FALSE)

    DATA[["Raw"]][[DataFileNames$Subdirectory[i]]][[DataFileNames$FactorLab[i]]][[DataFileNames$MaxDLab[i]]][["Cen"]] <- read.table(paste0(WORK_DIR, "/", DataFileNames$Subdirectory[i], "/Extra/", FIL_NAME, ".Centromere"), sep="", col.names=c("Position", "NumSNP", "Length", "StateCopyNumber", "ID", "StartSNP", "EndSNP", "Confidence"), skip=0, blank.lines.skip=TRUE, header=FALSE, stringsAsFactors=FALSE)

    DATA[["Raw"]][[DataFileNames$Subdirectory[i]]][[DataFileNames$FactorLab[i]]][[DataFileNames$MaxDLab[i]]][["Imu"]] <- read.table(paste0(WORK_DIR, "/", DataFileNames$Subdirectory[i], "/Extra/", FIL_NAME, ".Imumune"), sep="", col.names=c("Position", "NumSNP", "Length", "StateCopyNumber", "ID", "StartSNP", "EndSNP", "Confidence"), skip=0, blank.lines.skip=TRUE, header=FALSE, stringsAsFactors=FALSE)

    DATA[["Raw"]][[DataFileNames$Subdirectory[i]]][[DataFileNames$FactorLab[i]]][[DataFileNames$MaxDLab[i]]][["Tel"]] <- read.table(paste0(WORK_DIR, "/", DataFileNames$Subdirectory[i], "/Extra/", FIL_NAME, ".Telomere"), sep="", col.names=c("Position", "NumSNP", "Length", "StateCopyNumber", "ID", "StartSNP", "EndSNP", "Confidence"), skip=0, blank.lines.skip=TRUE, header=FALSE, stringsAsFactors=FALSE)

    DATA[["Raw"]][[DataFileNames$Subdirectory[i]]][[DataFileNames$FactorLab[i]]][[DataFileNames$MaxDLab[i]]][["SegDup"]] <- read.table(paste0(WORK_DIR, "/", DataFileNames$Subdirectory[i], "/Extra/", FIL_NAME, ".SegDup"), sep="", col.names=c("Position", "NumSNP", "Length", "StateCopyNumber", "ID", "StartSNP", "EndSNP", "Confidence"), skip=0, blank.lines.skip=TRUE, header=FALSE, stringsAsFactors=FALSE)
}

# CHECKPOINT SAVE
save(DATA, file=paste0(WORK_DIR,"/Analysis/Validation_Import_082924.RData"))
load(paste0(WORK_DIR,"/Analysis/Validation_Import_082924.RData"))
```
### Prepare Samples for Analysis
Locally, run the following code to prepare samples for the analysis.

```bash
# SET UP WORKING DIRECTORY
WKDIR="..."
cd $WKDIR

# CONVERT PLINK FILES INTO BFILES
cd $WKDIR/TS_TAAICG_GSA_1
./plink --file TS_TAAICG_GSA-MDv1_Wave1_C --make-bed --out GSA_S01

cd $WKDIR/TS_TAAICG_GSA_2
./plink --file TS_TAAICG_GSA-MDv1_Wave2_C --make-bed --out GSA_S02

cd $WKDIR/TS_TAAICG_GSA_3
./plink --file TS_TAAICG_GSA-MDv1_Wave3_C --make-bed --out GSA_S03

cd $WKDIR/TS_TAAICG_OEE_1
./plink --file 2024_OEE_033s --make-bed --out OEE_S01

cd $WKDIR/TS_TAAICG_OEE_2
./plink --file 2024_OEE_328s --make-bed --out OEE_S02

# MOVE ALL CONVERTED FILES INTO A SINGLE PROCESSING FOLDER
mkdir $WKDIR/PROCESSING
mv $WKDIR/TS_TAAICG_GSA_1/GSA_S01*  $WKDIR/PROCESSING
mv $WKDIR/TS_TAAICG_GSA_2/GSA_S02*  $WKDIR/PROCESSING
mv $WKDIR/TS_TAAICG_GSA_3/GSA_S03*  $WKDIR/PROCESSING
mv $WKDIR/TS_TAAICG_OEE_1/OEE_S01*  $WKDIR/PROCESSING
mv $WKDIR/TS_TAAICG_OEE_2/OEE_S02*  $WKDIR/PROCESSING

# COMBINE GSA FILES THEN OEE FILES
cd $WKDIR/PROCESSING
./plink --merge-list GSA_LIST.txt --make-bed --out GSA
./plink --merge-list OEE_LIST.txt --make-bed --out OEE

# PERFORM BASIC QC ON THE PLINK DATASETS
./plink --bfile GSA --geno 0.02 --maf 0.01 --hwe 0.000001 --make-bed --out GSA_QC
./plink --bfile OEE --geno 0.02 --maf 0.01 --hwe 0.000001 --make-bed --out OEE_QC

# WRITE OUT MAF INFO FOR GSA AND OEE
./plink --bfile GSA --freqx --out GSA_MAF
./plink --bfile OEE --freqx --out OEE_MAF
```

Open R, and use following code to determine which alleles to keep from both OEE and GSA sets for further analysis.
```R
# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)

# IMPORT AND PREPARE DATA
GSA_SNP <- read_tsv("GSA.bim", col_names=FALSE) %>%
  rename(Chr=X1, SNP=X2, Start=X3, End=X4, A1=X5, A2=X6)

OEE_SNP <- read_tsv("OEE.bim", col_names=FALSE) %>%
  rename(Chr=X1, SNP=X2, Start=X3, End=X4, A1=X5, A2=X6)

GSA_MAF <- read_tsv("GSA_MAF.frqx", col_names=TRUE) %>%
  mutate(C_ALL = 2 * (`C(HOM A1)` + `C(HET)` + `C(HOM A2)` + `C(HAP A1)` + `C(HAP A2)` + `C(MISSING)`)) %>%
  mutate(A1_F = (2 * `C(HOM A1)` + `C(HET)`)/C_ALL,
         A2_F = (2 * `C(HOM A2)` + `C(HET)`)/C_ALL) %>%
  select(c(CHR, SNP, A1, A2, A1_F, A2_F))

OEE_MAF <- read_tsv("OEE_MAF.frqx", col_names=TRUE) %>%
  mutate(C_ALL = 2 * (`C(HOM A1)` + `C(HET)` + `C(HOM A2)` + `C(HAP A1)` + `C(HAP A2)` + `C(MISSING)`)) %>%
  mutate(A1_F = (2 * `C(HOM A1)` + `C(HET)`)/C_ALL,
         A2_F = (2 * `C(HOM A2)` + `C(HET)`)/C_ALL) %>%
  select(c(CHR, SNP, A1, A2, A1_F, A2_F))

GSA <- full_join(GSA_SNP, GSA_MAF, by=c("Chr"="CHR", "SNP"="SNP", "A1"="A1", "A2"="A2"))
OEE <- full_join(OEE_SNP, OEE_MAF, by=c("Chr"="CHR", "SNP"="SNP", "A1"="A1", "A2"="A2"))

rm(GSA_MAF, GSA_SNP, OEE_MAF, OEE_SNP)

# GET A LIST OF SNPS WITH SAME POSITIONS AND SAME IDS AND SAME ALLELES
SNP_IID_IP_NF <- OEE %>%
  inner_join(GSA, by=c("Chr", "SNP", "A1", "A2", "Start", "End"), keep=TRUE, suffix=c(".OEE", ".GSA"))

TOSS_SNP <- SNP_IID_IP_NF$SNP.OEE

# GET A LIST OF SNPS WITH SAME POSITIONS AND SAME IDS AND FLPPED ALLELES
SNP_IID_IP_YF <- OEE %>%
  inner_join(GSA, by=c("Chr", "SNP", "A1"="A2", "A2"="A1", "Start", "End"), keep=TRUE, suffix=c(".OEE", ".GSA"))

TOSS_SNP <- c(TOSS_SNP, SNP_IID_IP_YF$SNP.OEE)

# GET A LIST OF SNPs WITH SAME ID BUT NOT FLIPPED ALLELES
SNP_IID_NF <- OEE %>%
  inner_join(GSA, by=c("Chr", "SNP", "A1", "A2"), keep=TRUE, suffix=c(".OEE", ".GSA")) %>%
  filter(!SNP.OEE %in% TOSS_SNP)

TOSS_SNP <- c(TOSS_SNP, SNP_IID_NF$SNP.OEE)

# GET A LIST OF SNPs WITH SAME ID BUT FLIPPED ALLELES
SNP_IID_YF <- OEE %>%
  inner_join(GSA, by=c("Chr", "SNP", "A1"="A2", "A2"="A1"), keep=TRUE, suffix=c(".OEE", ".GSA")) %>%
  filter(!SNP.OEE %in% TOSS_SNP)

TOSS_SNP <- c(TOSS_SNP, SNP_IID_YF$SNP.OEE)

# GET A LIST OF SNPS WITH DIFFERNT IDs SAME POSITIONS AND SAME ALLELES
SNP_DID_IP_NF <- OEE %>%
  inner_join(GSA, by=c("Chr", "Start", "End", "A1", "A2"), keep=TRUE, suffix=c(".OEE", ".GSA")) %>%
  filter(!SNP.OEE %in% TOSS_SNP) %>%
  mutate(SNP.GSA.V2=gsub("GSA-", "", SNP.GSA)) %>%
  filter(SNP.OEE == SNP.GSA.V2)

# GET A LIST OF SNPS WITH DIFFERNT IDs SAME POSITIONS AND FLIPPED ALLELES
SNP_DID_IP_YF <- OEE %>%
  inner_join(GSA, by=c("Chr", "Start", "End", "A1"="A2", "A2"="A1"), keep=TRUE, suffix=c(".OEE", ".GSA")) %>%
  filter(!SNP.OEE %in% TOSS_SNP) %>%
  mutate(SNP.GSA.V2=gsub("GSA-", "", SNP.GSA)) %>%
  filter(SNP.OEE == SNP.GSA.V2)

# USE THE 160,487 SNPs WITH SAME ID AND SAME POSITION AND SAME (BUT PERHPAS FLIPPED) ALLELES
# WRITE NECESSARY FILES (TO FILTER, UPDATE ID, OR FLIP ALLELES)
bind_rows(SNP_IID_IP_NF, SNP_IID_IP_YF, SNP_DID_IP_NF, SNP_DID_IP_YF) %>%
  select(SNP.OEE) %>%
  distinct() %>%
  write_delim("KEEP_SNP_OEE.txt", col_names=FALSE)

bind_rows(bind_rows(SNP_IID_IP_NF, SNP_IID_IP_YF) %>% select(SNP.GSA),
  bind_rows(SNP_DID_IP_NF, SNP_DID_IP_YF) %>% select(SNP.GSA.V2) %>% rename(SNP.GSA=SNP.GSA.V2)) %>%
  distinct() %>%
  write_delim("KEEP_SNP_GSA.txt", col_names=FALSE)

bind_rows(SNP_DID_IP_NF, SNP_DID_IP_YF) %>% 
  select(c(SNP.GSA, SNP.GSA.V2)) %>% 
  write_delim("UPDATE_SNPID_GSA.txt", col_names=FALSE)

bind_rows(SNP_IID_IP_YF %>% select(SNP.GSA),
          SNP_DID_IP_YF %>% select(SNP.GSA.V2) %>% rename(SNP.GSA=SNP.GSA.V2)) %>%
  distinct() %>%
  write_delim("FLIP_SNP_GSA.txt", col_names=FALSE)
```

Back in terminal, use Plink to create a merged OEE/GSA dataset.

```bash
# SET UP WORKING DIRECTORY
WKDIR="..."
cd $WKDIR

# UPDATE IDS FOR GSA
./plink --bfile GSA --update-name UPDATE_SNPID_GSA.txt --make-bed --out GSA_UID

# FLIP STRANDS FOR GSA
./plink --bfile GSA_UID --flip FLIP_SNP_GSA.txt --make-bed --out GSA_FLIP

# FILTER SNPS FROM BOTH GSA AND OEE
./plink --bfile GSA_FLIP --extract KEEP_SNP_GSA.txt --make-bed --out GSA_KEEP
./plink --bfile OEE --extract KEEP_SNP_OEE.txt --make-bed --out OEE_KEEP

# TRIAL MERGE GSA AND OEE
./plink --bfile OEE_KEEP --bmerge GSA_KEEP.bed GSA_KEEP.bim GSA_KEEP.fam --make-bed --out GSAOEE_TRIAL

# EXCLUDE PROBLEMATIC SNPS
./plink --bfile GSA_KEEP --exclude GSAOEE_TRIAL-merge.missnp --make-bed --out GSA_EXCL
./plink --bfile OEE_KEEP --exclude GSAOEE_TRIAL-merge.missnp --make-bed --out OEE_EXCL

# CROSS-CONTAMINATION AND TWIN CHECK IN INDIVIDUAL DATASETS
./plink --bfile GSA_EXCL  --genome --min 0.15 --out GSA_TWINS
./plink --bfile OEE_EXCL  --genome --min 0.15 --out OEE_TWINS
```

Back in R, use the output of Plink `--genom` option to identify samples in GSA and OEE that need to be excluded from the analysis due to cross-contamination and/or twin status.

```R
# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)

# IMPORT DATASET
GENO_GSA <- read.csv("GSA_TWINS.genome", sep="", header=TRUE)
GENO_OEE <- read.csv("OEE_TWINS.genome", sep="", header=TRUE)
  
## GENO STRUCTURE
# FID1	Family ID for first sample
# IID1	Individual ID for first sample
# FID2	Family ID for second sample
# IID2	Individual ID for second sample
# RT	Relationship type inferred from .fam/.ped file
# EZ	IBD sharing expected value, based on just .fam/.ped relationship
# Z0	P(IBD=0)
# Z1	P(IBD=1)
# Z2	P(IBD=2)
# PI_HAT	Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
# PHE	Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)
# DST	IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)
# PPC	IBS binomial test
# RATIO	HETHET : IBS0 SNP ratio (expected value 2)

# CHECK FOR CROSS-CONTAMINATION
COUNTS_GSA <- as.data.frame.table(table(c(GENO_GSA$IID1, GENO_GSA$IID2))) %>% arrange(desc(Freq))
COUNTS_OEE <- as.data.frame.table(table(c(GENO_OEE$IID1, GENO_OEE$IID2))) %>% arrange(desc(Freq))

# REMOVE ANYONE WITH 30+ matches
REMOVE_GSA <- COUNTS_GSA %>%
  filter(Freq >= 30) %>%
  pull(Var1) %>%
  as.character(.)

REMOVE_OEE <- COUNTS_OEE %>%
  filter(Freq >= 30) %>%
  pull(Var1) %>%
  as.character(.)

# IDENTIFY TWINS
TWINS_GSA <- GENO_GSA %>%
  filter(Z2 > 0.9) %>%
  unique() %>%
  filter(!IID1 %in% REMOVE_GSA | !IID2 %in% REMOVE_GSA)

TWINS_GSA <- c(TWINS_GSA$IID1, TWINS_GSA$IID2) %>%
  unique()

TWINS_OEE <- GENO_OEE %>%
  filter(Z2 > 0.9) %>%
  unique() %>%
  filter(!IID1 %in% REMOVE_OEE | !IID2 %in% REMOVE_OEE)

TWINS_OEE <- c(TWINS_OEE$IID1, TWINS_OEE$IID2) %>%
  unique()

# SAVE THE FILE FOR FILTERING
bind_rows(GENO_GSA %>% select(c(FID1, IID1)),
          GENO_GSA %>% select(c(FID2, IID2)) %>% rename(FID1=FID2, IID1=IID2)) %>%
  filter(IID1 %in% c(REMOVE_GSA, TWINS_GSA)) %>%
  distinct() %>%
  write_tsv("REMOVE_SAMPLE_GSA.txt", col_names=FALSE)

bind_rows(GENO_OEE %>% select(c(FID1, IID1)),
          GENO_OEE %>% select(c(FID2, IID2)) %>% rename(FID1=FID2, IID1=IID2)) %>%
  filter(IID1 %in% c(REMOVE_OEE, TWINS_OEE)) %>%
  distinct() %>%
  write_tsv("REMOVE_SAMPLE_OEE.txt", col_names=FALSE)
```

Back in terminal, use Plink to remove twins and potential cross-contaminations, and identify duplicate samples across GSA and OEE.

```bash
# SET UP WORKING DIRECTORY
WKDIR="..."
cd $WKDIR

# REMOVE CROSS-CONTAMINATED SAMPLES AND TWINS
./plink --bfile GSA_EXCL  --remove REMOVE_SAMPLE_GSA.txt --make-bed --out GSA_CLN
./plink --bfile OEE_EXCL  --remove REMOVE_SAMPLE_OEE.txt --make-bed --out OEE_CLN

# MERGE GSA AND OEE
./plink --bfile OEE_CLN --bmerge GSA_CLN.bed GSA_CLN.bim GSA_CLN.fam --make-bed --out GSAOEE

# PERFORM BASIC QC ON GSAOEE
./plink --bfile GSAOEE --geno 0.02 --maf 0.01 --hwe 0.000001 --make-bed --out GSAOEE_QC

# IDENTIFY GENETIC RELATIONSHIPS
./plink --bfile GSAOEE_QC --genome --min 0.9 --out GSAOEE_RELATEDNESS
```

Back in R, use the output of Plink `--genom` option to identify samples that have been genotyped at both GSA and OEE platforms.

```R
# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)

# IMPORT DATASET

GENO <- read.csv("GSAOEE_RELATEDNESS.genome", sep="", header=TRUE) %>%
  filter(Z2 > 0.99)

FAM_GSA <- read.csv("GSA_CLN.fam", sep="", header=FALSE) %>%
  select(V2)

FAM_OEE <- read.csv("OEE_CLN.fam", sep="", header=FALSE) %>%
  select(V2)

## GENO STRUCTURE
# FID1	Family ID for first sample
# IID1	Individual ID for first sample
# FID2	Family ID for second sample
# IID2	Individual ID for second sample
# RT	Relationship type inferred from .fam/.ped file
# EZ	IBD sharing expected value, based on just .fam/.ped relationship
# Z0	P(IBD=0)
# Z1	P(IBD=1)
# Z2	P(IBD=2)
# PI_HAT	Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
# PHE	Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)
# DST	IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)
# PPC	IBS binomial test
# RATIO	HETHET : IBS0 SNP ratio (expected value 2)

# CREATE ANNOTATION TABLE
ANNOT <- bind_rows(
  FAM_GSA %>%
    rename(ID=V2) %>%
    mutate(ARRAY="GSA"),
  FAM_OEE %>%
    rename(ID=V2) %>%
    mutate(ARRAY="OEE")) %>%
  filter(ID %in% c(GENO$IID1, GENO$IID2)) %>%
  select(c(ID, ARRAY))

# ANNOTATE GENO WITH ARRAY INFORMATION
DUPLICATES <- GENO %>%
  select(c(IID1, IID2)) %>%
  left_join(ANNOT, by=c("IID1"="ID")) %>%
  rename(ARR1=ARRAY) %>%
  left_join(ANNOT, by=c("IID2"="ID")) %>%
  rename(ARR2=ARRAY) %>%
  mutate(ID_GSA=ifelse(ARR1=="GSA", IID1, IID2),
         ID_OEE=ifelse(ARR1=="OEE", IID1, IID2))

# WRITE TABLE FOR FUTURE USE
DUPLICATES %>%
  select(c(ID_GSA, ID_OEE)) %>%
  write_tsv("TS_DUPLICATES.tsv", col_names=TRUE)
```

### Prepare CNV Data for Final Analyses
Download the `Validation_Import_082924.RData` from the Cloud and import it into local R session, analyzing the outputs.

```R
# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)

# LOAD DATA
load("Validation_Import_082924.RData")

# PROCESS IMPORTED DATA
QC_Curate <- "NOTICE: quality summary for /Adjusted/ASD/|NOTICE: quality summary for /Adjusted/OEE/|NOTICE: quality summary for /Adjusted/TS|_Omni2.5_FinalReport.txt:|LRR_mean=|LRR_median=|LRR_SD=|BAF_mean=|BAF_median=|BAF_SD=|BAF_DRIFT=|WF=|GCWF=|_Omni2.5_FinalReport.txt|/|numsnp=|length=|state|cn=|startsnp=|endsnp=|conf=|/Adjusted/ASD/|/Adjusted/OEE/|/Adjusted/TS/|.adjusted:|.adjusted"
QC_ColNames <- c("ID", "LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_drift", "WF", "GCWF")
CNV_ColNames <- c("Position", "N_SNP", "LEN", "State", "CN", "ID", "start_snp", "end_snp", "CONF")

# LOAD DATA
load("Validation_Import_082924.RData")

# PROCESS IMPORTED DATA
QC_Curate <- "NOTICE: quality summary for /blue/carolmathews/njofrica/Adjusted/ASD/|NOTICE: quality summary for /blue/carolmathews/njofrica/Adjusted/OEE/|NOTICE: quality summary for /blue/carolmathews/njofrica/Adjusted/TS|_Omni2.5_FinalReport.txt:|LRR_mean=|LRR_median=|LRR_SD=|BAF_mean=|BAF_median=|BAF_SD=|BAF_DRIFT=|WF=|GCWF=|_Omni2.5_FinalReport.txt|/|numsnp=|length=|state|cn=|startsnp=|endsnp=|conf=|/blue/carolmathews/njofrica/Adjusted/ASD/|/blue/carolmathews/njofrica/Adjusted/OEE/|/blue/carolmathews/njofrica/Adjusted/TS/|.adjusted:|.adjusted"
QC_ColNames <- c("ID", "LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_drift", "WF", "GCWF")
CNV_ColNames <- c("Position", "N_SNP", "LEN", "State", "CN", "ID", "start_snp", "end_snp", "CONF")

# CREATE DATAFRAME OF INPUT LISTS
DataFileNames <- data.frame(
  Subdirectory=c(rep("SSC", 50), rep("GSA", 6), rep("OEE", 6)),
  Factor=c("FullSet", "PerfectMatch", rep("BAF", 12), rep("LRRmean", 12), rep("LRRsd", 12), rep("Pos", 12), "GSA_FullSet", "GSA_PerfectMatch", "GSA_BAF", "GSA_LRRmean", "GSA_LRRsd", "GSA_Pos", "OEE_FullSet", "OEE_PerfectMatch", "OEE_BAF", "OEE_LRRmean", "OEE_LRRsd", "OEE_Pos"),
  MaxD=c(NA, NA, rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4), NA, "0", rep("10000", 4),  NA, "0", rep("10000", 4)),
  FactorLab=c("FullSet", "PerfectMatch", rep("BAF", 12), rep("LRRmean", 12), rep("LRRsd", 12), rep("Pos", 12), rep(c("FullSet", "PerfectMatch", "BAF", "LRRmean", "LRRsd", "Pos"), 2)),
  MaxDLab=c("0", "0", rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4), "0", "0", rep("10000", 4),  "0", "0", rep("10000", 4)),
  stringsAsFactors=FALSE)

Process_QC <- function(DF) {
  DF %>%
    filter(grepl("quality summary", lines)) %>%
    mutate(lines=gsub(QC_Curate, "", lines)) %>%
    separate(lines, sep="\\s", into=QC_ColNames) %>%
    mutate_at(vars(LRR_mean:GCWF), as.numeric)
}

Process_CNV <- function(DF) {
  DF %>%
    mutate_at(.vars=1:8, .fun=gsub, pattern=QC_Curate, replacement="") %>%
    rowwise() %>%
    mutate(ID=ID,
           CHR=gsub("chr","",str_split(Position, ":|-")[[1]][1]),
           START=str_split(Position, ":|-")[[1]][2],
           END=str_split(Position, ":|-")[[1]][3],
           STATE=str_split(StateCopyNumber, ",")[[1]][1],
           CN=str_split(StateCopyNumber, ",")[[1]][2],
           .keep="unused") %>%
    ungroup() %>%
    mutate(CNV_ID=paste0(ID, "-", CHR, "-", START, "-", END, "-", STATE, "-", CN))
}

for(h in 1:nrow(DataFileNames)){
  k <- DataFileNames$Subdirectory[h]
  i <- DataFileNames$FactorLab[h]
  j <- DataFileNames$MaxDLab[h]
  
  message(k, "\t", i, "\t", j)
  
  DATA[["Raw"]][[k]][[i]][[j]][["QC"]] <- DATA[["Raw"]][[k]][[i]][[j]][["QC"]] %>% Process_QC(.)
  DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] <- DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] %>% Process_CNV(.)
}

for(h in 1:nrow(DataFileNames)){
  k <- DataFileNames$Subdirectory[h]
  i <- DataFileNames$FactorLab[h]
  j <- DataFileNames$MaxDLab[h]
  
  message(k, "\t", i, "\t", j)
  
  if (nrow(DATA[["Raw"]][[k]][[i]][[j]][["Cen"]]) > 0) {
    DATA[["Raw"]][[k]][[i]][[j]][["Cen"]] <- DATA[["Raw"]][[k]][[i]][[j]][["Cen"]] %>% Process_CNV(.) %>% pull(CNV_ID)
  } else {cat(DataFileNames$Factor[h], "_", DataFileNames$MaxDLab[h], " has no overlaps with centromeres.\n")}
  
  if (nrow(DATA[["Raw"]][[k]][[i]][[j]][["Tel"]]) > 0) {
    DATA[["Raw"]][[k]][[i]][[j]][["Tel"]] <- DATA[["Raw"]][[k]][[i]][[j]][["Tel"]] %>% Process_CNV(.) %>% pull(CNV_ID)
  } else {cat(DataFileNames$Factor[h], "_", DataFileNames$MaxDLab[h], " has no overlaps with telomeres.\n")}
  
  if (nrow(DATA[["Raw"]][[k]][[i]][[j]][["Imu"]]) > 0) {
    DATA[["Raw"]][[k]][[i]][[j]][["Imu"]] <- DATA[["Raw"]][[k]][[i]][[j]][["Imu"]] %>% Process_CNV(.) %>% pull(CNV_ID)
  } else {cat(DataFileNames$Factor[h], "_", DataFileNames$MaxDLab[h], " has no overlaps with immune regions.\n")}
  
  if (nrow(DATA[["Raw"]][[k]][[i]][[j]][["SegDup"]]) > 0) {
    DATA[["Raw"]][[k]][[i]][[j]][["SegDup"]] <- DATA[["Raw"]][[k]][[i]][[j]][["SegDup"]] %>% Process_CNV(.) %>% pull(CNV_ID)
  } else {cat(DataFileNames$Factor[h], "_", DataFileNames$MaxDLab[h], " has no overlaps with segmental duplications.\n")}
}

# CHECKPOINT SAVE
save(DATA, file="Validation_Chekpoint1_082924.RData")
load("Validation_Chekpoint1_082924.RData")

# CONFIDENCE SCORES
CONF_SCORE <- as.numeric(NULL)

for(h in 1:nrow(DataFileNames)){
  k <- DataFileNames$Subdirectory[h]
  i <- DataFileNames$FactorLab[h]
  j <- DataFileNames$MaxDLab[h]
  
  # Annotate overlaps
  CONF_SCORE <- c(CONF_SCORE, DATA[["Raw"]][[k]][[i]][[j]][["CNV"]]$Confidence)
}

quantile(as.numeric(CONF_SCORE), probs=seq(0.1,1, by=0.1))

# STRATIFIERS
for(h in 1:nrow(DataFileNames)){
  k <- DataFileNames$Subdirectory[h]
  i <- DataFileNames$FactorLab[h]
  j <- DataFileNames$MaxDLab[h] 
  
  # Annotate overlaps
  DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] <- DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] %>%
    mutate(Confidence=as.numeric(Confidence)) %>%
    mutate(Tel=ifelse(CNV_ID %in% DATA[["Raw"]][[k]][[i]][[j]][["Tel"]], 1, 0),
           Cen=ifelse(CNV_ID %in% DATA[["Raw"]][[k]][[i]][[j]][["Cen"]], 1, 0),
           Imu=ifelse(CNV_ID %in% DATA[["Raw"]][[k]][[i]][[j]][["Imu"]], 1, 0),
           SegDup=ifelse(CNV_ID %in% DATA[["Raw"]][[k]][[i]][[j]][["SegDup"]], 1, 0),
           CONF_DECILE=case_when(
             Confidence < 11.016 ~ 0,
             Confidence >= 11.016 & Confidence < 13.302 ~ 10,
             Confidence >= 13.302 & Confidence < 14.797 ~ 20,
             Confidence >= 14.797 & Confidence < 16.291 ~ 30,
             Confidence >= 16.291 & Confidence < 18.452 ~ 40,
             Confidence >= 18.452 & Confidence < 21.317 ~ 50,
             Confidence >= 21.317 & Confidence < 25.883 ~ 60,
             Confidence >= 25.883 & Confidence < 33.956 ~ 70,
             Confidence >= 33.956 & Confidence < 54.361 ~ 80,
             Confidence >= 54.361 ~ 90,
             TRUE ~ NA_integer_))
}

# MEDIAN ABSOLUTE DEVIATION
MAD <- function(X=X, direction=direction, scale=1){
  if(direction=="+") {
    return(median(X)+mad(X, constant=scale*1.4826))
  } else if (direction=="-") {
    return(median(X)-mad(X, constant=scale*1.4826))
  } else {
    message("Direction has to be + or -.")
  }
}

# QUALITY CONTROL
for(h in 1:nrow(DataFileNames)){
  k <- DataFileNames$Subdirectory[h]
  i <- DataFileNames$FactorLab[h]
  j <- DataFileNames$MaxDLab[h] 
  
  # Initialize data frames
  DATA[["QCd"]][[k]][[i]][[j]][["CNV"]] <- data.frame(NULL)
  DATA[["QCd"]][[k]][[i]][[j]][["QC"]] <- data.frame(NULL)
  
  # Calculate stats
  LRR_MEAN_UP <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$LRR_mean, direction="+", scale=4)
  LRR_MEAN_DO <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$LRR_mean, direction="-", scale=4)
  LRR_SDEV_UP <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$LRR_SD, direction="+", scale=4)
  BAF_MEAN_UP <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$BAF_mean, direction="+", scale=4)
  BAF_MEAN_DO <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$BAF_mean, direction="-", scale=4)
  BAF_SDEV_UP <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$BAF_SD, direction="+", scale=4)
  BAF_DRIF_UP <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$BAF_drift, direction="+", scale=4)
  GCWF_UP <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$GCWF, direction="+", scale=4)
  GCWF_DO <- MAD(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$GCWF, direction="-", scale=4)
  
  # Get a list of IDs to remove from the dataset from sample QC in Logs
  KEEP_LOG_SAM <- DATA[["Raw"]][[k]][[i]][[j]][["QC"]] %>%
    filter(LRR_mean >= LRR_MEAN_DO & LRR_mean <= LRR_MEAN_UP) %>%
    filter(LRR_SD <= LRR_SDEV_UP) %>%
    filter(BAF_mean >= BAF_MEAN_DO &  BAF_mean <= BAF_MEAN_UP) %>%
    filter(BAF_SD <= BAF_SDEV_UP) %>%
    filter(BAF_drift <= BAF_DRIF_UP) %>%
    filter(GCWF >= GCWF_DO & GCWF <= GCWF_UP) %>%
    pull(ID)
  
  REMOVE_LOG_SAM <- setdiff(DATA[["Raw"]][[k]][[i]][[j]][["QC"]]$ID, KEEP_LOG_SAM)
  
  # Calculate CNV Counts
  CNV_COUNT <- data.frame(table(ID=pull(DATA[["Raw"]][[k]][[i]][[j]][["CNV"]], ID))) %>% 
    arrange(desc(Freq))
  
  CNV_COUNT_LIM <- MAD(CNV_COUNT$Freq, direction="+", scale=4)
  
  # Get a list of IDs to remove from the dataset from sample QC in Raw CNV calls
  REMOVE_CNV_SAM <- CNV_COUNT %>%
    filter(Freq > CNV_COUNT_LIM) %>%
    pull(ID) %>%
    as.character(.)
  
  # Form a list of samples to remove
  REMOVE_SAM <- union(REMOVE_CNV_SAM, REMOVE_LOG_SAM)
  
  # Filter CNVs and QCs to only include the samples passing QC above, additionally, 
  # filter out CNVs spanning less than 10 SNPs and less than 20kb
  DATA[["QCd"]][[k]][[i]][[j]][["CNV"]] <- DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] %>%
    filter(!ID %in% REMOVE_SAM) %>%
    filter(NumSNP >= 10) %>%
    filter(Length >= 20000)
  
  DATA[["QCd"]][[k]][[i]][[j]][["QC"]] <- DATA[["Raw"]][[k]][[i]][[j]][["QC"]] %>%
    filter(!ID %in% REMOVE_SAM)
}

# CHECKPOINT SAVE
save(DATA, file="Validation_Checkpoint2_082924.RData")
load("Validation_Checkpoint2_082924.RData")

# IMPORT SAMPLES TO RENAME/KEEP IN GSA/OEE
SAM_NAM <- read_tsv("TS_DUPLICATES.tsv", col_names=TRUE)

# QUALITY CONTROL
for(h in 1:nrow(DataFileNames)){
  k <- DataFileNames$Subdirectory[h]
  i <- DataFileNames$FactorLab[h]
  j <- DataFileNames$MaxDLab[h]
  
  if(k=="SSC") {next}
  if(k=="OEE") {
    DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] <- DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] %>%
      filter(ID %in% SAM_NAM$ID_OEE)
    DATA[["Raw"]][[k]][[i]][[j]][["QC"]] <- DATA[["Raw"]][[k]][[i]][[j]][["QC"]] %>%
      filter(ID %in% SAM_NAM$ID_OEE)
    DATA[["QCd"]][[k]][[i]][[j]][["CNV"]] <- DATA[["QCd"]][[k]][[i]][[j]][["CNV"]] %>%
      filter(ID %in% SAM_NAM$ID_OEE)
    DATA[["QCd"]][[k]][[i]][[j]][["QC"]] <- DATA[["QCd"]][[k]][[i]][[j]][["QC"]] %>%
      filter(ID %in% SAM_NAM$ID_OEE)
  }
  
  if(k=="GSA") {
    DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] <- DATA[["Raw"]][[k]][[i]][[j]][["CNV"]] %>%
      filter(ID %in% SAM_NAM$ID_GSA) %>%
      left_join(SAM_NAM, by=c("ID"="ID_GSA")) %>%
      select(-ID) %>%
      rename(ID=ID_OEE)
    DATA[["Raw"]][[k]][[i]][[j]][["QC"]] <- DATA[["Raw"]][[k]][[i]][[j]][["QC"]] %>%
      filter(ID %in% SAM_NAM$ID_GSA) %>%
      left_join(SAM_NAM, by=c("ID"="ID_GSA")) %>%
      select(-ID) %>%
      rename(ID=ID_OEE)
    DATA[["QCd"]][[k]][[i]][[j]][["CNV"]] <- DATA[["QCd"]][[k]][[i]][[j]][["CNV"]] %>%
      filter(ID %in% SAM_NAM$ID_GSA) %>%
      left_join(SAM_NAM, by=c("ID"="ID_GSA")) %>%
      select(-ID) %>%
      rename(ID=ID_OEE)
    DATA[["QCd"]][[k]][[i]][[j]][["QC"]] <- DATA[["QCd"]][[k]][[i]][[j]][["QC"]] %>%
      filter(ID %in% SAM_NAM$ID_GSA) %>%
      left_join(SAM_NAM, by=c("ID"="ID_GSA")) %>%
      select(-ID) %>%
      rename(ID=ID_OEE)
  }
}

# CHECKPOINT SAVE
save(DATA, file="Validation_Final_082924.RData")
load("Validation_Final_082924.RData")
```

### Step-One Validation Analyses
Using the `Validation_Final_082924.RData` in local R session, validate CNV calls in simulated MarkerMatch scenario.

Code for step-one genome-wide performance metrics, manuscript **Figure 8-28** is available at [Plot8-28.R](Scripts/Plot8-28.R).
Code for step-one regional performance metrics, manuscript **Figure 29-35** is available at [Plot29-35.R](Scripts/Plot29-35.R).

### Step-Two Validation Analyses
Using the `Validation_Final_082924.RData` in local R session, validate CNV calls in simulated MarkerMatch scenario.

Code for step-two genome-wide performance metrics, manuscript **Figure 37-57** is available at [Plot37-57.R](Scripts/Plot37-57.R).
Code for step-two regional performance metrics, manuscript **Figure 58-64** is available at [Plot58-64.R](Scripts/Plot58-64.R).
Code for step-two determination of optimal LEN and N_SNP cutoffs, manuscript **Figure 65-67** is available at [Plot65-67.R](Scripts/Plot65-67.R), and the output of model summary is shown below.
```
Call:
lm(formula = PPV ~ LEN_Cutoff + N_SNP_Cutoff + Matching_Method, 
    data = filter(ANALYSIS_STEP2_CUTOFFS, !Matching_Method %in% 
        c("FullSet", "PerfectMatch")))

Residuals:
     Min       1Q   Median       3Q      Max 
-0.26643 -0.20000 -0.04473  0.19982  0.23138 

Coefficients:
                         Estimate Std. Error t value Pr(>|t|)    
(Intercept)             7.232e-01  1.679e-02  43.065  < 2e-16 ***
LEN_Cutoff              5.054e-08  3.274e-08   1.544  0.12282    
N_SNP_Cutoff            1.174e-03  3.657e-04   3.209  0.00135 ** 
Matching_MethodLRRmean  1.056e-02  1.335e-02   0.791  0.42920    
Matching_MethodLRRsd   -5.762e-03  1.335e-02  -0.431  0.66617    
Matching_MethodPos     -3.533e-03  1.335e-02  -0.265  0.79136    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2003 on 1794 degrees of freedom
Multiple R-squared:  0.007985,	Adjusted R-squared:  0.005221 
F-statistic: 2.888 on 5 and 1794 DF,  p-value: 0.01328
```
