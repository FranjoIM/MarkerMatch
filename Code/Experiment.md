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
In R, ran locally. Code for creating marker density plots, manuscript **Figure 3**. First, the data cen be prepared using [UsedSNPs.R](Scripts/UsedSNPs.R) script, and then plotted using [Figure3.R](Scripts/Figure3.R).

### Graph Gaps Plots
In R, ran locally. Code for creating marker gaps plots, manuscript **Figure 4**.

### Graph BAF Plots

### Graph LRR-SD Plots

### Graph LRR-Mean Plots

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


