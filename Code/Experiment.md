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

### Complete the Manifest Files

In R, run on a local machine, use thr [ManifSum.R](ManifSum.R) function to suummarize the copy-number metrics (CNM) export into copu-number summary (CNS) file.

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
In R, ran locally. Code for creating marker density plots, manuscript **Figure 3**.

```R

```
