# Set working directory
W_DIR <- "..."
setwd(W_DIR)

# Load packages
library(tidyverse)

# Import array manifests
OMNI_MAN <- read_delim("DATA/OMNI_MAN.txt", delim="\t")
GSA_MAN <- read_delim("DATA/GSA_MAN.txt", delim="\t")
OEE_MAN <- read_delim("DATA/OEE_MAN.txt", delim="\t")

# Process OMNI data
FACTORS <- c("Pos", "BAF", "LRRsd", "LRRmean")
MAXDS <- c("10", "50", "100", "500", "1000", "5000", "10000", "50000", 
           "100000", "500000", "1000000", "5000000")

PerfectMatch <- read_delim("DATA/PerfectMatch_OMNI_MatSelect", delim="\t", col_names=F) %>%
  mutate(PerfectMatch_0=1)

USEDSNPS <- OMNI_MAN %>%
  select(Name) %>%
  mutate(FullSet_0=1) %>%
  left_join(PerfectMatch, by=c("Name"="X1"))

for(i in FACTORS){
  for(j in MAXDS){
    Name <- quo_name(paste0(i, "_", j))
    
    DF_MAT <- read_delim(paste0("DATA/", i, "_", j, "_MatSelect.csv"),
                         delim="\t", col_names=F) %>%
      distinct(.) %>%
      mutate(X2=1) %>%
      rename(!!Name := X2)
    
    USEDSNPS <- USEDSNPS %>%
      left_join(DF_MAT, by=c("Name"="X1"))
  }
}

USEDSNPS %>%
  write_delim("DATA/UsedSNPs_OMNI.txt", delim="\t", col_names=T)

# Process GSA data
FACTORS <- c("Pos", "BAF", "LRRsd", "LRRmean")
MAXDS <- c("10000")

PerfectMatch <- read_delim("DATA/PerfectMatch_GSA_MatSelect", delim="\t", col_names=F) %>%
  mutate(PerfectMatch_0=1)

USEDSNPS <- GSA_MAN %>%
  select(Name) %>%
  mutate(FullSet_0=1) %>%
  left_join(PerfectMatch, by=c("Name"="X1"))

for(i in FACTORS){
  for(j in MAXDS){
    Name <- quo_name(paste0(i, "_", j))
    
    DF_MAT <- read_delim(paste0("DATA/OEE_", i, "_", j, "_RefSelect.csv"),
                         delim="\t", col_names=F) %>%
      distinct(.) %>%
      mutate(X2=1) %>%
      rename(!!Name := X2)
    
    USEDSNPS <- USEDSNPS %>%
      left_join(DF_MAT, by=c("Name"="X1"))
  }
}

USEDSNPS %>%
  write_delim("DATA/UsedSNPs_GSA.txt", delim="\t", col_names=T)

# Process OEE data
FACTORS <- c("Pos", "BAF", "LRRsd", "LRRmean")
MAXDS <- c("10000")

PerfectMatch <- read_delim("DATA/PerfectMatch_OEE_MatSelect", delim="\t", col_names=F) %>%
  mutate(PerfectMatch_0=1)

USEDSNPS <- OEE_MAN %>%
  select(Name) %>%
  mutate(FullSet_0=1) %>%
  left_join(PerfectMatch, by=c("Name"="X1"))

for(i in FACTORS){
  for(j in MAXDS){
    Name <- quo_name(paste0(i, "_", j))
    
    DF_MAT <- read_delim(paste0("DATA/OEE_", i, "_", j, "_MatSelect.csv"),
                         delim="\t", col_names=F) %>%
      distinct(.) %>%
      mutate(X2=1) %>%
      rename(!!Name := X2)
    
    USEDSNPS <- USEDSNPS %>%
      left_join(DF_MAT, by=c("Name"="X1"))
  }
}

USEDSNPS %>%
  write_delim("DATA/UsedSNPs_OEE.txt", delim="\t", col_names=T)
