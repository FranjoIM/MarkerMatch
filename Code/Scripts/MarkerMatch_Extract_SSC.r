#!/usr/bin/env Rscript

# Import tidyverse
library(tidyverse)

# Define working directory
WK_DIR <- "..."

# Import OMNI PFB file and selection file
PFB <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/Omni_FullSet.pfb"), delim="\t")

# Crete Pos10 PFB
Pos_10 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_10_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_10$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_10.pfb"), delim="\t")

# Crete BAF10 PFB
BAF_10 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_10_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_10$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_10.pfb"), delim="\t")

# Crete LRRsd10 PFB
LRRsd_10 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_10_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_10$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_10.pfb"), delim="\t")

# Crete LRRmean10 PFB
LRRmean_10 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_10_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_10$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_10.pfb"), delim="\t")

# Crete Pos50 PFB
Pos_50 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_50_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_50$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_50.pfb"), delim="\t")

# Crete BAF50 PFB
BAF_50 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_50_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_50$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_50.pfb"), delim="\t")

# Crete LRRsd50 PFB
LRRsd_50 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_50_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_50$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_50.pfb"), delim="\t")

# Crete LRRmean50 PFB
LRRmean_50 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_50_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_50$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_50.pfb"), delim="\t")

# Crete Pos100 PFB
Pos_100 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_100_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_100$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_100.pfb"), delim="\t")

# Crete BAF100 PFB
BAF_100 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_100_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_100$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_100.pfb"), delim="\t")

# Crete LRRsd100 PFB
LRRsd_100 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_100_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_100$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_100.pfb"), delim="\t")

# Crete LRRmean100 PFB
LRRmean_100 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_100_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_100$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_100.pfb"), delim="\t")

# Crete Pos500 PFB
Pos_500 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_500_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_500$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_500.pfb"), delim="\t")

# Crete BAF500 PFB
BAF_500 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_500_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_500$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_500.pfb"), delim="\t")

# Crete LRRsd500 PFB
LRRsd_500 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_500_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_500$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_500.pfb"), delim="\t")

# Crete LRRmean500 PFB
LRRmean_500 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_500_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_500$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_500.pfb"), delim="\t")

# Crete Pos1000 PFB
Pos_1000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_1000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_1000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_1000.pfb"), delim="\t")

# Crete BAF1000 PFB
BAF_1000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_1000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_1000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_1000.pfb"), delim="\t")

# Crete LRRsd1000 PFB
LRRsd_1000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_1000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_1000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_1000.pfb"), delim="\t")

# Crete LRRmean1000 PFB
LRRmean_1000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_1000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_1000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_1000.pfb"), delim="\t")

# Crete Pos5000 PFB
Pos_5000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_5000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_5000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_5000.pfb"), delim="\t")

# Crete BAF5000 PFB
BAF_5000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_5000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_5000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_5000.pfb"), delim="\t")

# Crete LRRsd5000 PFB
LRRsd_5000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_5000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_5000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_5000.pfb"), delim="\t")

# Crete LRRmean5000 PFB
LRRmean_5000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_5000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_5000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_5000.pfb"), delim="\t")

# Crete Pos10000 PFB
Pos_10000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_10000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_10000.pfb"), delim="\t")

# Crete BAF10000 PFB
BAF_10000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_10000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_10000.pfb"), delim="\t")

# Crete LRRsd10000 PFB
LRRsd_10000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_10000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_10000.pfb"), delim="\t")

# Crete LRRmean10000 PFB
LRRmean_10000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_10000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_10000.pfb"), delim="\t")

# Crete Pos50000 PFB
Pos_50000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_50000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_50000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_50000.pfb"), delim="\t")

# Crete BAF50000 PFB
BAF_50000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_50000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_50000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_50000.pfb"), delim="\t")

# Crete LRRsd50000 PFB
LRRsd_50000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_50000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_50000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_50000.pfb"), delim="\t")

# Crete LRRmean50000 PFB
LRRmean_50000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_50000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_50000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_50000.pfb"), delim="\t")

# Crete Pos100000 PFB
Pos_100000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_100000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_100000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_100000.pfb"), delim="\t")

# Crete BAF100000 PFB
BAF_100000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_100000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_100000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_100000.pfb"), delim="\t")

# Crete LRRsd100000 PFB
LRRsd_100000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_100000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_100000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_100000.pfb"), delim="\t")

# Crete LRRmean100000 PFB
LRRmean_100000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_100000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_100000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_100000.pfb"), delim="\t")

# Crete Pos500000 PFB
Pos_500000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_500000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_500000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_500000.pfb"), delim="\t")

# Crete BAF500000 PFB
BAF_500000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_500000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_500000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_500000.pfb"), delim="\t")

# Crete LRRsd500000 PFB
LRRsd_500000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_500000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_500000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_500000.pfb"), delim="\t")

# Crete LRRmean500000 PFB
LRRmean_500000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_500000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_500000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_500000.pfb"), delim="\t")

# Crete Pos1000000 PFB
Pos_1000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_1000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_1000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_1000000.pfb"), delim="\t")

# Crete BAF1000000 PFB
BAF_1000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_1000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_1000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_1000000.pfb"), delim="\t")

# Crete LRRsd1000000 PFB
LRRsd_1000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_1000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_1000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_1000000.pfb"), delim="\t")

# Crete LRRmean1000000 PFB
LRRmean_1000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_1000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_1000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_1000000.pfb"), delim="\t")

# Crete Pos5000000 PFB
Pos_5000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_5000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% Pos_5000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Pos_5000000.pfb"), delim="\t")

# Crete BAF5000000 PFB
BAF_5000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_5000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% BAF_5000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/BAF_5000000.pfb"), delim="\t")

# Crete LRRsd5000000 PFB
LRRsd_5000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_5000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRsd_5000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRsd_5000000.pfb"), delim="\t")

# Crete LRRmean5000000 PFB
LRRmean_5000000 <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_5000000_MatSelect.csv"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% LRRmean_5000000$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/LRRmean_5000000.pfb"), delim="\t")

# Crete PerfectMatch PFB
PerfectMatch <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/PerfectMatch_OMNI_MatSelect"), delim="\t", col_names=F)
PFB %>% 
    filter(Name %in% PerfectMatch$X1)%>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/PerfectMatch.pfb"), delim="\t")
