#!/usr/bin/env Rscript

# Import tidyverse
library(tidyverse)

# Define working directory
WK_DIR <- "..."

# Import OMNI PFB file and selection file
PFB_OEE <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/OEE_FullSet.pfb"), delim="\t")
PFB_GSA <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/GSA1_FullSet.pfb"), delim="\t")

# Create Pos 10k OEE
Pos_10000_OEE <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_Pos_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB_OEE %>%
    filter(Name %in% Pos_10000_OEE$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_OEE_Pos_10000.pfb"), delim="\t")

# Create Pos 10k GSA
Pos_10000_GSA <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_Pos_10000_RefSelect.csv"), delim="\t", col_names=F)
PFB_GSA %>%
    filter(Name %in% Pos_10000_GSA$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_GSA_Pos_10000.pfb"), delim="\t")

# Create BAF 10k OEE
BAF_10000_OEE <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_BAF_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB_OEE %>%
    filter(Name %in% BAF_10000_OEE$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_OEE_BAF_10000.pfb"), delim="\t")

# Create BAF 10k GSA
BAF_10000_GSA <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_BAF_10000_RefSelect.csv"), delim="\t", col_names=F)
PFB_GSA %>%
    filter(Name %in% BAF_10000_GSA$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_GSA_BAF_10000.pfb"), delim="\t")

# Create LRRmean 10k OEE
LRRmean_10000_OEE <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_LRRmean_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB_OEE %>%
    filter(Name %in% LRRmean_10000_OEE$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_OEE_LRRmean_10000.pfb"), delim="\t")

# Create Perfect Match OEE
PerfectMatch_OEE <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/PerfectMatch_OEE_MatSelect"), delim="\t", col_names=F)
PFB_OEE %>%
    filter(Name %in% PerfectMatch_OEE$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_OEE_PerfectMatch_0.pfb"), delim="\t")

# Create LRRmean 10k GSA
LRRmean_10000_GSA <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_LRRmean_10000_RefSelect.csv"), delim="\t", col_names=F)
PFB_GSA %>%
    filter(Name %in% LRRmean_10000_GSA$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_GSA_LRRmean_10000.pfb"), delim="\t")

# Create LRRsd 10k OEE
LRRsd_10000_OEE <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_LRRsd_10000_MatSelect.csv"), delim="\t", col_names=F)
PFB_OEE %>%
    filter(Name %in% LRRsd_10000_OEE$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_OEE_LRRsd_10000.pfb"), delim="\t")

# Create LRRsd 10k GSA
LRRsd_10000_GSA <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/OEE_LRRsd_10000_RefSelect.csv"), delim="\t", col_names=F)
PFB_GSA %>%
    filter(Name %in% LRRsd_10000_GSA$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_GSA_LRRsd_10000.pfb"), delim="\t")

# Create LRRsd Perfect Match GSA
PerfectMatch_GSA <- readr::read_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/PerfectMatch_GSA_MatSelect"), delim="\t", col_names=F)
PFB_GSA %>%
    filter(Name %in% PerfectMatch_GSA$X1) %>%
    readr::write_delim(paste0(WK_DIR, "/SupportingFiles/MarkerMatch/Validation2_GSA_PerfectMatch_0.pfb"), delim="\t")
