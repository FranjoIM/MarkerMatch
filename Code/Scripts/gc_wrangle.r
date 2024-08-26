#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)

# Update SNP Pos files with GC content percentage and export
snppos <- read.table(args[1], header=TRUE)
gcrsid <- read.table(args[2], header=TRUE)

df <- snppos %>%
    left_join(gcrsid, by="Name") %>%
    filter(!Chr %in% c("MT", "X", "Y", "XY", "0"))

write.table(df, args[3], col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", na="")
