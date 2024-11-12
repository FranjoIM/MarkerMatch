# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)

# LOAD DATA
load("Validation_Final_103024.RData")

# PREPARE A DATA FILE NAMES FOR SSC
DataFileNames <- data.frame(
  Factor=c("FullSet", "PerfectMatch", rep("BAF", 12), rep("LRRmean", 12), rep("LRRsd", 12), rep("Pos", 12)),
  D_MAX=c(NA, NA, rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4)),
  D_MAXLab=c("0", "0", rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4)),
  stringsAsFactors=FALSE)

# SUMMARIZE CNV CALLSETS
CALLSET_SUMMARIZER <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$D_MAXLab[h] 
  
  for(o in c("Raw", "QCd")){
    
    # PULL DATA FRAMES
    DF_CNV <- DATA[[o]][["SSC"]][[i]][[j]][["CNV"]] %>%
      mutate(CNSize=case_when(
        LEN<100000 ~ "Small",
        LEN>=100000 & LEN<500000 ~ "Medium",
        LEN>=500000 & LEN<1000000 ~ "Large",
        LEN>=1000000 & LEN<5000000 ~ "Very Large",
        LEN>=5000000 ~ "Ultra Large",
        TRUE ~ "Edge Case"),
        CNType=case_when(
          CN < 2 ~ "Deletion",
          CN > 2 ~ "Duplication",
          TRUE ~ "Edge Case"))
    

    for(p in c("All", "Deletions", "Duplications")){
      for(q in c("All", "Small", "Medium", "Large", "Very Large", "Ultra Large")){
        
        if(p=="Deletions"){
          FILT <- DF_CNV %>%
            filter(CN < 2)
        } else if (p=="Duplications") {
          FILT <- DF_CNV %>%
            filter(CN > 2)
        } else if (p=="All") {
          FILT <- DF_CNV
        }
        
        if(q!="All"){
          FILT <- FILT %>%
            filter(CNSize==q)
        }
        
        ROW <- data.frame(NULL)
        
        ROW[1, "Factor"] <- i
        ROW[1, "D_MAX"] <- j
        ROW[1, "QC"] <- o
        ROW[1, "Size"] <- q
        ROW[1, "Type"] <- p
        
        COUNTS <- FILT %>%
          count(ID, name="N")
        
        ROW[1, "N_CNV"] <- nrow(FILT)
        ROW[1, "PER_SAMPLE_N_CNV_MEAN"] <- mean(COUNTS$N)
        ROW[1, "PER_SAMPLE_N_CNV_SD"] <- sd(COUNTS$N)
        ROW[1, "PER_SAMPLE_N_CNV_MEDIAN"] <- median(COUNTS$N)
        ROW[1, "PER_SAMPLE_N_CNV_IQR"] <- IQR(COUNTS$N)
        ROW[1, "PER_SAMPLE_N_CNV_MIN"] <- min(COUNTS$N)
        ROW[1, "PER_SAMPLE_N_CNV_MAX"] <- max(COUNTS$N)
        
        ROW[1, "CNV_SIZE_MEAN"] <- mean(FILT$LEN)
        ROW[1, "CNV_SIZE_SD"] <- sd(FILT$LEN)
        ROW[1, "CNV_SIZE_MEDIAN"] <- median(FILT$LEN)
        ROW[1, "CNV_SIZE_IQR"] <- IQR(FILT$LEN)
        ROW[1, "CNV_SIZE_MIN"] <- min(FILT$LEN)
        ROW[1, "CNV_SIZE_MAX"] <- max(FILT$LEN)
        
        ROW[1, "CNV_CONF_MEAN"] <- mean(FILT$Confidence)
        ROW[1, "CNV_CONF_SD"] <- sd(FILT$Confidence)
        ROW[1, "CNV_CONF_MEDIAN"] <- median(FILT$Confidence)
        ROW[1, "CNV_CONF_IQR"] <- IQR(FILT$Confidence)
        ROW[1, "CNV_CONF_MIN"] <- min(FILT$Confidence)
        ROW[1, "CNV_CONF_MAX"] <- max(FILT$Confidence)
        
        CALLSET_SUMMARIZER <- rbind(CALLSET_SUMMARIZER, ROW) 
        
        rm(ROW)
      }
    }
  }
}

CALLSET_SUMMARIZER %>%
  rename(FactorS=Factor) %>%
  mutate(Factor=case_when(
    FactorS=="FullSet" ~ "Full Set",
    FactorS=="PerfectMatch" ~ "Perfect Match",
    FactorS=="LRRmean" ~ "LRR mean",
    FactorS=="LRRsd" ~ "LRR sd",
    FactorS=="Pos" ~ "Distance",
    TRUE ~ FactorS), .before="FactorS") %>%
  select(-FactorS) %>%
  write_tsv("TABLES/TableS1G.tsv", col_names=TRUE)

# SUMMARIZE SAMPLE CALLSETS
SAMPLE_SUMMARIZER <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$D_MAXLab[h] 
  
  for(o in c("Raw", "QCd")){
    
    # PULL DATA FRAMES
    DF_SAM <- DATA[[o]][["SSC"]][[i]][[j]][["QC"]]
    
    ROW <- data.frame(NULL)
    
    ROW[1, "Factor"] <- i
    ROW[1, "D_MAX"] <- j
    ROW[1, "QC"] <- o

    ROW[1, "N_SAMPLE"] <- nrow(DF_SAM)
    
    ROW[1, "LRR_MEAN_MEAN"] <- mean(DF_SAM$LRR_mean)
    ROW[1, "LRR_MEAN_SD"] <- sd(DF_SAM$LRR_mean)
    
    ROW[1, "LRR_SD_MEAN"] <- mean(DF_SAM$LRR_SD)
    ROW[1, "LRR_SD_SD"] <- sd(DF_SAM$LRR_SD)
    
    ROW[1, "BAF_MEAN_MEAN"] <- mean(DF_SAM$BAF_mean)
    ROW[1, "BAF_MEAN_SD"] <- sd(DF_SAM$BAF_mean)
    
    ROW[1, "BAF_SD_MEAN"] <- mean(DF_SAM$BAF_SD)
    ROW[1, "BAF_SD_SD"] <- sd(DF_SAM$BAF_SD)
    
    ROW[1, "BAF_Drift_MEAN"] <- mean(DF_SAM$BAF_drift)
    ROW[1, "BAF_Drift_SD"] <- sd(DF_SAM$BAF_drift)
    
    ROW[1, "WF_MEAN"] <- mean(DF_SAM$WF)
    ROW[1, "WF_SD"] <- sd(DF_SAM$WF)
    
    ROW[1, "GCWF_MEAN"] <- mean(DF_SAM$GCWF)
    ROW[1, "GCWF_SD"] <- sd(DF_SAM$GCWF)
    
    SAMPLE_SUMMARIZER <- rbind(SAMPLE_SUMMARIZER, ROW) 
    
    rm(ROW)
  }
}

SAMPLE_SUMMARIZER %>%
  rename(FactorS=Factor) %>%
  mutate(Factor=case_when(
    FactorS=="FullSet" ~ "Full Set",
    FactorS=="PerfectMatch" ~ "Perfect Match",
    FactorS=="LRRmean" ~ "LRR mean",
    FactorS=="LRRsd" ~ "LRR sd",
    FactorS=="Pos" ~ "Distance",
    TRUE ~ FactorS), .before="FactorS") %>%
  select(-FactorS) %>%
  write_tsv("TABLES/TableS1H.tsv", col_names=TRUE)

# INITIATE HOLDING DF FOR VALIDATION
ANALYSIS_SSC <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$D_MAXLab[h] 
  
  # Pull out CNVs overlapping in full set and partial sets
  for(o in c("Raw", "QCd")){
    REF <- DATA[["Raw"]][["SSC"]][["FullSet"]][["0"]][["CNV"]]
    MAT <- DATA[[o]][["SSC"]][[i]][[j]][["CNV"]]

    # Keep only same for the analysis
    KEEP_IDs <- intersect(unique(REF$ID), unique(MAT$ID))
    REF <- REF %>%
      filter(ID %in% KEEP_IDs)
    MAT <- MAT %>%
      filter(ID %in% KEEP_IDs)

    # Remove telomeric, centromeric, and immunoglobulin regions from the QCd callset
    if (o=="QCd") {
      MAT <- MAT %>%
        filter(Tel==0 & Cen==0 & Imu==0)}

    # Join the CNVs into ID, CHR, STATE, and CN matched DF
    # Classify CNVs by size and by FP/TP/FN/TN status
    COM <- full_join(REF, MAT, 
                     by=c("ID", "CHR", "STATE", "CN"), 
                     relationship="many-to-many") %>%
      mutate(cIDx=CNV_ID.x,
             cIDy=CNV_ID.y,
             ConfM=case_when(
               is.na(START.y) ~ "FN",
               is.na(START.x) ~ "FP", 
               START.x<=END.y & END.x>=START.y ~ "TP",
               START.x>END.y | END.x<START.y ~ "FP",
               TRUE ~ "Edge Case")) %>%
      mutate(CNSize=case_when(
        is.na(LEN.y) & LEN.x<100000 ~ "Small",
        is.na(LEN.y) & LEN.x>=100000 & LEN.x<500000 ~ "Medium",
        is.na(LEN.y) & LEN.x>=500000 & LEN.x<1000000 ~ "Large",
        is.na(LEN.y) & LEN.x>=1000000 & LEN.x<5000000 ~ "Very Large",
        is.na(LEN.y) & LEN.x>=5000000 ~ "Ultra Large",
        !is.na(LEN.y) & LEN.y<100000 ~ "Small",
        !is.na(LEN.y) & LEN.y>=100000 & LEN.y<500000 ~ "Medium",
        !is.na(LEN.y) & LEN.y>=500000 & LEN.y<1000000 ~ "Large",
        !is.na(LEN.y) & LEN.y>=1000000 & LEN.y<5000000 ~ "Very Large",
        !is.na(LEN.y) & LEN.y>=5000000 ~ "Ultra Large",
        TRUE ~ "Edge Case"))
    
    COM_TEMP1 <- COM %>%
      filter(ConfM == "TP")
    
    COM_TP <- COM_TEMP1 %>%
      pull(CNV_ID.y) %>%
      unique(.)
    
    COM_TEMP2 <- COM %>%
      filter(!CNV_ID.y %in% COM_TP) %>%
      filter(ConfM == "FN")
    
    COM_FN <- COM_TEMP2 %>%
      pull(CNV_ID.x) %>%
      unique(.)
    
    COM_TEMP3 <- COM %>%
      filter(!CNV_ID.y %in% COM_TP) %>%
      filter(!CNV_ID.x %in% COM_FN)
    
    COM_FP <- COM_TEMP3 %>%
      pull(CNV_ID.y) %>%
      unique(.)
    
    COM_TEMP4 <- COM %>%
      filter(!CNV_ID.y %in% COM_TP) %>%
      filter(!CNV_ID.x %in% COM_FN) %>%
      filter(!CNV_ID.y %in% COM_FP)
    
    if(nrow(COM_TEMP4) > 0) {message("Potential problem in confusion matrix.\n")}
    
    COM <- rbind(COM_TEMP1, COM_TEMP2, COM_TEMP3, COM_TEMP4)
    
    for(p in c("All", "Deletions", "Duplications")){
      for(q in c("All", "Small", "Medium", "Large", "Very Large", "Ultra Large")){
        
        if(p=="Deletions"){
          FILT <- COM %>%
            filter(CN < 2)
        } else if (p=="Duplications") {
          FILT <- COM %>%
            filter(CN > 2)
        } else if (p=="All") {
          FILT <- COM
        }
        
        if(q!="All"){
          FILT <- FILT %>%
            filter(CNSize==q)
        }
        
        TP <- FILT %>%
          filter(ConfM=="TP") %>%
          pull(cIDy) %>%
          unique(.)
        
        FP <- FILT %>%
          filter(ConfM=="FP") %>%
          filter(!cIDy %in% TP) %>%
          pull(cIDy) %>%
          unique(.)
        
        FN <- FILT %>%
          filter(ConfM=="FN") %>%
          pull(cIDx) %>%
          unique(.)
        
        NEW_ROW <- data.frame(
          Matching_Method=i,
          Matching_Distance=j,
          Matching_Type=o,
          CNV_Type=p,
          CNV_Size=q,
          TP=length(TP),
          FP=length(FP),
          FN=length(FN))
        
        ANALYSIS_SSC <- rbind(ANALYSIS_SSC, NEW_ROW)
        
        timestamp()
        message("Finished ", i, "-", j, "-", o, "-", p, "-", q, ".\n")
        
      }
    }
  }
}

# CHECKPOINT
save(ANALYSIS_SSC, file="Analysis_SSC_103024.RData")
load("Analysis_SSC_103024.RData")

# TIDY UP ANALYSIS
ANALYSIS_SSC <- ANALYSIS_SSC %>%
  mutate(Sensitivity=round(TP/(TP+FN), digits=3),
         PPV=round(TP/(FP+TP), digits=3),
         FNR=round(FN/(FN+TP), digits=3),
         FDR=round(FP/(FP+TP), digits=3),
         F1=round((2*TP)/(2*TP+FP+FN), digits=3),
         FMI=round(sqrt((TP/(FP+TP))*(TP/(TP+FN))), digits=3),
         JI=round(TP/(TP+FN+FP), digits=3)) %>%
  mutate(Matching_Distance=as.numeric(Matching_Distance)) %>%
  mutate(D_MAX_LOG=log10(Matching_Distance))

# DEFINE PLOTTING FUNCTION
MetricPlot <- function(a, b, c){
  H1 <- ANALYSIS_SSC %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Type==a & CNV_Size==b & Matching_Type=="Raw") %>%
    pull(.data[[c]])
  
  H2 <- ANALYSIS_SSC %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Type==a & CNV_Size==b & Matching_Type=="QCd") %>%
    pull(.data[[c]])
  
  H3 <- ANALYSIS_SSC %>%
    filter(Matching_Method=="FullSet" & CNV_Type==a & CNV_Size==b & Matching_Type=="Raw") %>%
    pull(.data[[c]])
  
  H4 <- ANALYSIS_SSC %>%
    filter(Matching_Method=="FullSet" & CNV_Type==a & CNV_Size==b & Matching_Type=="QCd") %>%
    pull(.data[[c]])
  
  TITLE <- case_when(
    b == "All" ~ "All CNVs",
    b == "Small" ~ "CNV < 100kb",
    b == "Medium" ~ "100kb < CNV < 500kb",
    b == "Large" ~ "500kb < CNV < 1Mb",
    b == "Very Large" ~ "1Mb < CNV < 5Mb",
    b == "Ultra Large" ~ "5Mb < CNV",
    TRUE ~ NA_character_)
  
  PLOT <- ANALYSIS_SSC %>%
    filter(!Matching_Method %in% c("PerfectMatch", "FullSet")) %>%
    filter(CNV_Type==a & CNV_Size==b) %>%
    ggplot(aes(x=D_MAX_LOG, y=.data[[c]], linetype=Matching_Type, color=Matching_Method)) +
    geom_hline(aes(yintercept=H1, color="Perfect Match", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H2, color="Perfect Match", linetype="QCd"), linewidth=1) +
    geom_hline(aes(yintercept=H3, color="Reference", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H4, color="Reference", linetype="QCd"), linewidth=1) +
    geom_line(linewidth=1) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRRmean", "LRRsd", "Pos", "Perfect Match", "Reference")) +
    labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
         y=toupper(c),
         linetype="CNV CALLSET QC",
         color="MATCHING METHOD",
         subtitle=TITLE) +
    ylim(0, 1) +
    theme_bw() + 
    theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold")) +
    guides(color=guide_legend(title.position="top", nrow=1),
           linetype=guide_legend(title.position="top"))
}

# PREPARE FOR PLOTTING
A <- c("All", "Deletions", "Duplications")
A <- set_names(A)

B <- c("All", "Small", "Medium", "Large", "Very Large", "Ultra Large")
B <- set_names(B)

C <- c("Sensitivity", "PPV", "FNR", "FDR", "F1", "FMI", "JI")
C <- set_names(C)

# PLOT METRICS
PLOTS <- map(A, function(x) map(B, function(y) map(C, function(z) MetricPlot(a=x, b=y, c=z))))
                                            
# SENSITIVITY, ALL (PLOT 8)
ggarrange(PLOTS$All$All$Sensitivity + rremove("xlab"), 
          PLOTS$All$Small$Sensitivity + rremove("xlab") + rremove("ylab"),
          PLOTS$All$Medium$Sensitivity + rremove("xlab") + rremove("ylab"), 
          PLOTS$All$Large$Sensitivity,
          PLOTS$All$`Very Large`$Sensitivity + rremove("ylab"),
          PLOTS$All$`Ultra Large`$Sensitivity + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot8.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")
                                                   
# SENSITIVITY, DELETIONS (PLOT 9)
ggarrange(PLOTS$Deletions$All$Sensitivity + rremove("xlab"), 
          PLOTS$Deletions$Small$Sensitivity + rremove("xlab") + rremove("ylab"),
          PLOTS$Deletions$Medium$Sensitivity + rremove("xlab") + rremove("ylab"), 
          PLOTS$Deletions$Large$Sensitivity,
          PLOTS$Deletions$`Very Large`$Sensitivity + rremove("ylab"),
          PLOTS$Deletions$`Ultra Large`$Sensitivity + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot9.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# SENSITIVITY, DUPLICATIONS (PLOT 10)
ggarrange(PLOTS$Duplications$All$Sensitivity + rremove("xlab"), 
          PLOTS$Duplications$Small$Sensitivity + rremove("xlab") + rremove("ylab"),
          PLOTS$Duplications$Medium$Sensitivity + rremove("xlab") + rremove("ylab"), 
          PLOTS$Duplications$Large$Sensitivity,
          PLOTS$Duplications$`Very Large`$Sensitivity + rremove("ylab"),
          PLOTS$Duplications$`Ultra Large`$Sensitivity + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot10.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")
                                                   
# PPV, ALL (PLOT 11)
ggarrange(PLOTS$All$All$PPV + rremove("xlab") + ylim(0.75, 1), 
          PLOTS$All$Small$PPV + rremove("xlab") + rremove("ylab") + ylim(0.75, 1),
          PLOTS$All$Medium$PPV + rremove("xlab") + rremove("ylab") + ylim(0.75, 1), 
          PLOTS$All$Large$PPV + ylim(0.75, 1),
          PLOTS$All$`Very Large`$PPV + rremove("ylab") + ylim(0.75, 1),
          PLOTS$All$`Ultra Large`$PPV + rremove("ylab") + ylim(0.75, 1),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top")  %>%
  ggsave(filename="FIGURES/Plot11.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# PPV, DELETIONS (PLOT 12)
ggarrange(PLOTS$Deletions$All$PPV + rremove("xlab") + ylim(0.75, 1), 
          PLOTS$Deletions$Small$PPV + rremove("xlab") + rremove("ylab") + ylim(0.75, 1),
          PLOTS$Deletions$Medium$PPV + rremove("xlab") + rremove("ylab") + ylim(0.75, 1), 
          PLOTS$Deletions$Large$PPV + ylim(0.75, 1),
          PLOTS$Deletions$`Very Large`$PPV + rremove("ylab") + ylim(0.75, 1),
          PLOTS$Deletions$`Ultra Large`$PPV + rremove("ylab") + ylim(0.75, 1),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot12.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# PPV, DUPLICATIONS (PLOT 13)
ggarrange(PLOTS$Duplications$All$PPV + rremove("xlab") + ylim(0.7, 1), 
          PLOTS$Duplications$Small$PPV + rremove("xlab") + rremove("ylab") + ylim(0.7, 1),
          PLOTS$Duplications$Medium$PPV + rremove("xlab") + rremove("ylab") + ylim(0.7, 1), 
          PLOTS$Duplications$Large$PPV + ylim(0.7, 1),
          PLOTS$Duplications$`Very Large`$PPV + rremove("ylab") + ylim(0.7, 1),
          PLOTS$Duplications$`Ultra Large`$PPV + rremove("ylab") + ylim(0.7, 1),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot13.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# FNR, ALL (PLOT 14)
ggarrange(PLOTS$All$All$FNR + rremove("xlab"), 
          PLOTS$All$Small$FNR + rremove("xlab") + rremove("ylab"),
          PLOTS$All$Medium$FNR + rremove("xlab") + rremove("ylab"), 
          PLOTS$All$Large$FNR,
          PLOTS$All$`Very Large`$FNR + rremove("ylab"),
          PLOTS$All$`Ultra Large`$FNR + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot14.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# FNR, DELETIONS (PLOT 15)
ggarrange(PLOTS$Deletions$All$FNR + rremove("xlab"), 
          PLOTS$Deletions$Small$FNR + rremove("xlab") + rremove("ylab"),
          PLOTS$Deletions$Medium$FNR + rremove("xlab") + rremove("ylab"), 
          PLOTS$Deletions$Large$FNR,
          PLOTS$Deletions$`Very Large`$FNR + rremove("ylab"),
          PLOTS$Deletions$`Ultra Large`$FNR + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot15.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")                                               
                                                   
# FNR, DUPLICATIONS (PLOT 16)
ggarrange(PLOTS$Duplications$All$FNR + rremove("xlab"), 
          PLOTS$Duplications$Small$FNR + rremove("xlab") + rremove("ylab"),
          PLOTS$Duplications$Medium$FNR + rremove("xlab") + rremove("ylab"), 
          PLOTS$Duplications$Large$FNR,
          PLOTS$Duplications$`Very Large`$FNR + rremove("ylab"),
          PLOTS$Duplications$`Ultra Large`$FNR + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot16.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")    

# FDR, ALL (PLOT 17)
ggarrange(PLOTS$All$All$FDR + rremove("xlab") + ylim(0, 0.25), 
          PLOTS$All$Small$FDR + rremove("xlab") + rremove("ylab") + ylim(0, 0.25),
          PLOTS$All$Medium$FDR + rremove("xlab") + rremove("ylab") + ylim(0, 0.25), 
          PLOTS$All$Large$FDR + ylim(0, 0.25),
          PLOTS$All$`Very Large`$FDR + rremove("ylab") + ylim(0, 0.25),
          PLOTS$All$`Ultra Large`$FDR + rremove("ylab") + ylim(0, 0.25),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot17.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white") 

# FDR, DELETIONS (PLOT 18)
ggarrange(PLOTS$Deletions$All$FDR + rremove("xlab") + ylim(0, 0.25), 
          PLOTS$Deletions$Small$FDR + rremove("xlab") + rremove("ylab") + ylim(0, 0.25),
          PLOTS$Deletions$Medium$FDR + rremove("xlab") + rremove("ylab") + ylim(0, 0.25), 
          PLOTS$Deletions$Large$FDR + ylim(0, 0.25),
          PLOTS$Deletions$`Very Large`$FDR + rremove("ylab") + ylim(0, 0.25),
          PLOTS$Deletions$`Ultra Large`$FDR + rremove("ylab") + ylim(0, 0.25),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot18.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white") 

# FDR, DUPLICATIONS (PLOT 19)
ggarrange(PLOTS$Duplications$All$FDR + rremove("xlab") + ylim(0, 0.3), 
          PLOTS$Duplications$Small$FDR + rremove("xlab") + rremove("ylab") + ylim(0, 0.3),
          PLOTS$Duplications$Medium$FDR + rremove("xlab") + rremove("ylab") + ylim(0, 0.3), 
          PLOTS$Duplications$Large$FDR + ylim(0, 0.3),
          PLOTS$Duplications$`Very Large`$FDR + rremove("ylab") + ylim(0, 0.3),
          PLOTS$Duplications$`Ultra Large`$FDR + rremove("ylab") + ylim(0, 0.3),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot19.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white") 

# F1, ALL (PLOT 20)
ggarrange(PLOTS$All$All$F1 + rremove("xlab"), 
          PLOTS$All$Small$F1 + rremove("xlab") + rremove("ylab"),
          PLOTS$All$Medium$F1 + rremove("xlab") + rremove("ylab"), 
          PLOTS$All$Large$F1,
          PLOTS$All$`Very Large`$F1 + rremove("ylab"),
          PLOTS$All$`Ultra Large`$F1 + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot20.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# F1, DELETIONS (PLOT 21)
ggarrange(PLOTS$Deletions$All$F1 + rremove("xlab"), 
          PLOTS$Deletions$Small$F1 + rremove("xlab") + rremove("ylab"),
          PLOTS$Deletions$Medium$F1 + rremove("xlab") + rremove("ylab"), 
          PLOTS$Deletions$Large$F1,
          PLOTS$Deletions$`Very Large`$F1 + rremove("ylab"),
          PLOTS$Deletions$`Ultra Large`$F1 + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot21.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# F1, DUPLICATIONS (PLOT 22)
ggarrange(PLOTS$Duplications$All$F1 + rremove("xlab"), 
          PLOTS$Duplications$Small$F1 + rremove("xlab") + rremove("ylab"),
          PLOTS$Duplications$Medium$F1 + rremove("xlab") + rremove("ylab"), 
          PLOTS$Duplications$Large$F1,
          PLOTS$Duplications$`Very Large`$F1 + rremove("ylab"),
          PLOTS$Duplications$`Ultra Large`$F1 + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot22.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FMI, ALL (PLOT 23)
ggarrange(PLOTS$All$All$FMI + rremove("xlab"), 
          PLOTS$All$Small$FMI + rremove("xlab") + rremove("ylab"),
          PLOTS$All$Medium$FMI + rremove("xlab") + rremove("ylab"), 
          PLOTS$All$Large$FMI,
          PLOTS$All$`Very Large`$FMI + rremove("ylab"),
          PLOTS$All$`Ultra Large`$FMI + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot23.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FMI, DELETIONS (PLOT 24)
ggarrange(PLOTS$Deletions$All$FMI + rremove("xlab"), 
          PLOTS$Deletions$Small$FMI + rremove("xlab") + rremove("ylab"),
          PLOTS$Deletions$Medium$FMI + rremove("xlab") + rremove("ylab"), 
          PLOTS$Deletions$Large$FMI,
          PLOTS$Deletions$`Very Large`$FMI + rremove("ylab"),
          PLOTS$Deletions$`Ultra Large`$FMI + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot24.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FMI, DUPLICATIONS (PLOT 25)
ggarrange(PLOTS$Duplications$All$FMI + rremove("xlab"), 
          PLOTS$Duplications$Small$FMI + rremove("xlab") + rremove("ylab"),
          PLOTS$Duplications$Medium$FMI + rremove("xlab") + rremove("ylab"), 
          PLOTS$Duplications$Large$FMI,
          PLOTS$Duplications$`Very Large`$FMI + rremove("ylab"),
          PLOTS$Duplications$`Ultra Large`$FMI + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot25.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# JI, ALL (PLOT 26)
ggarrange(PLOTS$All$All$JI + rremove("xlab"), 
          PLOTS$All$Small$JI + rremove("xlab") + rremove("ylab"),
          PLOTS$All$Medium$JI + rremove("xlab") + rremove("ylab"), 
          PLOTS$All$Large$JI,
          PLOTS$All$`Very Large`$JI + rremove("ylab"),
          PLOTS$All$`Ultra Large`$JI + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot26.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# JI, DELETIONS (PLOT 27)
ggarrange(PLOTS$Deletions$All$JI + rremove("xlab"), 
          PLOTS$Deletions$Small$JI + rremove("xlab") + rremove("ylab"),
          PLOTS$Deletions$Medium$JI + rremove("xlab") + rremove("ylab"), 
          PLOTS$Deletions$Large$JI,
          PLOTS$Deletions$`Very Large`$JI + rremove("ylab"),
          PLOTS$Deletions$`Ultra Large`$JI + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot27.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# JI, DUPLICATIONS (PLOT 28)
ggarrange(PLOTS$Duplications$All$JI + rremove("xlab"), 
          PLOTS$Duplications$Small$JI + rremove("xlab") + rremove("ylab"),
          PLOTS$Duplications$Medium$JI + rremove("xlab") + rremove("ylab"), 
          PLOTS$Duplications$Large$JI,
          PLOTS$Duplications$`Very Large`$JI + rremove("ylab"),
          PLOTS$Duplications$`Ultra Large`$JI + rremove("ylab"),
          align="hv", labels=c("A", "B", "C", "D", "E", "F"), common.legend=T,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot28.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# WRITE TABLE DATA
ANALYSIS_SSC %>%
  rename(FactorS=Matching_Method,
    D_MAX=Matching_Distance,
    QC=Matching_Type) %>%
  select(-D_MAX_LOG) %>%
  mutate(Factor=case_when(
    FactorS=="FullSet" ~ "Full Set",
    FactorS=="PerfectMatch" ~ "Perfect Match",
    FactorS=="LRRmean" ~ "LRR mean",
    FactorS=="LRRsd" ~ "LRR sd",
    FactorS=="Pos" ~ "Distance",
    TRUE ~ FactorS), .before="FactorS") %>%
  select(-FactorS) %>%
  write_tsv("TABLES/TableS1I.tsv", col_names=TRUE)
