# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)

# LOAD DATA
load("Validation_Final_082924.RData")

# PREPARE A DATA FILE NAMES FOR SSC
DataFileNames <- data.frame(
  Factor=c("FullSet", "PerfectMatch", rep("BAF", 12), rep("LRRmean", 12), rep("LRRsd", 12), rep("Pos", 12)),
  MaxD=c(NA, NA, rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4)),
  MaxDLab=c("0", "0", rep(c("10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000", "500000", "1000000", "5000000"), 4)),
  stringsAsFactors=FALSE)

# INITIATE HOLDING DF
ANALYSIS_REGIONAL_SSC <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$MaxDLab[h] 
  
  # Pull out CNVs overlapping in full set and partial sets
  for(o in c("Raw", "QCd")){
    REF <- DATA[["Raw"]][["SSC"]][["FullSet"]][["0"]][["CNV"]] %>%
      mutate(N_SNP = as.numeric(NumSNP),
             LEN = as.numeric(gsub(",", "", Length)),
             CHR = as.numeric(CHR),
             START = as.numeric(START),
             END = as.numeric(END),
             STATE = as.numeric(STATE),
             CN = as.numeric(CN)) %>%
      select(-c(StartSNP, EndSNP, NumSNP, Length))
    MAT <- DATA[[o]][["SSC"]][[i]][[j]][["CNV"]] %>%
      mutate(N_SNP = as.numeric(NumSNP),
             LEN = as.numeric(gsub(",", "", Length)),
             CHR = as.numeric(CHR),
             START = as.numeric(START),
             END = as.numeric(END),
             STATE = as.numeric(STATE),
             CN = as.numeric(CN)) %>%
      select(-c(StartSNP, EndSNP, NumSNP, Length))
    
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
      filter(!CNV_ID.x %in% COM_FN) %>%
      filter(ConfM == "FP")
    
    COM_FP <- COM_TEMP3 %>%
      pull(CNV_ID.y) %>%
      unique(.)
    
    COM_TEMP4 <- COM %>%
      filter(!CNV_ID.y %in% COM_TP) %>%
      filter(!CNV_ID.x %in% COM_FN) %>%
      filter(!CNV_ID.y %in% COM_FP)
    
    if(nrow(COM_TEMP4) > 0) {message("Potential problem in confusion matrix.\n")}
    
    COM <- rbind(COM_TEMP1, COM_TEMP2, COM_TEMP3, COM_TEMP4)
    
    for(p in c("Tel", "Cen", "SegDup")){
        
      if(p=="Tel"){
        p_Lab <- "Telomeric"
        p_Regx <- "Tel.x"
        p_Regy <- "Tel.y"
      } else if (p=="Cen") {
        p_Lab <- "Centromeric"
        p_Regx <- "Cen.x"
        p_Regy <- "Cen.y"
      } else if (p=="SegDup") {
        p_Lab <- "Segmental Duplication"
        p_Regx <- "SegDup.x"
        p_Regy <- "SegDup.y"
      }
      
      TP <- COM %>%
        filter(.[[p_Regy]]==1) %>%
        filter(ConfM=="TP") %>%
        pull(cIDy) %>%
        unique(.)
      
      FP <- COM %>%
        filter(.[[p_Regy]]==1) %>%
        filter(ConfM=="FP") %>%
        filter(!cIDy %in% TP) %>%
        pull(cIDy) %>%
        unique(.)
      
      FN <- COM %>%
        filter(.[[p_Regx]]==1) %>%
        filter(ConfM=="FN") %>%
        pull(cIDx) %>%
        unique(.)
      
      NEW_ROW <- data.frame(
        Matching_Method=i,
        Matching_Distance=j,
        Matching_Type=o,
        CNV_Region=p,
        TP=length(TP),
        FP=length(FP),
        FN=length(FN))
      
      ANALYSIS_REGIONAL_SSC <- rbind(ANALYSIS_REGIONAL_SSC, NEW_ROW)
      
      timestamp()
      message("Finished ", i, "-", j, "-", o, "-", p, ".\n")
    }
  }
}

# TIDY UP ANALYSIS
ANALYSIS_REGIONAL_SSC <- ANALYSIS_REGIONAL_SSC %>%
  mutate(Sensitivity=round(TP/(TP+FN), digits=3),
         PPV=round(TP/(FP+TP), digits=3),
         FNR=round(FN/(FN+TP), digits=3),
         FDR=round(FP/(FP+TP), digits=3),
         F1=round((2*TP)/(2*TP+FP+FN), digits=3),
         FMI=round(sqrt((TP/(FP+TP))*(TP/(TP+FN))), digits=3),
         JI=round(TP/(TP+FN+FP), digits=3)) %>%
  mutate(Matching_Distance=as.numeric(Matching_Distance)) %>%
  mutate(MaxD_LOG=log10(Matching_Distance))

# DEFINE PLOTTING FUNCTION
MetricPlot <- function(a, b){
  H1 <- ANALYSIS_REGIONAL_SSC %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Region==a & Matching_Type=="Raw") %>%
    pull(.data[[b]])
  
  H2 <- ANALYSIS_REGIONAL_SSC %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Region==a & Matching_Type=="QCd") %>%
    pull(.data[[b]])
  
  H3 <- ANALYSIS_REGIONAL_SSC %>%
    filter(Matching_Method=="FullSet" & CNV_Region==a & Matching_Type=="Raw") %>%
    pull(.data[[b]])
  
  H4 <- ANALYSIS_REGIONAL_SSC %>%
    filter(Matching_Method=="FullSet" & CNV_Region==a & Matching_Type=="QCd") %>%
    pull(.data[[b]])
  
  TITLE <- case_when(
    a == "Tel" ~ "Telomeric",
    a == "Cen" ~ "Centromeric",
    a == "SegDup" ~ "Segmental Duplications",
    TRUE ~ NA_character_)
  
  PLOT <- ANALYSIS_REGIONAL_SSC %>%
    filter(!Matching_Method %in% c("PerfectMatch", "FullSet")) %>%
    filter(CNV_Region==a) %>%
    ggplot(aes(x=MaxD_LOG, y=.data[[b]], linetype=Matching_Type, color=Matching_Method)) +
    geom_hline(aes(yintercept=H1, color="Perfect Match", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H2, color="Perfect Match", linetype="QCd"), linewidth=1) +
    geom_hline(aes(yintercept=H3, color="Reference", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H4, color="Reference", linetype="QCd"), linewidth=1) +
    geom_line(linewidth=1) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRRmean", "LRRsd", "Pos", "Perfect Match", "Reference")) +
    labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
         y=toupper(b),
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
A <- c("Tel", "Cen", "SegDup")
A <- set_names(A)

B <- c("Sensitivity", "PPV", "FNR", "FDR", "F1", "FMI", "JI")
B <- set_names(B)

# PLOT METRICS
PLOTS <- map(A, function(x) map(B, function(y) MetricPlot(a=x, b=y)))

# SENSITIVITY, ALL (PLOT 29)
ggarrange(PLOTS$Tel$Sensitivity, 
          PLOTS$Cen$Sensitivity + rremove("ylab"),
          PLOTS$SegDup$Sensitivity + rremove("ylab"), 
          align="hv", labels=c("A", "B", "C"), common.legend=T, nrow=1,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot29.png",
         device="png",
         width=11,
         height=4,
         units="in",
         dpi=350,
         bg="white")

# PPV, ALL (PLOT 30)
ggarrange(PLOTS$Tel$PPV + ylim(0.65, 1), 
          PLOTS$Cen$PPV + rremove("ylab") + ylim(0.65, 1),
          PLOTS$SegDup$PPV + rremove("ylab") + ylim(0.65, 1), 
          align="hv", labels=c("A", "B", "C"), common.legend=T, nrow=1,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot30.png",
         device="png",
         width=11,
         height=4,
         units="in",
         dpi=350,
         bg="white")

# FNR, ALL (PLOT 31)
ggarrange(PLOTS$Tel$FNR, 
          PLOTS$Cen$FNR + rremove("ylab"),
          PLOTS$SegDup$FNR + rremove("ylab"), 
          align="hv", labels=c("A", "B", "C"), common.legend=T, nrow=1,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot31.png",
         device="png",
         width=11,
         height=4,
         units="in",
         dpi=350,
         bg="white")

# FDR, ALL (PLOT 32)
ggarrange(PLOTS$Tel$FDR + ylim(0, 0.35), 
          PLOTS$Cen$FDR + rremove("ylab") + ylim(0, 0.35),
          PLOTS$SegDup$FDR + rremove("ylab") + ylim(0, 0.35), 
          align="hv", labels=c("A", "B", "C"), common.legend=T, nrow=1,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot32.png",
         device="png",
         width=11,
         height=4,
         units="in",
         dpi=350,
         bg="white")

# F1, ALL (PLOT 33)
ggarrange(PLOTS$Tel$F1, 
          PLOTS$Cen$F1 + rremove("ylab"),
          PLOTS$SegDup$F1 + rremove("ylab"), 
          align="hv", labels=c("A", "B", "C"), common.legend=T, nrow=1,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot33.png",
         device="png",
         width=11,
         height=4,
         units="in",
         dpi=350,
         bg="white")

# FMI, ALL (PLOT 34)
ggarrange(PLOTS$Tel$FMI, 
          PLOTS$Cen$FMI + rremove("ylab"),
          PLOTS$SegDup$FMI + rremove("ylab"), 
          align="hv", labels=c("A", "B", "C"), common.legend=T, nrow=1,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot34.png",
         device="png",
         width=11,
         height=4,
         units="in",
         dpi=350,
         bg="white")

# JI, ALL (PLOT 35)
ggarrange(PLOTS$Tel$JI, 
          PLOTS$Cen$JI + rremove("ylab"),
          PLOTS$SegDup$JI + rremove("ylab"), 
          align="hv", labels=c("A", "B", "C"), common.legend=T, nrow=1,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot35.png",
         device="png",
         width=11,
         height=4,
         units="in",
         dpi=350,
         bg="white")

# WRITE TABLE DATA
ANALYSIS_REGIONAL_SSC %>%
  rename(Method=Matching_Method,
         MaxD=Matching_Distance,
         QC=Matching_Type) %>%
  select(-MaxD_LOG) %>%
  write_tsv("TABLES/Table10.tsv", col_names=TRUE)                               
