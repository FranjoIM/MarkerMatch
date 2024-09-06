# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)

# LOAD DATA
load("Validation_Final_082924.RData")

# PREPARE A DATA FILE NAMES FOR SSC
DataFileNames <- data.frame(
  Factor=rep(c("FullSet", "PerfectMatch", "BAF", "LRRmean", "LRRsd", "Pos"), 2),
  MaxD=rep(c(NA, NA, rep("10000", 4)), 2),
  MaxDLab=rep(c("0", "0", rep(c("10000"), 4)), 2),
  Ref=c(rep("GSA", 6), rep("OEE", 6)),
  Mat=c(rep("OEE", 6), rep("GSA", 6)),
  stringsAsFactors=FALSE)

# INITIATE HOLDING DF
ANALYSIS_STEP2_REGIONAL <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$MaxDLab[h]
  k <- DataFileNames$Ref[h]
  l <- DataFileNames$Mat[h]
  
  # Pull out CNVs overlapping in full set and partial sets
  for(o in c("Raw", "QCd")){
    REF <- DATA[["Raw"]][[k]][["FullSet"]][["0"]][["CNV"]] %>%
      mutate(N_SNP = as.numeric(NumSNP),
             LEN = as.numeric(gsub(",", "", Length)),
             CHR = as.numeric(CHR),
             START = as.numeric(START),
             END = as.numeric(END),
             STATE = as.numeric(STATE),
             CN = as.numeric(CN)) %>%
      select(-c(StartSNP, EndSNP, NumSNP, Length))
    MAT <- DATA[[o]][[l]][[i]][[j]][["CNV"]] %>%
      mutate(N_SNP = as.numeric(NumSNP),
             LEN = as.numeric(gsub(",", "", Length)),
             CHR = as.numeric(CHR),
             START = as.numeric(START),
             END = as.numeric(END),
             STATE = as.numeric(STATE),
             CN = as.numeric(CN)) %>%
      select(-c(StartSNP, EndSNP, NumSNP, Length))
    
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
        !is.na(LEN.y) & LEN.y>=500000 ~ "Large",
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
        Ref=k,
        Mat=l,
        Matching_Type=o,
        CNV_Region=p,
        TP=length(TP),
        FP=length(FP),
        FN=length(FN))
      
      ANALYSIS_STEP2_REGIONAL <- rbind(ANALYSIS_STEP2_REGIONAL, NEW_ROW)
      
      timestamp()
      message("Finished ", i, "-", j, "-", k, "-", l,  "-", o, "-", p, ".\n")
    }
  }
}

# TIDY UP ANALYSIS
ANALYSIS_STEP2_REGIONAL <- ANALYSIS_STEP2_REGIONAL %>%
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
MetricPlot <- function(a, b, c){
  H1 <- ANALYSIS_STEP2_REGIONAL %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Region==a & Mat==c & Matching_Type=="Raw") %>%
    pull(.data[[b]])
  
  H2 <- ANALYSIS_STEP2_REGIONAL %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Region==a & Mat==c & Matching_Type=="QCd") %>%
    pull(.data[[b]])
  
  H3 <- ANALYSIS_STEP2_REGIONAL %>%
    filter(Matching_Method=="FullSet" & CNV_Region==a & Mat==c & Matching_Type=="Raw") %>%
    pull(.data[[b]])
  
  H4 <- ANALYSIS_STEP2_REGIONAL %>%
    filter(Matching_Method=="FullSet" & CNV_Region==a & Mat==c & Matching_Type=="QCd") %>%
    pull(.data[[b]])
  
  TITLE <- case_when(
    a == "Tel" ~ "Telomeric",
    a == "Cen" ~ "Centromeric",
    a == "SegDup" ~ "Segmental Duplications",
    TRUE ~ NA_character_)
  
  PLOT <- ANALYSIS_STEP2_REGIONAL %>%
    filter(Mat==c) %>%
    filter(!Matching_Method %in% c("PerfectMatch", "FullSet")) %>%
    filter(CNV_Region==a) %>%
    ggplot() +
    geom_hline(aes(yintercept=H1, color="Perfect Match", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H2, color="Perfect Match", linetype="QCd"), linewidth=1) +
    geom_hline(aes(yintercept=H3, color="Full Set", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H4, color="Full Set", linetype="QCd"), linewidth=1) +
    geom_point(aes(x=MaxD_LOG, y=.data[[b]], color=Matching_Method, shape="Raw"), 
               data = ~filter(.x, Matching_Type=="Raw"),
               position = position_dodge(0.01), size=3) +
    geom_point(aes(x=MaxD_LOG, y=.data[[b]], color=Matching_Method, shape="QCd"), 
               data = ~filter(.x, Matching_Type=="QCd"), 
               position = position_dodge(0.01), size=3.5) +
    scale_shape_manual(values=c(2, 17),
                       breaks=c("Raw", "QCd")) +
    scale_linetype_manual(values=c("dashed", "solid"),
                          breaks=c("Raw", "QCd")) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRRmean", "LRRsd", "Pos", "Perfect Match", "Full Set")) +
    labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
         y=toupper(b),
         linetype="CNV CALLSET QC",
         color="MATCHING METHOD",
         shape=" ",
         subtitle=TITLE) +
    scale_x_continuous(breaks=4, n.breaks=1, labels=4, limits=c(3.98, 4.02)) +
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
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           linetype=guide_legend(title.position="top", order=2),
           shape=guide_legend(title.position="top", order=3))
}

# PREPARE FOR PLOTTING
A <- c("Tel", "Cen", "SegDup")
A <- set_names(A)

B <- c("Sensitivity", "PPV", "FNR", "FDR", "F1", "FMI", "JI")
B <- set_names(B)

C <- c("GSA", "OEE")
C <- set_names(C)

# PLOT METRICS
PLOTS <- map(A, function(x) map(B, function(y) map(C, function(z) MetricPlot(a=x, b=y, c=z))))

# SENSITIVITY, ALL (PLOT 58)
ggarrange(PLOTS$Tel$Sensitivity$GSA + rremove("xlab"), 
          PLOTS$Cen$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),
          PLOTS$SegDup$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),
          
          PLOTS$Tel$Sensitivity$OEE, 
          PLOTS$Cen$Sensitivity$OEE + rremove("ylab"),
          PLOTS$SegDup$Sensitivity$OEE + rremove("ylab"),
          
          align="hv", labels=c("A", " ", " ", "B", " ", " "), common.legend=T, nrow=2, ncol=3,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot58.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

# PPV, ALL (PLOT 59)
ggarrange(PLOTS$Tel$PPV$GSA + rremove("xlab"), 
          PLOTS$Cen$PPV$GSA + rremove("xlab") + rremove("ylab"),
          PLOTS$SegDup$PPV$GSA + rremove("xlab") + rremove("ylab"),
          
          PLOTS$Tel$PPV$OEE, 
          PLOTS$Cen$PPV$OEE + rremove("ylab"),
          PLOTS$SegDup$PPV$OEE + rremove("ylab"),
          
          align="hv", labels=c("A", " ", " ", "B", " ", " "), common.legend=T, nrow=2, ncol=3,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot59.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")


# FNR, ALL (PLOT 60)
ggarrange(PLOTS$Tel$FNR$GSA + rremove("xlab"), 
          PLOTS$Cen$FNR$GSA + rremove("xlab") + rremove("ylab"),
          PLOTS$SegDup$FNR$GSA + rremove("xlab") + rremove("ylab"),
          
          PLOTS$Tel$FNR$OEE, 
          PLOTS$Cen$FNR$OEE + rremove("ylab"),
          PLOTS$SegDup$FNR$OEE + rremove("ylab"),
          
          align="hv", labels=c("A", " ", " ", "B", " ", " "), common.legend=T, nrow=2, ncol=3,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot60.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

# FDR, ALL (PLOT 61)
ggarrange(PLOTS$Tel$FDR$GSA + rremove("xlab"), 
          PLOTS$Cen$FDR$GSA + rremove("xlab") + rremove("ylab"),
          PLOTS$SegDup$FDR$GSA + rremove("xlab") + rremove("ylab"),
          
          PLOTS$Tel$FDR$OEE, 
          PLOTS$Cen$FDR$OEE + rremove("ylab"),
          PLOTS$SegDup$FDR$OEE + rremove("ylab"),
          
          align="hv", labels=c("A", " ", " ", "B", " ", " "), common.legend=T, nrow=2, ncol=3,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot61.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

# F1, ALL (PLOT 62)
ggarrange(PLOTS$Tel$F1$GSA + rremove("xlab"), 
          PLOTS$Cen$F1$GSA + rremove("xlab") + rremove("ylab"),
          PLOTS$SegDup$F1$GSA + rremove("xlab") + rremove("ylab"),
          
          PLOTS$Tel$F1$OEE, 
          PLOTS$Cen$F1$OEE + rremove("ylab"),
          PLOTS$SegDup$F1$OEE + rremove("ylab"),
          
          align="hv", labels=c("A", " ", " ", "B", " ", " "), common.legend=T, nrow=2, ncol=3,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot62.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

# FMI, ALL (PLOT 63)
ggarrange(PLOTS$Tel$FMI$GSA + rremove("xlab"), 
          PLOTS$Cen$FMI$GSA + rremove("xlab") + rremove("ylab"),
          PLOTS$SegDup$FMI$GSA + rremove("xlab") + rremove("ylab"),
          
          PLOTS$Tel$FMI$OEE, 
          PLOTS$Cen$FMI$OEE + rremove("ylab"),
          PLOTS$SegDup$FMI$OEE + rremove("ylab"),
          
          align="hv", labels=c("A", " ", " ", "B", " ", " "), common.legend=T, nrow=2, ncol=3,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot63.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

# JI, ALL (PLOT 64)
ggarrange(PLOTS$Tel$JI$GSA + rremove("xlab"), 
          PLOTS$Cen$JI$GSA + rremove("xlab") + rremove("ylab"),
          PLOTS$SegDup$JI$GSA + rremove("xlab") + rremove("ylab"),
          
          PLOTS$Tel$JI$OEE, 
          PLOTS$Cen$JI$OEE + rremove("ylab"),
          PLOTS$SegDup$JI$OEE + rremove("ylab"),
          
          align="hv", labels=c("A", " ", " ", "B", " ", " "), common.legend=T, nrow=2, ncol=3,
          legend="top") %>%
  ggsave(filename="FIGURES/Plot64.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

# WRITE TABLE DATA
ANALYSIS_STEP2_REGIONAL %>%
  rename(Method=Matching_Method,
         MaxD=Matching_Distance,
         QC=Matching_Type) %>%
  select(-MaxD_LOG) %>%
  write_tsv("TABLES/Table12.tsv", col_names=TRUE)  
