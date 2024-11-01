# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)

# LOAD DATA
load("Validation_Final_103024.RData")

# PREPARE A DATA FILE NAMES FOR SSC
DataFileNames <- data.frame(
  Factor=rep(c("FullSet", "PerfectMatch", "BAF", "LRRmean", "LRRsd", "Pos"), 2),
  MaxD=rep(c(NA, NA, rep("10000", 4)), 2),
  MaxDLab=rep(c("0", "0", rep(c("10000"), 4)), 2),
  Ref=c(rep("GSA", 6), rep("OEE", 6)),
  Mat=c(rep("OEE", 6), rep("GSA", 6)),
  stringsAsFactors=FALSE)

# INITIATE HOLDING DF
ANALYSIS_STEP3_SAMPLES <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$MaxDLab[h]
  k <- DataFileNames$Ref[h]
  l <- DataFileNames$Mat[h]
  
  for(o in CNV_Types){
    
    # Pull out CNVs overlapping in full set and partial sets
    REF <- DATA[["Raw"]][[k]][["FullSet"]][["0"]][["CNV"]]
    MAT <- DATA[[o]][[l]][[i]][[j]][["CNV"]]
    
    # Keep only same samples for the analysis
    KEEP_IDs <- intersect(unique(REF$ID), unique(MAT$ID))
    
    if (o=="QCd") {
      MAT <- MAT %>%
        filter(Tel==0 & Cen==0 & Imu==0)}
    
    for(m in KEEP_IDs){
      # Join the CNVs into ID, CHR, STATE, and CN matched DF
      # Classify CNVs by size and by FP/TP/FN/TN status
      COM <- full_join(REF, MAT, 
                       by=c("ID", "CHR", "STATE", "CN"), 
                       relationship="many-to-many") %>%
        filter(ID==m) %>%
        mutate(cIDx=CNV_ID.x,
               cIDy=CNV_ID.y,
               ConfM=case_when(
                 is.na(START.y) ~ "FN",
                 is.na(START.x) ~ "FP", 
                 START.x<=END.y & END.x>=START.y ~ "TP",
                 START.x>END.y | END.x<START.y ~ "FP",
                 TRUE ~ "Edge Case"))
      
      if(nrow(COM)==0) {message("Sample has no CNV calls.\n")}
      
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
      
      FILT <- rbind(COM_TEMP1, COM_TEMP2, COM_TEMP3, COM_TEMP4)
      
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
        ID=m,
        Matching_Method=i,
        Matching_Distance=j,
        Matching_Type=o,
        Ref=k,
        Mat=l,
        TP=length(TP),
        FP=length(FP),
        FN=length(FN))
      
      ANALYSIS_STEP3_SAMPLES <- bind_rows(ANALYSIS_STEP3_SAMPLES, NEW_ROW)
      
      timestamp()
      message("Finished ", i, "-", j, "-", l, "-", m,"-", o, ".\n")
    }
  }
}

# DEFINE PLOTTING FUNCTION
MetricPlot <- function(a, b, c, d){
  
  TITLE <- paste0(a)
  
  a <- gsub(" ", "", a)
  
  PLOT <- ANALYSIS_STEP3_SAMPLES %>%
    filter(Matching_Method==a) %>%
    filter(Matching_Distance==b) %>%
    filter(Matching_Type==c) %>%
    filter(Ref==d) %>%
    mutate(PPV=round(TP/(FP+TP), digits=3)) %>%
    arrange(desc(PPV)) %>%
    mutate(ID2 = row_number()) %>%
    pivot_longer(c(TP, FP), names_to="Call_Type", values_to="Call_Type_Count") %>%
    ggplot(aes(x=as.factor(ID2), y=Call_Type_Count, color=Call_Type, fill=Call_Type)) +
    geom_bar(position="fill", stat="identity") +
    labs(x=NULL,
         y=NULL,
         color="Call Type",
         fill="Call Type",
         subtitle=TITLE) +
    scale_color_manual(values=c("red3", "steelblue3"),
                       breaks=c("FP", "TP"),
                       labels=c("False Positive", "True Positive")) +
    scale_fill_manual(values=c("red3", "steelblue3"),
                      breaks=c("FP", "TP"),
                      labels=c("False Positive", "True Positive")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold"),
          plot.caption=element_text(size=12, face="bold.italic"),
          axis.ticks=element_blank(),
          panel.grid=element_blank()) +
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           fill=guide_legend(title.position="top", nrow=1, order=1))
}

# PREPARE FOR PLOTTING
A <- c("Full Set", "Perfect Match", "BAF", "LRR mean", "LRR sd", "Pos")
A <- set_names(A)

B <- c("0", "10000")
B <- set_names(B)

C <- c("Raw", "QCd")
C <- set_names(C)

D <- c("GSA", "OEE")
D <- set_names(D)

# PLOT METRICS
PLOTS <- map(A, function(x) map(B, function(y) map(C, function(z) map(D, function(w) MetricPlot(a=x, b=y, c=z, d=w)))))

# PLOT RAW CNV CALLS
ggarrange(
  PLOTS$`Full Set`$`0`$Raw$OEE + rremove("xlab"), 
  PLOTS$`Full Set`$`0`$Raw$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`BAF`$`10000`$Raw$OEE + rremove("xlab"),
  PLOTS$`BAF`$`10000`$Raw$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR mean`$`10000`$Raw$OEE + rremove("xlab"),
  PLOTS$`LRR mean`$`10000`$Raw$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR sd`$`10000`$Raw$OEE + rremove("xlab"),
  PLOTS$`LRR sd`$`10000`$Raw$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`Pos`$`10000`$Raw$OEE + rremove("xlab"),
  PLOTS$`Pos`$`10000`$Raw$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`Perfect Match`$`0`$Raw$OEE,
  PLOTS$`Perfect Match`$`0`$Raw$GSA +  rremove("ylab"),
  align="hv", labels=c("A", "B", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "), common.legend=T,
  legend="bottom", nrow=6, ncol=2) %>%
  ggsave(filename="FIGURES/Plot68.png",
         device="png",
         width=7,
         height=11,
         units="in",
         dpi=350,
         bg="white")

# PLOT QCD CNV CALLS
ggarrange(
  PLOTS$`Full Set`$`0`$QCd$OEE + rremove("xlab"), 
  PLOTS$`Full Set`$`0`$QCd$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`BAF`$`10000`$QCd$OEE + rremove("xlab"),
  PLOTS$`BAF`$`10000`$QCd$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR mean`$`10000`$QCd$OEE + rremove("xlab"),
  PLOTS$`LRR mean`$`10000`$QCd$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR sd`$`10000`$QCd$OEE + rremove("xlab"),
  PLOTS$`LRR sd`$`10000`$QCd$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`Pos`$`10000`$QCd$OEE + rremove("xlab"),
  PLOTS$`Pos`$`10000`$QCd$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$`Perfect Match`$`0`$QCd$OEE,
  PLOTS$`Perfect Match`$`0`$QCd$GSA +  rremove("ylab"),
  align="hv", labels=c("A", "B", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "), common.legend=T,
  legend="bottom", nrow=6, ncol=2) %>%
  ggsave(filename="FIGURES/Plot69.png",
         device="png",
         width=7,
         height=11,
         units="in",
         dpi=350,
         bg="white")
