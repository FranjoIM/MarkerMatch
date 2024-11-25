# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)

# LOAD DATA
load("Validation_Final_103024.RData")

# PREPARE A DATA FILE NAMES FOR ANALYSIS
DataFileNames <- data.frame(
  Factor=rep(c("FullSet", "PerfectMatch", "BAF", "LRRmean", "LRRsd", "Pos"), 3),
  MaxD=rep(c(NA, NA, rep("10000", 4)), 3),
  MaxDLab=rep(c("0", "0", rep(c("10000"), 4)), 3),
  Ref=c(rep("SSC", 6), rep("GSA", 6), rep("OEE", 6)),
  Mat=c(rep("SSC", 6), rep("OEE", 6), rep("GSA", 6)),
  stringsAsFactors=FALSE)

# INITIATE HOLDING DF
ANALYSIS_STEP3_SAMPLES <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$MaxDLab[h]
  k <- DataFileNames$Ref[h]
  l <- DataFileNames$Mat[h]
  
  for(o in c("Raw", "QCd")){
    
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
        Factor=i,
        D_MAX=j,
        QC=o,
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

# CHECKPOINT
save(ANALYSIS_STEP3_SAMPLES, file="Analysis_Samplewise_112524.RData")
load("Analysis_Samplewise_112524.RData")

# DEFINE PLOTTING FUNCTION
MetricPlot <- function(a, b, c, d){
  
  TITLE <- paste0(a)
  
  a <- gsub(" ", "", a)
  
  PLOT <- ANALYSIS_STEP3_SAMPLES %>%
    filter(Factor==a) %>%
    filter(D_MAX==b) %>%
    filter(QC==c) %>%
    filter(Mat==d) %>%
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
          plot.subtitle=element_text(size=12, hjust=0, vjust=0, face="bold.italic"),
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

D <- c("GSA", "OEE", "SSC")
D <- set_names(D)

# PLOT METRICS
PLOTS <- map(A, function(x) map(B, function(y) map(C, function(z) map(D, function(w) MetricPlot(a=x, b=y, c=z, d=w)))))

# PLOT RAW CNV CALLS
ggarrange(
  PLOTS$`Full Set`$`0`$QCd$GSA + rremove("xlab"), 
  PLOTS$`Full Set`$`0`$QCd$OEE + rremove("xlab") + rremove("ylab"),
  PLOTS$`Full Set`$`0`$QCd$SSC + rremove("xlab") + rremove("ylab"),
  PLOTS$`BAF`$`10000`$QCd$GSA + rremove("xlab"),
  PLOTS$`BAF`$`10000`$QCd$OEE + rremove("xlab") + rremove("ylab"),
  PLOTS$`BAF`$`10000`$QCd$SSC + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR mean`$`10000`$QCd$GSA + rremove("xlab"),
  PLOTS$`LRR mean`$`10000`$QCd$OEE + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR mean`$`10000`$QCd$SSC + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR sd`$`10000`$QCd$GSA + rremove("xlab"),
  PLOTS$`LRR sd`$`10000`$QCd$OEE + rremove("xlab") + rremove("ylab"),
  PLOTS$`LRR sd`$`10000`$QCd$SSC + rremove("xlab") + rremove("ylab"),
  PLOTS$`Pos`$`10000`$QCd$GSA + rremove("xlab"),
  PLOTS$`Pos`$`10000`$QCd$OEE + rremove("xlab") + rremove("ylab"),
  PLOTS$`Pos`$`10000`$QCd$SSC + rremove("xlab") + rremove("ylab"),
  PLOTS$`Perfect Match`$`0`$QCd$GSA,
  PLOTS$`Perfect Match`$`0`$QCd$OEE +  rremove("ylab"),
  PLOTS$`Perfect Match`$`0`$QCd$SSC +  rremove("ylab"),
  align="hv", labels=c("A", "B", "C", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "), 
  common.legend=T, legend="bottom", nrow=6, ncol=3) %>%
  ggsave(filename="FIGURES/Plot80.png",
         device="png",
         width=7,
         height=7,
         units="in",
         dpi=350,
         bg="white")

ANALYSIS_PPV <- ANALYSIS_STEP3_SAMPLES %>%
  mutate(METRIC=round(TP/(FP+TP), digits=3)) %>%
  arrange(desc(METRIC)) %>%
  filter(!(TP==0 & FP==0)) %>%
  select(-c(ID, D_MAX, Ref, TP, FP, FN)) %>%
  mutate(Tens=case_when(
    METRIC == 0.0 ~ "0.0",
    METRIC >  0.0 & METRIC < 0.1 ~ "(0.0, 0.1)",
    METRIC >= 0.1 & METRIC < 0.2 ~ "[0.1, 0.2)",
    METRIC >= 0.2 & METRIC < 0.3 ~ "[0.2, 0.3)",
    METRIC >= 0.3 & METRIC < 0.4 ~ "[0.3, 0.4)",
    METRIC >= 0.4 & METRIC < 0.5 ~ "[0.4, 0.5)",
    METRIC >= 0.5 & METRIC < 0.6 ~ "[0.5, 0.6)",
    METRIC >= 0.6 & METRIC < 0.7 ~ "[0.6, 0.7)",
    METRIC >= 0.7 & METRIC < 0.8 ~ "[0.7, 0.8)",
    METRIC >= 0.8 & METRIC < 0.9 ~ "[0.8, 0.9)",
    METRIC >= 0.9 & METRIC < 1.0 ~ "[0.9,  1.0)",
    METRIC == 1.0 ~ "1.0",
    TRUE ~ NA_character_)) %>%
  group_by(Mat, Factor, QC, Tens) %>%
  summarize(PPV_C=n())
  
ANALYSIS_F1 <- ANALYSIS_STEP3_SAMPLES %>%
  mutate(METRIC=round((2*TP)/(2*TP+FP+FN), digits=3)) %>%
  arrange(desc(METRIC)) %>%
  filter(!(TP==0 & FP==0 & FN==0)) %>%
  select(-c(ID, D_MAX, Ref, TP, FP, FN)) %>%
  mutate(Tens=case_when(
    METRIC == 0.0 ~ "0.0",
    METRIC >  0.0 & METRIC < 0.1 ~ "(0.0, 0.1)",
    METRIC >= 0.1 & METRIC < 0.2 ~ "[0.1, 0.2)",
    METRIC >= 0.2 & METRIC < 0.3 ~ "[0.2, 0.3)",
    METRIC >= 0.3 & METRIC < 0.4 ~ "[0.3, 0.4)",
    METRIC >= 0.4 & METRIC < 0.5 ~ "[0.4, 0.5)",
    METRIC >= 0.5 & METRIC < 0.6 ~ "[0.5, 0.6)",
    METRIC >= 0.6 & METRIC < 0.7 ~ "[0.6, 0.7)",
    METRIC >= 0.7 & METRIC < 0.8 ~ "[0.7, 0.8)",
    METRIC >= 0.8 & METRIC < 0.9 ~ "[0.8, 0.9)",
    METRIC >= 0.9 & METRIC < 1.0 ~ "[0.9,  1.0)",
    METRIC == 1.0 ~ "1.0",
    TRUE ~ NA_character_)) %>%
  group_by(Mat, Factor, QC, Tens) %>%
  summarize(F1_C=n())

full_join(ANALYSIS_PPV, ANALYSIS_F1, by=c("Mat", "Factor", "QC", "Tens")) %>%
  mutate(PPV_C=ifelse(is.na(PPV_C), 0, PPV_C),
         F1_C=ifelse(is.na(F1_C), 0, F1_C)) %>%
  pivot_wider(names_from=c("Mat", "Factor"), values_from=c("PPV_C", "F1_C"), values_fill=0) %>%
  mutate(QC=ifelse(QC=="Raw", "Low-stringency QC", "Medium-stringency QC")) %>%
  write_tsv("TABLES/TableS1P.tsv", col_names=TRUE)

# PREPARE CUMULATIVE CALCULATIONS OF PPV
ITER_DF <- ANALYSIS_STEP3_SAMPLES %>%
  select(c(Factor, QC, Mat)) %>%
  distinct() %>%
  filter(Factor!="FullSet")

CUM_ANALYSIS <- data.frame(NULL)

for(i in 1:nrow(ITER_DF)){
  DEN <- ANALYSIS_STEP3_SAMPLES %>%
    filter(Factor == ITER_DF$Factor[i],
           QC == ITER_DF$QC[i],
           Mat == ITER_DF$Mat[i]) %>%
    summarize(SumTP = sum(TP),
              SumFP = sum(FP),
              SumFN = sum(FN))
  
  for(j in seq(0.0, 1.0, by=c(0.05))) {
      
      ROW <- data.frame(NULL)
      
      ROW[1, "Factor"] <- ITER_DF$Factor[i]
      ROW[1, "QC"] <- ITER_DF$QC[i]
      ROW[1, "Mat"] <- ITER_DF$Mat[i]
      
      FIL <- ANALYSIS_STEP3_SAMPLES %>%
        filter(Factor == ITER_DF$Factor[i] &
                 QC == ITER_DF$QC[i] &
                 Mat == ITER_DF$Mat[i]) %>%
        mutate(PPV=round(TP/(FP+TP), digits=3)) %>%
        filter(PPV <= j)
      
      NUM <- FIL %>%
        summarize(SumTP = sum(TP),
                  SumFP = sum(FP),
                  SumFN = sum(FN))
      
      COM <- ANALYSIS_STEP3_SAMPLES %>%
        filter(Factor == ITER_DF$Factor[i] &
                 QC == ITER_DF$QC[i] &
                 Mat == ITER_DF$Mat[i]) %>%
        filter(!ID %in% FIL$ID) %>%
        summarize(SumTP = sum(TP),
                  SumFP = sum(FP),
                  SumFN = sum(FN))
      
      ROW[1, "PPV Cutoff"] <- j
      ROW[1, "Samples"] <- nrow(FIL)
      ROW[1, "TP"] <- NUM$SumTP
      ROW[1, "FP"] <- NUM$SumFP
      ROW[1, "FN"] <- NUM$SumFN
      ROW[1, "Under_PPV"] <- round(NUM$SumTP / (NUM$SumFP + NUM$SumTP), digits=3)
      ROW[1, "Over_PPV"] <- round(COM$SumTP / (COM$SumFP + COM$SumTP), digits=3)
      ROW[1, "Under_F1"] <- round(2*NUM$SumTP / (NUM$SumFP + 2*NUM$SumTP + NUM$SumFN), digits=3)
      ROW[1, "Over_F1"] <- round(2*COM$SumTP / (COM$SumFP + 2*COM$SumTP + NUM$SumFN), digits=3)
      
      CUM_ANALYSIS <- rbind(CUM_ANALYSIS, ROW)
  }
}

CUM_ANALYSIS <- CUM_ANALYSIS %>%
   mutate(QC=ifelse(QC=="Raw", "Low-stringency QC", "Medium-stringency QC"),
         Factor=ifelse(Factor=="FullSet", "Full Set",
                ifelse(Factor=="LRRmean", "LRR mean",
                ifelse(Factor=="LRRsd", "LRR sd",
                ifelse(Factor=="Pos", "Distance",
                ifelse(Factor=="PerfectMatch", "Perfect Match", Factor)))))) 

write_tsv(CUM_ANALYSIS, "TABLES/TableS1Q.tsv", col_names=TRUE)
  
# PLOT THE CURVES
CUM_ANALYSIS$FactorF <- factor(CUM_ANALYSIS$Factor,
                               levels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match"),
                               labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match"))

(CUM_ANALYSIS %>%
  filter(QC=="Medium-stringency QC") %>%
  ggplot(aes(x=`PPV Cutoff`, color=FactorF)) +
  geom_line(aes(y=Under_PPV, linetype="Fail"), linewidth=1) +
  geom_line(aes(y=Over_PPV, linetype="Pass"), linewidth=1) +
  labs(x="PPV Cutoff",
       y="PPV",
       linetype="PPV QC",
       color="FACTOR") +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  facet_wrap(. ~ Mat, ncol=3) +
  theme_bw() +
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        plot.subtitle=element_text(size=12, hjust=0, vjust=0, face="bold.italic"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.caption=element_text(size=12, face="bold.italic"),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff")) +
  guides(color=guide_legend(title.position="top", nrow=1, order=1),
         linetype=guide_legend(title.position="top", nrow=1, order=1))) %>%
  ggsave(filename="FIGURES/Plot81.png",
         device="png",
         width=7,
         height=7,
         units="in",
         dpi=350,
         bg="white")

(CUM_ANALYSIS %>%
    filter(QC=="Medium-stringency QC") %>%
    ggplot(aes(x=`PPV Cutoff`, color=FactorF)) +
    geom_line(aes(y=Under_F1, linetype="Fail"), linewidth=1) +
    geom_line(aes(y=Over_F1, linetype="Pass"), linewidth=1) +
    labs(x="PPV Cutoff",
         y="F1",
         linetype="PPV QC",
         color="FACTOR") +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
    facet_wrap(. ~ Mat, ncol=3) +
    theme_bw() +
    theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          plot.subtitle=element_text(size=12, hjust=0, vjust=0, face="bold.italic"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold"),
          plot.caption=element_text(size=12, face="bold.italic"),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff")) +
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           linetype=guide_legend(title.position="top", nrow=1, order=1))) %>%
  ggsave(filename="FIGURES/Plot82.png",
         device="png",
         width=7,
         height=7,
         units="in",
         dpi=350,
         bg="white")
