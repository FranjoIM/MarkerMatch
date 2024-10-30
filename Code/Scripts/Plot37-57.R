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
ANALYSIS_STEP2 <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$MaxDLab[h]
  k <- DataFileNames$Ref[h]
  l <- DataFileNames$Mat[h]
  
  # Pull out CNVs overlapping in full set and partial sets
  for(o in c("Raw", "QCd")){
    REF <- DATA[["Raw"]][[k]][["FullSet"]][["0"]][["CNV"]]
    MAT <- DATA[[o]][[l]][[i]][[j]][["CNV"]]
    
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
    
    for(p in c("All", "Deletions", "Duplications")){
      for(q in c("All", "Small", "Medium", "Large")){
        
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
          Ref=k,
          Mat=l,
          CNV_Type=p,
          CNV_Size=q,
          TP=length(TP),
          FP=length(FP),
          FN=length(FN))
        
        ANALYSIS_STEP2 <- rbind(ANALYSIS_STEP2, NEW_ROW)
        
        timestamp()
        message("Finished ", i, "-", j, "-", o, "-", p, "-", q, ".\n")
        
      }
    }
  }
}

# TIDY UP ANALYSIS
ANALYSIS_STEP2 <- ANALYSIS_STEP2 %>%
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
MetricPlot <- function(a, b, c, d){
  H1 <- ANALYSIS_STEP2 %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Type==a & CNV_Size==b & Mat==d & Matching_Type=="Raw") %>%
    pull(.data[[c]])
  
  H2 <- ANALYSIS_STEP2 %>%
    filter(Matching_Method=="PerfectMatch" & CNV_Type==a & CNV_Size==b & Mat==d & Matching_Type=="QCd") %>%
    pull(.data[[c]])
  
  H3 <- ANALYSIS_STEP2 %>%
    filter(Matching_Method=="FullSet" & CNV_Type==a & CNV_Size==b & Mat==d & Matching_Type=="Raw") %>%
    pull(.data[[c]])
  
  H4 <- ANALYSIS_STEP2 %>%
    filter(Matching_Method=="FullSet" & CNV_Type==a & CNV_Size==b & Mat==d & Matching_Type=="QCd") %>%
    pull(.data[[c]])
  
  TITLE <- case_when(
    b == "All" ~ "All CNVs",
    b == "Small" ~ "CNV < 100kb",
    b == "Medium" ~ "100kb < CNV < 500kb",
    b == "Large" ~ "500kb < CNV",
    TRUE ~ NA_character_)
  
  PLOT <- ANALYSIS_STEP2 %>%
    filter(Mat==d) %>%
    filter(!Matching_Method %in% c("PerfectMatch", "FullSet")) %>%
    filter(CNV_Type==a & CNV_Size==b) %>%
    ggplot() +
    geom_hline(aes(yintercept=H1, color="Perfect Match", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H2, color="Perfect Match", linetype="QCd"), linewidth=1) +
    geom_hline(aes(yintercept=H3, color="Full Set", linetype="Raw"), linewidth=1) +
    geom_hline(aes(yintercept=H4, color="Full Set", linetype="QCd"), linewidth=1) +
    geom_point(aes(x=MaxD_LOG, y=.data[[c]], color=Matching_Method, shape="Raw"), 
               data = ~filter(.x, Matching_Type=="Raw"),
               position = position_dodge(0.01), size=3) +
    geom_point(aes(x=MaxD_LOG, y=.data[[c]], color=Matching_Method, shape="QCd"), 
               data = ~filter(.x, Matching_Type=="QCd"), 
               position = position_dodge(0.01), size=3.5) +
    scale_shape_manual(values=c(2, 17),
                       breaks=c("Raw", "QCd")) +
    scale_linetype_manual(values=c("dashed", "solid"),
                       breaks=c("Raw", "QCd")) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRRmean", "LRRsd", "Pos", "Perfect Match", "Full Set")) +
    labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
         y=toupper(c),
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
A <- c("All", "Deletions", "Duplications")
A <- set_names(A)

B <- c("All", "Small", "Medium", "Large")
B <- set_names(B)

C <- c("Sensitivity", "PPV", "FNR", "FDR", "F1", "FMI", "JI")
C <- set_names(C)

D <- c("GSA", "OEE")
D <- set_names(D)

# PLOT METRICS
PLOTS <- map(A, function(x) map(B, function(y) map(C, function(z) map(D, function(w) MetricPlot(a=x, b=y, c=z, d=w)))))

# SENSITIVITY (PLOT 37)
ggarrange(
  PLOTS$All$All$Sensitivity$GSA + rremove("xlab"), 
  PLOTS$All$Small$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$All$Medium$Sensitivity$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$All$Large$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),

  PLOTS$All$All$Sensitivity$OEE, 
  PLOTS$All$Small$Sensitivity$OEE + rremove("ylab"),
  PLOTS$All$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$All$Large$Sensitivity$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot37.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# SENSITIVITY DELETIONS (PLOT 38)
ggarrange(
  PLOTS$Deletions$All$Sensitivity$GSA + rremove("xlab"), 
  PLOTS$Deletions$Small$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivity$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Deletions$Large$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Deletions$All$Sensitivity$OEE, 
  PLOTS$Deletions$Small$Sensitivity$OEE + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Deletions$Large$Sensitivity$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot38.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# SENSITIVITY DELETIONS (PLOT 39)
ggarrange(
  PLOTS$Duplications$All$Sensitivity$GSA + rremove("xlab"), 
  PLOTS$Duplications$Small$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivity$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Duplications$Large$Sensitivity$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Duplications$All$Sensitivity$OEE, 
  PLOTS$Duplications$Small$Sensitivity$OEE + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Duplications$Large$Sensitivity$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot39.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# PPV (PLOT 40)
ggarrange(
  PLOTS$All$All$PPV$GSA + rremove("xlab"), 
  PLOTS$All$Small$PPV$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$All$Medium$PPV$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$All$Large$PPV$GSA + rremove("xlab") + rremove("ylab"),

  PLOTS$All$All$PPV$OEE, 
  PLOTS$All$Small$PPV$OEE + rremove("ylab"),
  PLOTS$All$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$All$Large$PPV$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot40.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# PPV DELETIONS (PLOT 41)
ggarrange(
  PLOTS$Deletions$All$PPV$GSA + rremove("xlab"), 
  PLOTS$Deletions$Small$PPV$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Deletions$Medium$PPV$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Deletions$Large$PPV$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Deletions$All$PPV$OEE, 
  PLOTS$Deletions$Small$PPV$OEE + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Deletions$Large$PPV$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot41.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# PPV DELETIONS (PLOT 42)
ggarrange(
  PLOTS$Duplications$All$PPV$GSA + rremove("xlab"), 
  PLOTS$Duplications$Small$PPV$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Duplications$Medium$PPV$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Duplications$Large$PPV$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Duplications$All$PPV$OEE, 
  PLOTS$Duplications$Small$PPV$OEE + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Duplications$Large$PPV$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot42.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FNR (PLOT 43)
ggarrange(
  PLOTS$All$All$FNR$GSA + rremove("xlab"), 
  PLOTS$All$Small$FNR$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$All$Medium$FNR$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$All$Large$FNR$GSA + rremove("xlab") + rremove("ylab"),

  PLOTS$All$All$FNR$OEE, 
  PLOTS$All$Small$FNR$OEE + rremove("ylab"),
  PLOTS$All$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$All$Large$FNR$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot43.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FNR DELETIONS (PLOT 44)
ggarrange(
  PLOTS$Deletions$All$FNR$GSA + rremove("xlab"), 
  PLOTS$Deletions$Small$FNR$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Deletions$Medium$FNR$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Deletions$Large$FNR$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Deletions$All$FNR$OEE, 
  PLOTS$Deletions$Small$FNR$OEE + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Deletions$Large$FNR$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot44.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FNR DELETIONS (PLOT 45)
ggarrange(
  PLOTS$Duplications$All$FNR$GSA + rremove("xlab"), 
  PLOTS$Duplications$Small$FNR$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Duplications$Medium$FNR$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Duplications$Large$FNR$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Duplications$All$FNR$OEE, 
  PLOTS$Duplications$Small$FNR$OEE + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Duplications$Large$FNR$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot45.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FDR (PLOT 46)
ggarrange(
  PLOTS$All$All$FDR$GSA + rremove("xlab"), 
  PLOTS$All$Small$FDR$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$All$Medium$FDR$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$All$Large$FDR$GSA + rremove("xlab") + rremove("ylab"),

  PLOTS$All$All$FDR$OEE, 
  PLOTS$All$Small$FDR$OEE + rremove("ylab"),
  PLOTS$All$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$All$Large$FDR$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot46.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FDR DELETIONS (PLOT 47)
ggarrange(
  PLOTS$Deletions$All$FDR$GSA + rremove("xlab"), 
  PLOTS$Deletions$Small$FDR$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Deletions$Medium$FDR$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Deletions$Large$FDR$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Deletions$All$FDR$OEE, 
  PLOTS$Deletions$Small$FDR$OEE + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Deletions$Large$FDR$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot47.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FDR DELETIONS (PLOT 48)
ggarrange(
  PLOTS$Duplications$All$FDR$GSA + rremove("xlab"), 
  PLOTS$Duplications$Small$FDR$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Duplications$Medium$FDR$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Duplications$Large$FDR$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Duplications$All$FDR$OEE, 
  PLOTS$Duplications$Small$FDR$OEE + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Duplications$Large$FDR$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot48.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# F1 (PLOT 49)
ggarrange(
  PLOTS$All$All$F1$GSA + rremove("xlab"), 
  PLOTS$All$Small$F1$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$All$Medium$F1$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$All$Large$F1$GSA + rremove("xlab") + rremove("ylab"),

  PLOTS$All$All$F1$OEE, 
  PLOTS$All$Small$F1$OEE + rremove("ylab"),
  PLOTS$All$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$All$Large$F1$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot49.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# F1 DELETIONS (PLOT 50)
ggarrange(
  PLOTS$Deletions$All$F1$GSA + rremove("xlab"), 
  PLOTS$Deletions$Small$F1$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Deletions$Medium$F1$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Deletions$Large$F1$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Deletions$All$F1$OEE, 
  PLOTS$Deletions$Small$F1$OEE + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Deletions$Large$F1$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot50.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# F1 DELETIONS (PLOT 51)
ggarrange(
  PLOTS$Duplications$All$F1$GSA + rremove("xlab"), 
  PLOTS$Duplications$Small$F1$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Duplications$Medium$F1$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Duplications$Large$F1$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Duplications$All$F1$OEE, 
  PLOTS$Duplications$Small$F1$OEE + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Duplications$Large$F1$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot51.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FMI (PLOT 52)
ggarrange(
  PLOTS$All$All$FMI$GSA + rremove("xlab"), 
  PLOTS$All$Small$FMI$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$All$Medium$FMI$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$All$Large$FMI$GSA + rremove("xlab") + rremove("ylab"),

  PLOTS$All$All$FMI$OEE, 
  PLOTS$All$Small$FMI$OEE + rremove("ylab"),
  PLOTS$All$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$All$Large$FMI$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot52.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FMI DELETIONS (PLOT 53)
ggarrange(
  PLOTS$Deletions$All$FMI$GSA + rremove("xlab"), 
  PLOTS$Deletions$Small$FMI$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Deletions$Medium$FMI$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Deletions$Large$FMI$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Deletions$All$FMI$OEE, 
  PLOTS$Deletions$Small$FMI$OEE + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Deletions$Large$FMI$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot53.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# FMI DELETIONS (PLOT 54)
ggarrange(
  PLOTS$Duplications$All$FMI$GSA + rremove("xlab"), 
  PLOTS$Duplications$Small$FMI$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Duplications$Medium$FMI$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Duplications$Large$FMI$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Duplications$All$FMI$OEE, 
  PLOTS$Duplications$Small$FMI$OEE + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Duplications$Large$FMI$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot54.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

 # JI (PLOT 55)
ggarrange(
  PLOTS$All$All$JI$GSA + rremove("xlab"), 
  PLOTS$All$Small$JI$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$All$Medium$JI$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$All$Large$JI$GSA + rremove("xlab") + rremove("ylab"),

  PLOTS$All$All$JI$OEE, 
  PLOTS$All$Small$JI$OEE + rremove("ylab"),
  PLOTS$All$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$All$Large$JI$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot55.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# JI DELETIONS (PLOT 56)
ggarrange(
  PLOTS$Deletions$All$JI$GSA + rremove("xlab"), 
  PLOTS$Deletions$Small$JI$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Deletions$Medium$JI$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Deletions$Large$JI$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Deletions$All$JI$OEE, 
  PLOTS$Deletions$Small$JI$OEE + rremove("ylab"),
  PLOTS$Deletions$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Deletions$Large$JI$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot56.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white") 

# JI DELETIONS (PLOT 57)
ggarrange(
  PLOTS$Duplications$All$JI$GSA + rremove("xlab"), 
  PLOTS$Duplications$Small$JI$GSA + rremove("xlab") + rremove("ylab"),
  PLOTS$Duplications$Medium$JI$GSA + rremove("xlab") + rremove("ylab"), 
  PLOTS$Duplications$Large$JI$GSA + rremove("xlab") + rremove("ylab"),
  
  PLOTS$Duplications$All$JI$OEE, 
  PLOTS$Duplications$Small$JI$OEE + rremove("ylab"),
  PLOTS$Duplications$Medium$Sensitivit$OEE + rremove("ylab"), 
  PLOTS$Duplications$Large$JI$OEE + rremove("ylab"),
  align="hv", labels=c("A", " ", " ", " ", "B", " ", " ", " "), common.legend=T,
  legend="top", nrow=2, ncol=4) %>%
  ggsave(filename="FIGURES/Plot57.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

# WRITE TABLE DATA
ANALYSIS_STEP2 %>%
  rename(Method=Matching_Method,
         MaxD=Matching_Distance,
         QC=Matching_Type) %>%
  select(-MaxD_LOG) %>%
  write_tsv("TABLES/Table11.tsv", col_names=TRUE)
