# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)

# LOAD DATA
load("Validation_Final_082924.RData")

DataFileNames <- data.frame(
  Factor=c("FullSet", "PerfectMatch", "BAF", "LRRmean", "LRRsd", "Pos"),
  MaxD=c(NA, NA, rep("10000", 4)),
  MaxDLab=c("0", "0", rep(c("10000"), 4)),
  Ref=rep("GSA", 6),
  Mat=rep("OEE", 6),
  stringsAsFactors=FALSE)

CNV_Lengths <- seq(from=20000, to=500000, by=20000)
SNP_Numbers <- seq(from=10, to=50, by=5)
CNV_Types <- c("Raw", "QCd")

ANALYSIS_STEP2_CUTOFFS <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$MaxDLab[h]
  k <- DataFileNames$Ref[h]
  l <- DataFileNames$Mat[h]
  
  for(m in CNV_Lengths){
    for(n in SNP_Numbers){
      for(o in CNV_Types){
       
        # Pull out CNVs overlapping in full set and partial sets
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
          filter(LEN >= m & N_SNP >= n) %>%
          select(-c(StartSNP, EndSNP, NumSNP, Length))
        
        # Keep only same samples for the analysis
        KEEP_IDs <- intersect(unique(REF$ID), unique(MAT$ID))
        REF <- REF %>%
          filter(ID %in% KEEP_IDs)
        MAT <- MAT %>%
          filter(ID %in% KEEP_IDs)
        
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
          Matching_Method=i,
          Matching_Distance=j,
          Matching_Type=o,
          Ref=k,
          Mat=l,
          LEN_Cutoff=m,
          N_SNP_Cutoff=n,
          TP=length(TP),
          FP=length(FP),
          FN=length(FN))
        
        ANALYSIS_STEP2_CUTOFFS <- bind_rows(ANALYSIS_STEP2_CUTOFFS, NEW_ROW)
        
        timestamp()
        message("Finished ", i, "-", j, "-", l, "-", m, "-", n, "-", o, ".\n")
         
      }
    }
  }
}

# TIDY UP ANALYSIS
ANALYSIS_STEP2_CUTOFFS <- ANALYSIS_STEP2_CUTOFFS %>%
  mutate(Sensitivity=round(TP/(TP+FN), digits=3),
         PPV=round(TP/(FP+TP), digits=3),
         FNR=round(FN/(FN+TP), digits=3),
         FDR=round(FP/(FP+TP), digits=3),
         F1=round((2*TP)/(2*TP+FP+FN), digits=3),
         FMI=round(sqrt((TP/(FP+TP))*(TP/(TP+FN))), digits=3),
         JI=round(TP/(TP+FN+FP), digits=3)) %>%
  mutate(Matching_Distance=as.numeric(Matching_Distance)) %>%
  mutate(MaxD_LOG=log10(Matching_Distance))

# MODEL PPV ON LEN AND N_SNP CUTOFFS
Model <- lm(PPV ~ LEN_Cutoff + N_SNP_Cutoff + Matching_Method, 
            data=filter(ANALYSIS_STEP2_CUTOFFS, !Matching_Method %in% c("FullSet", "PerfectMatch")))
summary(Model)


# PLOT LEN CUTOFFS, KEEP N_SNP=10
PLOT_65 <- ANALYSIS_STEP2_CUTOFFS %>%
  filter(N_SNP_Cutoff==10) %>%
  ggplot(aes(x=LEN_Cutoff, y=PPV, color=Matching_Method, linetype=Matching_Type)) +
  geom_line(linewidth=1) +
  labs(x="CNV Lenght Cutoff (bp)",
       y="PPV",
       linetype="CNV CALLSET QC",
       color="MATCHING METHOD",
       caption="N_SNP Cutoff = 10") +
  scale_linetype_manual(values=c("dashed", "solid"),
                        breaks=c("Raw", "QCd")) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRRmean", "LRRsd", "Pos", "PerfectMatch", "FullSet"),
                     labels=c("BAF", "LRRmean", "LRRsd", "Pos", "Perfect Match", "Full Set")) +
  scale_x_continuous(label=scales::comma) +
  scale_y_continuous(n.breaks=6) +
  theme_bw() + 
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.caption=element_text(size=12, face="bold.italic")) +
  guides(color=guide_legend(title.position="top", nrow=1, order=1),
         linetype=guide_legend(title.position="top", order=2))

ggsave(plot=PLOT_65,
       filename="FIGURES/Plot65.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white") 

# PLOT N_SNP CUTOFFS, KEEP LEN=20000
PLOT_66 <- ANALYSIS_STEP2_CUTOFFS %>%
  filter(LEN_Cutoff==20000) %>%
  ggplot(aes(x=N_SNP_Cutoff, y=PPV, color=Matching_Method, linetype=Matching_Type)) +
  geom_line(linewidth=1) +
  labs(x="N_SNP Cutoff",
       y="PPV",
       linetype="CNV CALLSET QC",
       color="MATCHING METHOD",
       caption="LEN Cutoff = 20,000") +
  scale_linetype_manual(values=c("dashed", "solid"),
                        breaks=c("Raw", "QCd")) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRRmean", "LRRsd", "Pos", "PerfectMatch", "FullSet"),
                     labels=c("BAF", "LRRmean", "LRRsd", "Pos", "Perfect Match", "Full Set")) +
  scale_x_continuous(label=scales::comma) +
  scale_y_continuous(n.breaks=6) +
  theme_bw() + 
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.caption=element_text(size=12, face="bold.italic")) +
  guides(color=guide_legend(title.position="top", nrow=1, order=1),
         linetype=guide_legend(title.position="top", order=2))

ggsave(plot=PLOT_66,
       filename="FIGURES/Plot66.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white") 

# PLOT LEN CUTOFFS, KEEP N_SNP=25
PLOT_67 <- ANALYSIS_STEP2_CUTOFFS %>%
  filter(N_SNP_Cutoff==25) %>%
  ggplot(aes(x=LEN_Cutoff, y=PPV, color=Matching_Method, linetype=Matching_Type)) +
  geom_line(linewidth=1) +
  labs(x="CNV Lenght Cutoff (bp)",
       y="PPV",
       linetype="CNV CALLSET QC",
       color="MATCHING METHOD",
       caption="N_SNP Cutoff = 25") +
  scale_linetype_manual(values=c("dashed", "solid"),
                        breaks=c("Raw", "QCd")) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRRmean", "LRRsd", "Pos", "PerfectMatch", "FullSet"),
                     labels=c("BAF", "LRRmean", "LRRsd", "Pos", "Perfect Match", "Full Set")) +
  scale_x_continuous(label=scales::comma) +
  scale_y_continuous(n.breaks=6) +
  theme_bw() + 
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.caption=element_text(size=12, face="bold.italic")) +
  guides(color=guide_legend(title.position="top", nrow=1, order=1),
         linetype=guide_legend(title.position="top", order=2))

ggsave(plot=PLOT_67,
       filename="FIGURES/Plot67.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white") 

# WRITE TABLE DATA
ANALYSIS_STEP2_CUTOFFS %>%
  rename(Method=Matching_Method,
         MaxD=Matching_Distance,
         QC=Matching_Type) %>%
  select(-MaxD_LOG) %>%
  write_tsv("TABLES/Table13.tsv", col_names=TRUE)
