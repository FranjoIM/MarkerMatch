# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)
library(betareg)
library(easystats)
library(marginaleffects)

# LOAD DATA
load("Validation_Final_103024.RData")

DataFileNames <- data.frame(
  Factor=rep(c("FullSet", "PerfectMatch", "BAF", "LRRmean", "LRRsd", "Pos"), 2),
  D_MAX=rep(c(NA, NA, rep("10000", 4)), 2),
  D_MAXLab=rep(c("0", "0", rep(c("10000"), 4)), 2),
  Ref=c(rep("GSA", 6), rep("OEE", 6)),
  Mat=c(rep("OEE", 6), rep("GSA", 6)),
  stringsAsFactors=FALSE)

CNV_Lengths <- seq(from=20000, to=500000, by=20000)
SNP_Numbers <- seq(from=10, to=50, by=5)
CNV_Types <- c("Raw", "QCd")

ANALYSIS_STEP2_CUTOFFS <- data.frame(NULL)

for(h in 1:nrow(DataFileNames)){
  i <- DataFileNames$Factor[h]
  j <- DataFileNames$D_MAXLab[h]
  k <- DataFileNames$Ref[h]
  l <- DataFileNames$Mat[h]
  
  for(m in CNV_Lengths){
    for(n in SNP_Numbers){
      for(o in CNV_Types){
        
        # Pull out CNVs overlapping in full set and partial sets
        REF <- DATA[["Raw"]][[k]][["FullSet"]][["0"]][["CNV"]] %>%
          filter(LEN>=m & N_SNP>=n)
        MAT <- DATA[[o]][[l]][[i]][[j]][["CNV"]] %>%
          filter(LEN>=m & N_SNP>=n)
        
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
          Factor=i,
          D_MAX=j,
          QC=o,
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
  mutate(D_MAX=as.numeric(D_MAX)) %>%
  mutate(D_MAX_LOG=log10(D_MAX)) %>%
  mutate(QC=ifelse(QC=="Raw", "Low-stringency QC", "Medium-stringency QC")) %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))

ANALYSIS_STEP2_CUTOFFS$FactorF <- factor(ANALYSIS_STEP2_CUTOFFS$FactorN,
    labels=c("Full Set", "Perfect Match", "LRR mean", "LRR sd", "BAF", "Distance"),
    levels=c("Full Set", "Perfect Match", "LRR mean", "LRR sd", "BAF", "Distance"))

# TRANSFORM PRIMES TO RANGE (0,1) from [0,1]
SV_Trans <- function(x){
  (x * ((nrow(ANALYSIS_STEP2_CUTOFFS)/2) - 1 ) + 0.5)/(nrow(ANALYSIS_STEP2_CUTOFFS)/2)
}

ANALYSIS_STEP2_CUTOFFS <- ANALYSIS_STEP2_CUTOFFS %>%
  mutate(across(c(Sensitivity, PPV, F1, FMI, JI), SV_Trans))

ANALYSIS_S2VS_QC <- ANALYSIS_STEP2_CUTOFFS %>%
  filter(QC=="Medium-stringency QC") %>%
  mutate(LEN_Cutoff=LEN_Cutoff/10000)

# MODEL METRICS ON LEN AND N_SNP CUTOFFS
MOD_SEN_GSA <- betareg(Sensitivity ~ 
                 LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                 data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="GSA"),
                 link="logit")

summary(MOD_SEN_GSA)

as.data.frame(model_parameters(MOD_SEN_GSA, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_SEN_OEE <- betareg(Sensitivity ~ 
                         LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                       data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="OEE"),
                       link="logit")

summary(MOD_SEN_OEE)

as.data.frame(model_parameters(MOD_SEN_OEE, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_PPV_GSA <- betareg(PPV ~ 
                         LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                       data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="GSA"),
                       link="logit")

summary(MOD_PPV_GSA)

as.data.frame(model_parameters(MOD_PPV_GSA, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_PPV_OEE <- betareg(PPV ~ 
                         LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                       data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="OEE"),
                       link="logit")

summary(MOD_PPV_OEE)

as.data.frame(model_parameters(MOD_PPV_OEE, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_F1_GSA <- betareg(F1 ~ 
                         LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                       data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="GSA"),
                       link="logit")

summary(MOD_F1_GSA)

as.data.frame(model_parameters(MOD_F1_GSA, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_F1_OEE <- betareg(F1 ~ 
                         LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                       data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="OEE"),
                       link="logit")

summary(MOD_F1_OEE)

as.data.frame(model_parameters(MOD_F1_OEE, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_FMI_GSA <- betareg(FMI ~ 
                        LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                      data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="GSA"),
                      link="logit")

summary(MOD_FMI_GSA)

as.data.frame(model_parameters(MOD_FMI_GSA, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_FMI_OEE <- betareg(FMI ~ 
                        LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                      data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="OEE"),
                      link="logit")

summary(MOD_FMI_OEE)

as.data.frame(model_parameters(MOD_FMI_OEE, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_JI_GSA <- betareg(JI ~ 
                        LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                      data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="GSA"),
                      link="logit")

summary(MOD_JI_GSA)

as.data.frame(model_parameters(MOD_JI_GSA, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

MOD_JI_OEE <- betareg(JI ~ 
                        LEN_Cutoff * N_SNP_Cutoff + FactorF | LEN_Cutoff + N_SNP_Cutoff + FactorF,
                      data=filter(ANALYSIS_S2VS_QC, !FactorN %in% c("Full Set", "Perfect Match") & Mat=="OEE"),
                      link="logit")

summary(MOD_JI_OEE)

as.data.frame(model_parameters(MOD_JI_OEE, component="conditional", exponentiate=TRUE)) %>%
  select(-c(CI, CI_low, CI_high, df_error)) %>%
  print(row.names=FALSE)

(ANALYSIS_S2VS_QC %>%
    filter(N_SNP_Cutoff!=50) %>%
    mutate(N_SNP_Cutoff=paste0("Cutoff = ", N_SNP_Cutoff)) %>%
    ggplot(aes(x=LEN_Cutoff*10, y=PPV, color=FactorF, linetype=Mat)) +
    geom_line(linewidth=1) +
    labs(x="CNV LENGTH CUTOFF (KB)",
         y="PPV",
         linetype="ARRAY",
         color="FACTOR",
         caption="Cutoff: marker coverage") +
    scale_linetype_manual(values=c("dashed", "solid"),
                          breaks=c("OEE", "GSA")) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set"),
                       labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
    scale_x_continuous(label=scales::comma, n.breaks=5) +
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
          plot.caption=element_text(size=12, face="bold.italic"),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),) +
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           linetype=guide_legend(title.position="top", order=2)) +
  facet_wrap(. ~ N_SNP_Cutoff, nrow=2)) %>%
  ggsave(filename="FIGURES/Plot64.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")
  
(ANALYSIS_S2VS_QC %>%
    filter(N_SNP_Cutoff!=50) %>%
    mutate(N_SNP_Cutoff=paste0("Cutoff = ", N_SNP_Cutoff)) %>%
    ggplot(aes(x=LEN_Cutoff*10, y=F1, color=FactorF, linetype=Mat)) +
    geom_line(linewidth=1) +
    labs(x="CNV LENGTH CUTOFF (KB)",
         y="F1",
         linetype="ARRAY",
         color="FACTOR",
         caption="Cutoff: marker coverage") +
    scale_linetype_manual(values=c("dashed", "solid"),
                          breaks=c("OEE", "GSA")) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set"),
                       labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
    scale_x_continuous(label=scales::comma, n.breaks=5) +
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
          plot.caption=element_text(size=12, face="bold.italic"),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),) +
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           linetype=guide_legend(title.position="top", order=2)) +
    facet_wrap(. ~ N_SNP_Cutoff, nrow=2)) %>%
    ggsave(filename="FIGURES/Plot65.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

(ANALYSIS_S2VS_QC %>%
    filter(LEN_Cutoff %in% c(2, 4, 6, 8, 10, 20, 30, 40)) %>%
    mutate(LEN_Cutoff=paste0("Cutoff = ", LEN_Cutoff, "0,000bp")) %>%
    ggplot(aes(x=N_SNP_Cutoff, y=PPV, color=FactorF, linetype=Mat)) +
    geom_line(linewidth=1) +
    labs(x="MARKER COVERAGE CUTOFF",
         y="PPV",
         linetype="ARRAY",
         color="FACTOR",
         caption="Cutoff: CNV length") +
    scale_linetype_manual(values=c("dashed", "solid"),
                          breaks=c("OEE", "GSA")) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set"),
                       labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
    scale_x_continuous(label=scales::comma, n.breaks=5) +
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
          plot.caption=element_text(size=12, face="bold.italic"),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),) +
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           linetype=guide_legend(title.position="top", order=2)) +
    facet_wrap(. ~ factor(LEN_Cutoff, levels=c("Cutoff = 20,000bp", "Cutoff = 40,000bp", "Cutoff = 60,000bp", "Cutoff = 80,000bp", "Cutoff = 100,000bp", "Cutoff = 200,000bp", "Cutoff = 300,000bp", "Cutoff = 400,000bp")), nrow=2)) %>%
    ggsave(filename="FIGURES/Plot66.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

(ANALYSIS_S2VS_QC %>%
    filter(LEN_Cutoff %in% c(2, 4, 6, 8, 10, 20, 30, 40)) %>%
    mutate(LEN_Cutoff=paste0("Cutoff = ", LEN_Cutoff, "0,000bp")) %>%
    ggplot(aes(x=N_SNP_Cutoff, y=F1, color=FactorF, linetype=Mat)) +
    geom_line(linewidth=1) +
    labs(x="MARKER COVERAGE CUTOFF",
         y="F1",
         linetype="ARRAY",
         color="FACTOR",
         caption="Cutoff: CNV length") +
    scale_linetype_manual(values=c("dashed", "solid"),
                          breaks=c("OEE", "GSA")) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set"),
                       labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
    scale_x_continuous(label=scales::comma, n.breaks=5) +
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
          plot.caption=element_text(size=12, face="bold.italic"),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),) +
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           linetype=guide_legend(title.position="top", order=2)) +
    facet_wrap(. ~ factor(LEN_Cutoff, levels=c("Cutoff = 20,000bp", "Cutoff = 40,000bp", "Cutoff = 60,000bp", "Cutoff = 80,000bp", "Cutoff = 100,000bp", "Cutoff = 200,000bp", "Cutoff = 300,000bp", "Cutoff = 400,000bp")), nrow=2)) %>%
  ggsave(filename="FIGURES/Plot67.png",
         device="png",
         width=11,
         height=7,
         units="in",
         dpi=350,
         bg="white")

ANALYSIS_S2VS_QC$FactorF <- factor(ANALYSIS_S2VS_QC$FactorN,
                                  labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Full Set", "Perfect Match"),
                                  levels=c("BAF", "LRR mean", "LRR sd", "Distance", "Full Set", "Perfect Match"))

(ANALYSIS_S2VS_QC %>%
  filter(N_SNP_Cutoff %in% c(10, 20, 30, 40) & LEN_Cutoff %in% c(2, 4, 8, 16, 32)) %>%
  mutate(N_SNP_Cutoff=paste0(N_SNP_Cutoff, " SNPs"),
         LEN_Cutoff=paste0(LEN_Cutoff, "0,000bp")) %>%
  ggplot(aes(y=PPV, x=FactorF, color=FactorF, shape=Mat)) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set"),
                     labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_shape_manual(values=c(2, 6), breaks=c("GSA", "OEE")) +
  geom_point(size=3) +
  scale_y_continuous(n.breaks=4, limits=c(0, 1)) +
  labs(x="FACTOR",
       y="PPV",
       shape="ARRAY",
       color="FACTOR") +
  theme_bw() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.caption=element_text(size=12, face="bold.italic"),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff")) +
  guides(color=guide_legend(title.position="top", nrow=1, order=1),
         shape=guide_legend(title.position="top", order=2)) +
  facet_grid(factor(LEN_Cutoff, levels=c("20,000bp", "40,000bp", "80,000bp", "160,000bp", "320,000bp")) ~ N_SNP_Cutoff)) %>%
  ggsave(filename="FIGURES/Plot68.png",
         device="png",
         width=9,
         height=10,
         units="in",
         dpi=350,
         bg="white")

(ANALYSIS_S2VS_QC %>%
    filter(N_SNP_Cutoff %in% c(10, 20, 30, 40) & LEN_Cutoff %in% c(2, 4, 8, 16, 32)) %>%
    mutate(N_SNP_Cutoff=paste0(N_SNP_Cutoff, " SNPs"),
           LEN_Cutoff=paste0(LEN_Cutoff, "0,000bp")) %>%
    ggplot(aes(y=F1, x=FactorF, color=FactorF, shape=Mat)) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set"),
                       labels=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
    scale_shape_manual(values=c(2, 6), breaks=c("GSA", "OEE")) +
    geom_point(size=3) +
    scale_y_continuous(n.breaks=4, limits=c(0, 1)) +
    labs(x="FACTOR",
         y="F1",
         shape="ARRAY",
         color="FACTOR") +
    theme_bw() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold"),
          plot.caption=element_text(size=12, face="bold.italic"),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff")) +
    guides(color=guide_legend(title.position="top", nrow=1, order=1),
           shape=guide_legend(title.position="top", order=2)) +
    facet_grid(factor(LEN_Cutoff, levels=c("20,000bp", "40,000bp", "80,000bp", "160,000bp", "320,000bp")) ~ N_SNP_Cutoff)) %>%
  ggsave(filename="FIGURES/Plot69.png",
         device="png",
         width=9,
         height=10,
         units="in",
         dpi=350,
         bg="white")
