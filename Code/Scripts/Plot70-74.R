# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)
library(ggpubr)

# CHECKPOINT
load("Analysis_SSC_111924.RData")

# TIDY UP ANALYSIS
ANALYSIS_SSC <- ANALYSIS_SSC %>%
  mutate(Sensitivity=round(TP/(TP+FN), digits=3),
         PPV=round(TP/(FP+TP), digits=3),
         FNR=round(FN/(FN+TP), digits=3),
         FDR=round(FP/(FP+TP), digits=3),
         F1=round((2*TP)/(2*TP+FP+FN), digits=3),
         FMI=round(sqrt((TP/(FP+TP))*(TP/(TP+FN))), digits=3),
         JI=round(TP/(TP+FN+FP), digits=3)) %>%
  mutate(D_MAX=as.numeric(D_MAX)) %>%
  mutate(D_MAX_LOG=log10(D_MAX)) %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_)) %>%
  mutate(QC=ifelse(Matching_Type=="Raw", "Low-stringency QC", "Medium-stringency QC"))

# KEEP ONLY MEDIUM-STRINGENCY QC
ANALYSIS_SSC <- ANALYSIS_SSC %>%
  filter(QC=="Medium-stringency QC") %>%
  select(-c(Matching_Type, Factor)) %>%
  filter(FactorN!="Perfect Match") %>%
  filter(!CNV_Size %in% c("Very Large", "Ultra Large")) %>%
  filter(D_MAX >= 100 | D_MAX==0)

# SCALE THE METRICS
ITER_DF <- ANALYSIS_SSC %>%
  select(c(D_MAX, FactorN, CNV_Type, CNV_Size)) %>%
  distinct() %>%
  filter(FactorN!="Full Set")

ANALYSIS_SSC_QC <- data.frame(NULL)

for(i in 1:nrow(ITER_DF)){
  ROW <- data.frame(NULL)
  
  ROW[1, "D_MAX"] <- ITER_DF$D_MAX[i]
  ROW[1, "Factor"] <- ITER_DF$FactorN[i]
  ROW[1, "CNV_Type"] <- ITER_DF$CNV_Type[i]
  ROW[1, "CNV_Size"] <- ITER_DF$CNV_Size[i]
  
  SCA <- ANALYSIS_SSC %>%
    filter(FactorN==ITER_DF$FactorN[i] & D_MAX==ITER_DF$D_MAX[i] & CNV_Type==ITER_DF$CNV_Type[i] & CNV_Size==ITER_DF$CNV_Size[i])
  
  REF <- ANALYSIS_SSC %>%
    filter(FactorN=="Full Set" & CNV_Type==ITER_DF$CNV_Type[i] & CNV_Size==ITER_DF$CNV_Size[i])
  
  ROW[1, "Sensitivity_p"] <- SCA$Sensitivity/REF$Sensitivity
  ROW[1, "PPV_p"] <- SCA$PPV/REF$PPV
  ROW[1, "F1_p"] <- SCA$F1/REF$F1
  ROW[1, "FMI_p"] <- SCA$FMI/REF$FMI
  ROW[1, "JI_p"] <- SCA$JI/REF$FMI
  
  ANALYSIS_SSC_QC <- rbind(ANALYSIS_SSC_QC, ROW)
  
  rm(ROW, SCA, REF)
}

ANALYSIS_SSC_QC$CNV_SizeF <- factor(ANALYSIS_SSC_QC$CNV_Size,
                                    levels = c("All", "Small", "Medium", "Large"),
                                    labels = c("All", "CNV < 100kb", "100kb < CNV < 500kb", "500kb < CNV < 1Mb"))

# PLOT THE SENSITIVITY OUTCOMES (PLOT 70)
(ggplot(ANALYSIS_SSC_QC, aes(x=log10(D_MAX), y=Sensitivity_p, color=Factor, linetype=CNV_Type)) +
  annotate("rect", xmin=3, xmax=5, ymin=-Inf, ymax=Inf, fill="steelblue3", alpha=0.1) +
  geom_vline(xintercept=4, color="steelblue3", linewidth=1.5) +
  geom_smooth(linewidth=1, method="loess", se=FALSE) +
  theme_bw() +
  facet_wrap(. ~ CNV_SizeF, nrow=4) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="SENSITIVITY '",
       linetype="CNV TYPE",
       color="FACTOR") +
  guides(color=guide_legend(title.position="top", nrow=1),
         linetype=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot70.png",
         device="png",
         width=7,
         height=11,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE PPV OUTCOMES (PLOT 71)
(ggplot(ANALYSIS_SSC_QC, aes(x=log10(D_MAX), y=PPV_p, color=Factor, linetype=CNV_Type)) +
  annotate("rect", xmin=3, xmax=5, ymin=-Inf, ymax=Inf, fill="steelblue3", alpha=0.1) +
  geom_vline(xintercept=4, color="steelblue3", linewidth=1.5) +
  geom_smooth(linewidth=1, method="loess", se=FALSE) +
  theme_bw() +
  facet_wrap(. ~ CNV_SizeF, nrow=4) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="PPV '",
       linetype="CNV TYPE",
       color="FACTOR") +
  guides(color=guide_legend(title.position="top", nrow=1),
         linetype=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot71.png",
         device="png",
         width=7,
         height=11,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE F1 OUTCOMES  (PLOT 72)
(ggplot(ANALYSIS_SSC_QC, aes(x=log10(D_MAX), y=F1_p, color=Factor, linetype=CNV_Type)) +
  annotate("rect", xmin=3, xmax=5, ymin=-Inf, ymax=Inf, fill="steelblue3", alpha=0.1) +
  geom_vline(xintercept=4, color="steelblue3", linewidth=1.5) +
  geom_smooth(linewidth=1, method="loess", se=FALSE) +
  theme_bw() +
  facet_wrap(. ~ CNV_SizeF, nrow=4) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="F1 '",
       linetype="CNV TYPE",
       color="FACTOR") +
  guides(color=guide_legend(title.position="top", nrow=1),
         linetype=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot72.png",
         device="png",
         width=7,
         height=11,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE FMI OUTCOMES  (PLOT 73)
(ggplot(ANALYSIS_SSC_QC, aes(x=log10(D_MAX), y=FMI_p, color=Factor, linetype=CNV_Type)) +
  annotate("rect", xmin=3, xmax=5, ymin=-Inf, ymax=Inf, fill="steelblue3", alpha=0.1) +
  geom_vline(xintercept=4, color="steelblue3", linewidth=1.5) +
  geom_smooth(linewidth=1, method="loess", se=FALSE) +
  theme_bw() +
  facet_wrap(. ~ CNV_SizeF, nrow=4) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="FMI '",
       linetype="CNV TYPE",
       color="FACTOR") +
  guides(color=guide_legend(title.position="top", nrow=1),
         linetype=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot73.png",
         device="png",
         width=7,
         height=11,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE JI OUTCOMES  (PLOT 74)
(ggplot(ANALYSIS_SSC_QC, aes(x=log10(D_MAX), y=JI_p, color=Factor, linetype=CNV_Type)) +
  annotate("rect", xmin=3, xmax=5, ymin=-Inf, ymax=Inf, fill="steelblue3", alpha=0.1) +
  geom_vline(xintercept=4, color="steelblue3", linewidth=1.5) +
  geom_smooth(linewidth=1, method="loess", se=FALSE) +
  theme_bw() +
  facet_wrap(. ~ CNV_SizeF, nrow=4) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="JI '",
       linetype="CNV TYPE",
       color="FACTOR") +
  guides(color=guide_legend(title.position="top", nrow=1),
         linetype=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot74.png",
         device="png",
         width=7,
         height=11,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE SENSITIVITY FACTOR OUTCOMES  (PLOT 75)
(ggplot(ANALYSIS_SSC_QC, aes(y=Sensitivity_p, color=Factor, fill=Factor)) +
  geom_boxplot(alpha=0.2) +
  theme_bw() +
  facet_wrap(CNV_Type ~ CNV_SizeF, nrow=3) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
  scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
        strip.background=element_rect(fill="#ffffff"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        axis.ticks.x=element_blank()) +
  labs(x=expression(bold("FACTOR")),
       y="SENSITIVITY '",
       fill="FACTOR",
       color="FACTOR") +
  guides(color=guide_legend(title.position="top", nrow=1),
         fill=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot75.png",
         device="png",
         width=11,
         height=9,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE PPV FACTOR OUTCOMES  (PLOT 76)
(ggplot(ANALYSIS_SSC_QC, aes(y=PPV_p, color=Factor, fill=Factor)) +
    geom_boxplot(alpha=0.2) +
    theme_bw() +
    facet_wrap(CNV_Type ~ CNV_SizeF, nrow=3) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                      breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold"),
          axis.ticks.x=element_blank()) +
    labs(x=expression(bold("FACTOR")),
         y="PPV '",
         fill="FACTOR",
         color="FACTOR") +
    guides(color=guide_legend(title.position="top", nrow=1),
           fill=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot76.png",
         device="png",
         width=11,
         height=9,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE F1 FACTOR OUTCOMES  (PLOT 77)
(ggplot(ANALYSIS_SSC_QC, aes(y=F1_p, color=Factor, fill=Factor)) +
    geom_boxplot(alpha=0.2) +
    theme_bw() +
    facet_wrap(CNV_Type ~ CNV_SizeF, nrow=3) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                      breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold"),
          axis.ticks.x=element_blank()) +
    labs(x=expression(bold("FACTOR")),
         y="F1 '",
         fill="FACTOR",
         color="FACTOR") +
    guides(color=guide_legend(title.position="top", nrow=1),
           fill=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot77.png",
         device="png",
         width=11,
         height=9,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE FMI FACTOR OUTCOMES  (PLOT 78)
(ggplot(ANALYSIS_SSC_QC, aes(y=FMI_p, color=Factor, fill=Factor)) +
    geom_boxplot(alpha=0.2) +
    theme_bw() +
    facet_wrap(CNV_Type ~ CNV_SizeF, nrow=3) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                      breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold"),
          axis.ticks.x=element_blank()) +
    labs(x=expression(bold("FACTOR")),
         y="FMI '",
         fill="FACTOR",
         color="FACTOR") +
    guides(color=guide_legend(title.position="top", nrow=1),
           fill=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot78.png",
         device="png",
         width=11,
         height=9,
         units="in",
         dpi=350,
         bg="white")

# PLOT THE JI FACTOR OUTCOMES  (PLOT 79)
(ggplot(ANALYSIS_SSC_QC, aes(y=JI_p, color=Factor, fill=Factor)) +
    geom_boxplot(alpha=0.2) +
    theme_bw() +
    facet_wrap(CNV_Type ~ CNV_SizeF, nrow=3) +
    scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                       breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                      breaks=c("BAF", "LRR mean", "LRR sd", "Distance")) +
    scale_y_continuous(breaks=scales::pretty_breaks(n = 3)) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
          axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
          axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
          strip.text=element_text(size=15, hjust=0, vjust=0.5, face="bold.italic"),
          strip.background=element_rect(fill="#ffffff"),
          legend.position="top",
          legend.justification="left",
          legend.title=element_text(size=12, face="bold"),
          axis.ticks.x=element_blank()) +
    labs(x=expression(bold("FACTOR")),
         y="JI '",
         fill="FACTOR",
         color="FACTOR") +
    guides(color=guide_legend(title.position="top", nrow=1),
           fill=guide_legend(title.position="top"))) %>%
  ggsave(filename="FIGURES/Plot79.png",
         device="png",
         width=11,
         height=9,
         units="in",
         dpi=350,
         bg="white")

# RUN A T-TEST ON THESE DATA
TEST_FACTOR <- data.frame(NULL)

for(i in unique(ANALYSIS_SSC_QC$CNV_Type)){
  for(j in unique(ANALYSIS_SSC_QC$CNV_Size)){
    for(m in c("Sensitivity_p", "PPV_p", "F1_p", "FMI_p", "JI_p")){
      for(k in unique(ANALYSIS_SSC_QC$Factor)){
      
        ROW <- data.frame(NULL)
        ROW[1, "CNV_Type"] <- i
        ROW[1, "CNV_Size"] <- j
        ROW[1, "Metric"] <- m
        ROW[1, "Factor"] <- k
        
        for(l in unique(ANALYSIS_SSC_QC$Factor)){
          
            A <- ANALYSIS_SSC_QC %>%
              filter(CNV_Type==i, CNV_Size==j, Factor==k) %>%
              pull(.data[[m]]) %>% unlist()
            
            B <- ANALYSIS_SSC_QC %>%
              filter(CNV_Type==i, CNV_Size==j, Factor==l)%>%
              pull(.data[[m]]) %>% unlist()
  
            ROW[1, l] <- t.test(A, B)$p.value
        }
        
        TEST_FACTOR <- rbind(TEST_FACTOR, ROW)
        
      }
    }
  }
}

TEST_FACTOR <- TEST_FACTOR %>%
  mutate(BAF=ifelse(Factor %in% c("BAF"), NA, BAF),
         `LRR mean`=ifelse(Factor %in% c("BAF", "LRR mean"), NA, `LRR mean`),
         `LRR sd`=ifelse(Factor %in% c("BAF", "LRR mean", "LRR sd"), NA, `LRR sd`),
         Distance=NA) %>%
  pivot_longer(c(BAF:Distance), names_to="Factor2", values_to="P") %>%
  drop_na(P) %>%
  mutate(Q=p.adjust(P, method="fdr")) %>% 
  pivot_wider(names_from="Factor2",
              values_from=c("P", "P.adj"),
              names_vary="slowest") %>%
  write_tsv("TABLES/TableS1K.tsv", col_names=TRUE)
