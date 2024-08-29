# Set working directory
W_DIR <- "..."
setwd(W_DIR)

# Load packages
library(tidyverse)
library(ggpubr)

# Import array manifests
OMNI_MAN <- read_delim("DATA/OMNI_MAN.txt", delim="\t")
GSA_MAN <- read_delim("DATA/GSA_MAN.txt", delim="\t")
OEE_MAN <- read_delim("DATA/OEE_MAN.txt", delim="\t")
OMNI_SNPS <- read_delim("DATA/UsedSNPs_OMNI.txt", delim="\t", col_names=T)
GSA_SNPS <- read_delim("DATA/UsedSNPs_GSA.txt", delim="\t", col_names=T)
OEE_SNPS <- read_delim("DATA/UsedSNPs_OEE.txt", delim="\t", col_names=T)

### LRRMEANS PLOTS FOR OMNI-GSA
LRRMEANS_OMNI <- data.frame()
SETS_OMNI <- colnames(OMNI_SNPS[2:51])

for(j in SETS_OMNI){
  FNAME <- OMNI_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
  LRRMEAN <- OMNI_MAN %>%
    filter(Name %in% FNAME) %>%
    pull(LRR_mean)
  
  a <- str_split(j, pattern="_")[[1]][1]
  b <- str_split(j, pattern="_")[[1]][2]
  
  ROW <- data.frame(Factor=a,
                    MaxD=b,
                    LRRMEAN=LRRMEAN)
  
  LRRMEANS_OMNI <- rbind(LRRMEANS_OMNI, ROW)
}

LRRMEANS_OMNI$MaxDF <- factor(LRRMEANS_OMNI$MaxD, 
                            levels=c("0", "10", "50", "100", "500", "1000", "5000", "10000", 
                                     "50000", "100000", "500000", "1000000", "5000000"),
                            labels=c("0", "10", "50", "100", "500", "1,000", "5,000", "10,000", 
                                     "50,000", "100,000", "500,000", "1,000,000", "5,000,000"))

LRRMEANS_OMNI <- LRRMEANS_OMNI %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_)) %>%
  drop_na(LRRMEAN)

LRRMEANS_OMNI$FactorF <- factor(LRRMEANS_OMNI$FactorN,
                              levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                       "LRR sd", "Distance"),
                              labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                       "LRR sd", "Distance"))

H1_OMNI <- LRRMEANS_OMNI %>%
  filter(Factor=="PerfectMatch") %>%
  pull(LRRMEAN) %>%
  median(., na.rm=TRUE)

H2_OMNI <- LRRMEANS_OMNI %>%
  filter(Factor=="FullSet") %>%
  pull(LRRMEAN) %>%
  median(., na.rm=TRUE)

Plot7_A <- LRRMEANS_OMNI %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=LRRMEAN)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_OMNI, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H2_OMNI, color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ LRR-mean ]")),
       subtitle="OMNI matched on GSA",
       color="MATCHING METHOD",
       fill="MATCHING METHOD") +
  theme_bw() + 
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(face="bold.italic", size=12),
        strip.background=element_rect(fill="transparent"),
        axis.ticks.x=element_blank()) +
  guides(color=guide_legend(title.position="top", nrow=1),
         fill=guide_legend(title.position="top", nrow=1)) + 
  facet_grid(. ~ MaxDF, scales="free_x", space="free_x", switch="both")

ggsave(plot=Plot7_A,
       filename="FIGURES/Plot7_A.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

ggsave(plot=Plot7_A + coord_cartesian(ylim=c(-0.02, 0.02)),
       filename="FIGURES/Plot7_A_ZOOM.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### LRRMEANS PLOTS FOR OEE-GSA
LRRMEANS_OEE <- data.frame()
SETS_OEE <- colnames(OEE_SNPS[2:7])

for(j in SETS_OEE){
  FNAME <- OEE_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
  LRRMEAN <- OEE_MAN %>%
    filter(Name %in% FNAME) %>%
    pull(LRR_mean)
  
  a <- str_split(j, pattern="_")[[1]][1]
  b <- str_split(j, pattern="_")[[1]][2]
  
  ROW <- data.frame(Factor=a,
                    MaxD=b,
                    LRRMEAN=LRRMEAN)
  
  LRRMEANS_OEE <- rbind(LRRMEANS_OEE, ROW)
}

LRRMEANS_OEE$MaxDF <- factor(LRRMEANS_OEE$MaxD, 
                           levels=c("0", "10000"),
                           labels=c("0", "10,000"))

LRRMEANS_OEE <- LRRMEANS_OEE %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_)) %>%
  drop_na(LRRMEAN)

LRRMEANS_OEE$FactorF <- factor(LRRMEANS_OEE$FactorN,
                             levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                      "LRR sd", "Distance"),
                             labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                      "LRR sd", "Distance"))

H1_OEE <- LRRMEANS_OEE %>%
  filter(Factor=="PerfectMatch") %>%
  pull(LRRMEAN) %>%
  median(., na.rm=TRUE)

H2_OEE <- LRRMEANS_OEE %>%
  filter(Factor=="FullSet") %>%
  pull(LRRMEAN) %>%
  median(., na.rm=TRUE)

Plot7_B <- LRRMEANS_OEE %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=LRRMEAN)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_OEE, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H2_OEE, color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ LRR-mean ]")),
       subtitle="OEE matched on GSA",
       color="MATCHING METHOD",
       fill="MATCHING METHOD") +
  theme_bw() + 
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(face="bold.italic", size=12),
        strip.background=element_rect(fill="transparent"),
        axis.ticks.x=element_blank()) +
  guides(color=guide_legend(title.position="top", nrow=1),
         fill=guide_legend(title.position="top", nrow=1)) + 
  facet_grid(. ~ MaxDF, scales="free_x", space="free_x", switch="both")

ggsave(plot=Plot7_B,
       filename="FIGURES/Plot7_B.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

ggsave(plot=Plot7_B + coord_cartesian(ylim=c(-0.02, 0.02)),
       filename="FIGURES/Plot7_B_ZOOM.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### LRRMEANS PLOTS FOR GSA-OEE
LRRMEANS_GSA <- data.frame()
SETS_GSA <- colnames(GSA_SNPS[2:7])

for(j in SETS_GSA){
  FNAME <- GSA_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
  LRRMEAN <- GSA_MAN %>%
    filter(Name %in% FNAME) %>%
    pull(LRR_mean)
  
  a <- str_split(j, pattern="_")[[1]][1]
  b <- str_split(j, pattern="_")[[1]][2]
  
  ROW <- data.frame(Factor=a,
                    MaxD=b,
                    LRRMEAN=LRRMEAN)
  
  LRRMEANS_GSA <- rbind(LRRMEANS_GSA, ROW)
}

LRRMEANS_GSA$MaxDF <- factor(LRRMEANS_GSA$MaxD, 
                           levels=c("0", "10000"),
                           labels=c("0", "10,000"))

LRRMEANS_GSA <- LRRMEANS_GSA %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_)) %>%
  drop_na(LRRMEAN)

LRRMEANS_GSA$FactorF <- factor(LRRMEANS_GSA$FactorN,
                             levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                      "LRR sd", "Distance"),
                             labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                      "LRR sd", "Distance"))

H1_GSA <- LRRMEANS_GSA %>%
  filter(Factor=="PerfectMatch") %>%
  pull(LRRMEAN) %>%
  median(., na.rm=TRUE)

H2_GSA <- LRRMEANS_GSA %>%
  filter(Factor=="FullSet") %>%
  pull(LRRMEAN) %>%
  median(., na.rm=TRUE)

Plot7_C <- LRRMEANS_GSA %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=LRRMEAN)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_GSA, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H2_GSA, color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ LRR-mean ]")),
       subtitle="GSA matched on OEE",
       color="MATCHING METHOD",
       fill="MATCHING METHOD") +
  theme_bw() + 
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(face="bold.italic", size=12),
        strip.background=element_rect(fill="transparent"),
        axis.ticks.x=element_blank()) +
  guides(color=guide_legend(title.position="top", nrow=1),
         fill=guide_legend(title.position="top", nrow=1)) + 
  facet_grid(. ~ MaxDF, scales="free_x", space="free_x", switch="both")

ggsave(plot=Plot7_C,
       filename="FIGURES/Plot7_C.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

ggsave(plot=Plot7_C + coord_cartesian(ylim=c(-0.02, 0.02)),
       filename="FIGURES/Plot7_C_ZOOM.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# SAVE ALL LRRMEANS FILES AS RDATA
save(LRRMEANS_OMNI, LRRMEANS_OEE, LRRMEANS_GSA, file="Table8.RData")

# TABULATE GAP SIZES ACROSS TABLES AND CONDITIONS
LRRMEANS_ALL <- bind_rows(
  LRRMEANS_GSA %>%
    select(c(FactorN, MaxD, LRRMEAN)) %>%
    mutate(Matching="GSA", Reference="OEE"),
  LRRMEANS_OEE %>%
    select(c(FactorN, MaxD, LRRMEAN)) %>%
    mutate(Matching="OEE", Reference="GSA"),
  LRRMEANS_OMNI %>%
    select(c(FactorN, MaxD, LRRMEAN)) %>%
    mutate(Matching="OMNI", Reference="GSA"))

LRRMEANS_ALL %>%
  group_by(Matching, Reference, FactorN, MaxD) %>%
  summarise(Min = min(LRRMEAN),
            Med = median(LRRMEAN),
            Mean = mean(LRRMEAN),
            Max = max(LRRMEAN),
            StDev = sd(LRRMEAN),
            IQR = IQR(LRRMEAN),
            .groups="keep") %>%
  write_tsv("TABLES/Table8.tsv", col_names=TRUE)
