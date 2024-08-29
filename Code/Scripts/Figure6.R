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

### LRRSDS PLOTS FOR OMNI-GSA
LRRSDS_OMNI <- data.frame()
SETS_OMNI <- colnames(OMNI_SNPS[2:51])

for(j in SETS_OMNI){
  FNAME <- OMNI_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
    LRRSD <- OMNI_MAN %>%
      filter(Name %in% FNAME) %>%
      pull(LRR_sd)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      LRRSD=LRRSD)
    
    LRRSDS_OMNI <- rbind(LRRSDS_OMNI, ROW)
}

LRRSDS_OMNI$MaxDF <- factor(LRRSDS_OMNI$MaxD, 
                            levels=c("0", "10", "50", "100", "500", "1000", "5000", "10000", 
                                     "50000", "100000", "500000", "1000000", "5000000"),
                            labels=c("0", "10", "50", "100", "500", "1,000", "5,000", "10,000", 
                                     "50,000", "100,000", "500,000", "1,000,000", "5,000,000"))

LRRSDS_OMNI <- LRRSDS_OMNI %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_)) %>%
  drop_na(LRRSD)

LRRSDS_OMNI$FactorF <- factor(LRRSDS_OMNI$FactorN,
                         levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                  "LRR sd", "Distance"),
                         labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                  "LRR sd", "Distance"))

H1_OMNI <- LRRSDS_OMNI %>%
  filter(Factor=="PerfectMatch") %>%
  pull(LRRSD) %>%
  median(., na.rm=TRUE)

H2_OMNI <- LRRSDS_OMNI %>%
  filter(Factor=="FullSet") %>%
  pull(LRRSD) %>%
  median(., na.rm=TRUE)

Plot6_A <- LRRSDS_OMNI %>%
  filter(LRRSD>0) %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=log10(LRRSD))) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=log10(H1_OMNI), color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=log10(H2_OMNI), color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ LRR-sd ]")),
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

ggsave(plot=Plot6_A,
       filename="FIGURES/Plot6_A.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### LRRSDS PLOTS FOR OEE-GSA
LRRSDS_OEE <- data.frame()
SETS_OEE <- colnames(OEE_SNPS[2:7])

for(j in SETS_OEE){
  FNAME <- OEE_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
    LRRSD <- OEE_MAN %>%
      filter(Name %in% FNAME) %>%
      pull(LRR_sd)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      LRRSD=LRRSD)
    
    LRRSDS_OEE <- rbind(LRRSDS_OEE, ROW)
}

LRRSDS_OEE$MaxDF <- factor(LRRSDS_OEE$MaxD, 
                       levels=c("0", "10000"),
                       labels=c("0", "10,000"))

LRRSDS_OEE <- LRRSDS_OEE %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_)) %>%
  drop_na(LRRSD)

LRRSDS_OEE$FactorF <- factor(LRRSDS_OEE$FactorN,
                         levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                  "LRR sd", "Distance"),
                         labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                  "LRR sd", "Distance"))

H1_OEE <- LRRSDS_OEE %>%
  filter(Factor=="PerfectMatch") %>%
  pull(LRRSD) %>%
  median(., na.rm=TRUE)

H2_OEE <- LRRSDS_OEE %>%
  filter(Factor=="FullSet") %>%
  pull(LRRSD) %>%
  median(., na.rm=TRUE)

Plot6_B <- LRRSDS_OEE %>%
  filter(LRRSD>0) %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=log10(LRRSD))) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=log10(H1_OEE), color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=log10(H2_OEE), color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ LRR-sd ]")),
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

ggsave(plot=Plot6_B,
       filename="FIGURES/Plot6_B.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### LRRSDS PLOTS FOR GSA-OEE
LRRSDS_GSA <- data.frame()
SETS_GSA <- colnames(GSA_SNPS[2:7])

for(j in SETS_GSA){
  FNAME <- GSA_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
    LRRSD <- GSA_MAN %>%
      filter(Name %in% FNAME) %>%
      pull(LRR_sd)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      LRRSD=LRRSD)
    
    LRRSDS_GSA <- rbind(LRRSDS_GSA, ROW)
}

LRRSDS_GSA$MaxDF <- factor(LRRSDS_GSA$MaxD, 
                       levels=c("0", "10000"),
                       labels=c("0", "10,000"))

LRRSDS_GSA <- LRRSDS_GSA %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_)) %>%
  drop_na(LRRSD)

LRRSDS_GSA$FactorF <- factor(LRRSDS_GSA$FactorN,
                         levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                  "LRR sd", "Distance"),
                         labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                  "LRR sd", "Distance"))

H1_GSA <- LRRSDS_GSA %>%
  filter(Factor=="PerfectMatch") %>%
  pull(LRRSD) %>%
  median(., na.rm=TRUE)

H2_GSA <- LRRSDS_GSA %>%
  filter(Factor=="FullSet") %>%
  pull(LRRSD) %>%
  median(., na.rm=TRUE)

Plot6_C <- LRRSDS_GSA %>%
  filter(LRRSD>0) %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=log10(LRRSD))) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=log10(H1_GSA), color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=log10(H2_GSA), color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ LRR-sd ]")),
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

ggsave(plot=Plot6_C,
       filename="FIGURES/Plot6_C.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# SAVE ALL LRRSDS FILES AS RDATA
save(LRRSDS_OMNI, LRRSDS_OEE, LRRSDS_GSA, file="Table7.RData")

# TABULATE GAP SIZES ACROSS TABLES AND CONDITIONS
LRRSDS_ALL <- bind_rows(
  LRRSDS_GSA %>%
    select(c(FactorN, MaxD, LRRSD)) %>%
    mutate(Matching="GSA", Reference="OEE"),
  LRRSDS_OEE %>%
    select(c(FactorN, MaxD, LRRSD)) %>%
    mutate(Matching="OEE", Reference="GSA"),
  LRRSDS_OMNI %>%
    select(c(FactorN, MaxD, LRRSD)) %>%
    mutate(Matching="OMNI", Reference="GSA"))

LRRSDS_ALL %>%
  group_by(Matching, Reference, FactorN, MaxD) %>%
  summarise(Min = min(LRRSD),
            Med = median(LRRSD),
            Mean = mean(LRRSD),
            Max = max(LRRSD),
            StDev = sd(LRRSD),
            IQR = IQR(LRRSD),
            .groups="keep") %>%
  write_tsv("TABLES/Table6.tsv", col_names=TRUE)
