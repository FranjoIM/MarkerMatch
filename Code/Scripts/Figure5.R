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

### BAFS PLOTS FOR OMNI-GSA
BAFS_OMNI <- data.frame()
SETS_OMNI <- colnames(OMNI_SNPS[2:51])

for(j in SETS_OMNI){
  FNAME <- OMNI_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
  BAF <- OMNI_MAN %>%
    filter(Name %in% FNAME) %>%
    pull(BAF)
  
  a <- str_split(j, pattern="_")[[1]][1]
  b <- str_split(j, pattern="_")[[1]][2]
  
  ROW <- data.frame(Factor=a,
                    MaxD=b,
                    BAF=BAF)
  
  BAFS_OMNI <- rbind(BAFS_OMNI, ROW)
}

BAFS_OMNI$MaxDF <- factor(BAFS_OMNI$MaxD, 
                          levels=c("0", "10", "50", "100", "500", "1000", "5000", "10000", 
                                   "50000", "100000", "500000", "1000000", "5000000"),
                          labels=c("0", "10", "50", "100", "500", "1,000", "5,000", "10,000", 
                                   "50,000", "100,000", "500,000", "1,000,000", "5,000,000"))

BAFS_OMNI <- BAFS_OMNI %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))  %>%
  drop_na(BAF)

BAFS_OMNI$FactorF <- factor(BAFS_OMNI$FactorN,
                            levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                     "LRR sd", "Distance"),
                            labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                     "LRR sd", "Distance"))

H1_OMNI <- BAFS_OMNI %>%
  filter(Factor=="PerfectMatch") %>%
  pull(BAF) %>%
  median(., na.rm=TRUE)

H2_OMNI <- BAFS_OMNI %>%
  filter(Factor=="FullSet") %>%
  pull(BAF) %>%
  median(., na.rm=TRUE)

Plot5_A <- BAFS_OMNI %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=BAF)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_OMNI, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H2_OMNI, color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="BAF",
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

ggsave(plot=Plot5_A,
       filename="FIGURES/Plot5_A.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### BAFS PLOTS FOR OEE-GSA
BAFS_OEE <- data.frame()
SETS_OEE <- colnames(OEE_SNPS[2:7])

for(j in SETS_OEE){
  FNAME <- OEE_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
    BAF <- OEE_MAN %>%
      filter(Name %in% FNAME) %>%
      pull(BAF)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      BAF=BAF)
    
    BAFS_OEE <- rbind(BAFS_OEE, ROW)
}

BAFS_OEE$MaxDF <- factor(BAFS_OEE$MaxD, 
                         levels=c("0", "10000"),
                         labels=c("0", "10,000"))

BAFS_OEE <- BAFS_OEE %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))  %>%
  drop_na(BAF)

BAFS_OEE$FactorF <- factor(BAFS_OEE$FactorN,
                           levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                    "LRR sd", "Distance"),
                           labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                    "LRR sd", "Distance"))

H1_OEE <- BAFS_OEE %>%
  filter(Factor=="PerfectMatch") %>%
  pull(BAF) %>%
  median(., na.rm=TRUE)

H2_OEE <- BAFS_OEE %>%
  filter(Factor=="FullSet") %>%
  pull(BAF) %>%
  median(., na.rm=TRUE)

Plot5_B <- BAFS_OEE %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=BAF)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_OEE, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H2_OEE, color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="BAF",
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

ggsave(plot=Plot5_B,
       filename="FIGURES/Plot5_B.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### BAFS PLOTS FOR GSA-OEE
BAFS_GSA <- data.frame()
SETS_GSA <- colnames(GSA_SNPS[2:7])

for(j in SETS_GSA){
  FNAME <- GSA_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
    BAF <- GSA_MAN %>%
      filter(Name %in% FNAME) %>%
      pull(BAF)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      BAF=BAF)
    
    BAFS_GSA <- rbind(BAFS_GSA, ROW)
}

BAFS_GSA$MaxDF <- factor(BAFS_GSA$MaxD, 
                         levels=c("0", "10000"),
                         labels=c("0", "10,000"))

BAFS_GSA <- BAFS_GSA %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))  %>%
  drop_na(BAF)

BAFS_GSA$FactorF <- factor(BAFS_GSA$FactorN,
                           levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                    "LRR sd", "Distance"),
                           labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                    "LRR sd", "Distance"))

H1_GSA <- BAFS_GSA %>%
  filter(Factor=="PerfectMatch") %>%
  pull(BAF) %>%
  median(., na.rm=TRUE)

H2_GSA <- BAFS_GSA %>%
  filter(Factor=="FullSet") %>%
  pull(BAF) %>%
  median(., na.rm=TRUE)

Plot5_C <- BAFS_GSA %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=BAF)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_GSA, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H2_GSA, color="Full Set"), size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="BAF",
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

ggsave(plot=Plot5_C,
       filename="FIGURES/Plot5_C.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# SAVE ALL BAF FILES AS RDATA
save(BAFS_OMNI, BAFS_OEE, BAFS_GSA, file="Table6.RData")

# TABULATE BAFS ACROSS TABLES AND CONDITIONS
BAFS_ALL <- bind_rows(
  BAFS_GSA %>%
    select(c(FactorN, MaxD, BAF)) %>%
    mutate(Matching="GSA", Reference="OEE"),
  BAFS_OEE %>%
    select(c(FactorN, MaxD, BAF)) %>%
    mutate(Matching="OEE", Reference="GSA"),
  BAFS_OMNI %>%
    select(c(FactorN, MaxD, BAF)) %>%
    mutate(Matching="OMNI", Reference="GSA"))

BAFS_ALL %>%
  group_by(Matching, Reference, FactorN, MaxD) %>%
  summarise(Min = min(BAF),
            Med = median(BAF),
            Mean = mean(BAF),
            Max = max(BAF),
            StDev = sd(BAF),
            IQR = IQR(BAF),
            .groups="keep") %>%
  write_tsv("TABLES/Table6.tsv", col_names=TRUE)
