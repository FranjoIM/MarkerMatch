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
                    D_MAX=b,
                    BAF=BAF)
  
  BAFS_OMNI <- rbind(BAFS_OMNI, ROW)
}

BAFS_OMNI$D_MAXF <- factor(BAFS_OMNI$D_MAX, 
                           levels=c("0", "10", "50", "100", "500", "1000", "5000", "10000", 
                                    "50000", "100000", "500000", "1000000", "5000000"),
                           labels=c("0", "10", "50", "100", "500", "1,000", "5,000", "10,000", 
                                    "50,000", "100,000", "500,000", "1,000,000", "5,000,000"))

BAFS_OMNI <- BAFS_OMNI %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Exact Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))  %>%
  drop_na(BAF)

BAFS_OMNI$FactorF <- factor(BAFS_OMNI$FactorN,
                            levels=c("Full Set", "Exact Match", "BAF", "LRR mean", 
                                     "LRR sd", "Distance"),
                            labels=c("Full Set", "Exact Match", "BAF", "LRR mean", 
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
  geom_hline(aes(yintercept=H1_OMNI, color="Exact Match"), linewidth=1) +
  geom_hline(aes(yintercept=H2_OMNI, color="Full Set"), linewidth=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Exact Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Exact Match", "Full Set")) +
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
  facet_grid(. ~ D_MAXF, scales="free_x", space="free_x", switch="both")

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
                    D_MAX=b,
                    BAF=BAF)
  
  BAFS_OEE <- rbind(BAFS_OEE, ROW)
}

BAFS_OEE$D_MAXF <- factor(BAFS_OEE$D_MAX, 
                          levels=c("0", "10000"),
                          labels=c("0", "10,000"))

BAFS_OEE <- BAFS_OEE %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Exact Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))  %>%
  drop_na(BAF)

BAFS_OEE$FactorF <- factor(BAFS_OEE$FactorN,
                           levels=c("Full Set", "Exact Match", "BAF", "LRR mean", 
                                    "LRR sd", "Distance"),
                           labels=c("Full Set", "Exact Match", "BAF", "LRR mean", 
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
  geom_hline(aes(yintercept=H1_OEE, color="Exact Match"), linewidth=1) +
  geom_hline(aes(yintercept=H2_OEE, color="Full Set"), linewidth=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Exact Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Exact Match", "Full Set")) +
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
  facet_grid(. ~ D_MAXF, scales="free_x", space="free_x", switch="both")

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
                    D_MAX=b,
                    BAF=BAF)
  
  BAFS_GSA <- rbind(BAFS_GSA, ROW)
}

BAFS_GSA$D_MAXF <- factor(BAFS_GSA$D_MAX, 
                          levels=c("0", "10000"),
                          labels=c("0", "10,000"))

BAFS_GSA <- BAFS_GSA %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Exact Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))  %>%
  drop_na(BAF)

BAFS_GSA$FactorF <- factor(BAFS_GSA$FactorN,
                           levels=c("Full Set", "Exact Match", "BAF", "LRR mean", 
                                    "LRR sd", "Distance"),
                           labels=c("Full Set", "Exact Match", "BAF", "LRR mean", 
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
  geom_hline(aes(yintercept=H1_GSA, color="Exact Match"), linewidth=1) +
  geom_hline(aes(yintercept=H2_GSA, color="Full Set"), linewidth=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Exact Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Exact Match", "Full Set")) +
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
  facet_grid(. ~ D_MAXF, scales="free_x", space="free_x", switch="both")

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

# SAVE THE PANNELED PLOT
ggarrange(ggarrange(Plot5_A + theme(legend.position="none", plot.subtitle=element_blank()),
                    align="hv", labels=c("A"), nrow=1, ncol=1, legend="top", common.legend=T),
          ggarrange(Plot5_B + theme(legend.position="none", plot.subtitle=element_blank()),
                    Plot5_C + theme(legend.position="none", plot.subtitle=element_blank()),
                    align="hv", labels=c("B", "C"), nrow=1, ncol=2),
          align="hv", nrow=2, ncol=1, legend="top", common.legend=T, heights=c(1, 1)) %>%
  ggsave(filename="FIGURES/Plot5.png",
         device="png",
         width=10,
         height=10,
         units="in",
         dpi=350,
         bg="white")

# TABULATE BAFS ACROSS TABLES AND CONDITIONS
BAFS_ALL <- bind_rows(
  BAFS_OMNI %>%
    select(c(FactorN, D_MAX, BAF)) %>%
    mutate(Figure_Reference="Figure S3A"),
  BAFS_OEE %>%
    select(c(FactorN, D_MAX, BAF)) %>%
    mutate(Figure_Reference="Figure S3B"),
  BAFS_GSA %>%
    select(c(FactorN, D_MAX, BAF)) %>%
    mutate(Figure_Reference="Figure S3C"))

BAFS_ALL %>%
  group_by(Figure_Reference, FactorN, D_MAX) %>%
  summarise(Min = min(BAF),
            Med = median(BAF),
            Mean = mean(BAF),
            Max = max(BAF),
            StDev = sd(BAF),
            IQR = IQR(BAF),
            .groups="keep") %>%
  write_tsv("TABLES/TableS1D.tsv", col_names=TRUE)
