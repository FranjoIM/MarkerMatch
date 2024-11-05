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

### GAPS PLOTS FOR OMNI-GSA
GAPS_OMNI <- data.frame()
SETS_OMNI <- colnames(OMNI_SNPS[2:51])

for(j in SETS_OMNI){
  FNAME <- OMNI_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
  for(k in seq(1:22)){
    GAP <- OMNI_MAN %>%
      filter(Chr==k) %>% 
      filter(Name %in% FNAME) %>%
      select(Position) %>%
      arrange(Position) %>%
      mutate(Lag=lag(Position)) %>%
      mutate(GAP=Position-Lag) %>%
      filter(!is.na(GAP)) %>%
      mutate(LOG10=ifelse(GAP==0, 0, log10(GAP))) %>%
      pull(LOG10)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      Chr=k,
                      Gaps=GAP)
    
    GAPS_OMNI <- rbind(GAPS_OMNI, ROW)
  }
}

GAPS_OMNI$MaxDF <- factor(GAPS_OMNI$MaxD, 
                     levels=c("0", "10", "50", "100", "500", "1000", "5000", "10000", 
                              "50000", "100000", "500000", "1000000", "5000000"),
                     labels=c("0", "10", "50", "100", "500", "1,000", "5,000", "10,000", 
                              "50,000", "100,000", "500,000", "1,000,000", "5,000,000"))

GAPS_OMNI <- GAPS_OMNI %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))

GAPS_OMNI$FactorF <- factor(GAPS_OMNI$FactorN,
                       levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                "LRR sd", "Distance"),
                       labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                "LRR sd", "Distance"))

H1_OMNI <- GAPS_OMNI %>%
  filter(Factor=="PerfectMatch") %>%
  pull(Gaps) %>%
  median(., na.rm=TRUE)

H2_OMNI <- GAPS_OMNI %>%
  filter(Factor=="FullSet") %>%
  pull(Gaps) %>%
  median(., na.rm=TRUE)

Plot4_A <- GAPS_OMNI %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=Gaps)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_OMNI, color="Perfect Match"), linewidth=1) +
  geom_hline(aes(yintercept=H2_OMNI, color="Full Set"), linewidth=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ GAP ]")),
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

ggsave(plot=Plot4_A,
       filename="FIGURES/Plot4_A.png",
       device="png",
       width=11,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### GAPS PLOTS FOR OEE-GSA
GAPS_OEE <- data.frame()
SETS_OEE <- colnames(OEE_SNPS[2:7])

for(j in SETS_OEE){
  FNAME <- OEE_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
  for(k in seq(1:22)){
    GAP <- OEE_MAN %>%
      filter(Chr==k) %>% 
      filter(Name %in% FNAME) %>%
      select(Position) %>%
      arrange(Position) %>%
      mutate(Lag=lag(Position)) %>%
      mutate(GAP=Position-Lag) %>%
      filter(!is.na(GAP)) %>%
      mutate(LOG10=ifelse(GAP==0, 0, log10(GAP))) %>%
      pull(LOG10)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      Chr=k,
                      Gaps=GAP)
    
    GAPS_OEE <- rbind(GAPS_OEE, ROW)
  }
}

GAPS_OEE$MaxDF <- factor(GAPS_OEE$MaxD, 
                          levels=c("0", "10000"),
                          labels=c("0", "10,000"))

GAPS_OEE <- GAPS_OEE %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))

GAPS_OEE$FactorF <- factor(GAPS_OEE$FactorN,
                            levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                     "LRR sd", "Distance"),
                            labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                     "LRR sd", "Distance"))

H1_OEE <- GAPS_OEE %>%
  filter(Factor=="PerfectMatch") %>%
  pull(Gaps) %>%
  median(., na.rm=TRUE)

H2_OEE <- GAPS_OEE %>%
  filter(Factor=="FullSet") %>%
  pull(Gaps) %>%
  median(., na.rm=TRUE)

Plot4_B <- GAPS_OEE %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=Gaps)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_OEE, color="Perfect Match"), linewidth=1) +
  geom_hline(aes(yintercept=H2_OEE, color="Full Set"), linewidth=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ GAP ]")),
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

print(Plot4_B)

ggsave(plot=Plot4_B,
       filename="FIGURES/Plot4_B.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

### GAPS PLOTS FOR GSA-OEE
GAPS_GSA <- data.frame()
SETS_GSA <- colnames(GSA_SNPS[2:7])

for(j in SETS_GSA){
  FNAME <- GSA_SNPS %>%
    select(c("Name", all_of(j))) %>%
    filter(!is.na(.[[2]])) %>%
    pull(Name)
  
  for(k in seq(1:22)){
    GAP <- GSA_MAN %>%
      filter(Chr==k) %>% 
      filter(Name %in% FNAME) %>%
      select(Position) %>%
      arrange(Position) %>%
      mutate(Lag=lag(Position)) %>%
      mutate(GAP=Position-Lag) %>%
      filter(!is.na(GAP)) %>%
      mutate(LOG10=ifelse(GAP==0, 0, log10(GAP))) %>%
      pull(LOG10)
    
    a <- str_split(j, pattern="_")[[1]][1]
    b <- str_split(j, pattern="_")[[1]][2]
    
    ROW <- data.frame(Factor=a,
                      MaxD=b,
                      Chr=k,
                      Gaps=GAP)
    
    GAPS_GSA <- rbind(GAPS_GSA, ROW)
  }
}

GAPS_GSA$MaxDF <- factor(GAPS_GSA$MaxD, 
                          levels=c("0", "10000"),
                          labels=c("0", "10,000"))

GAPS_GSA <- GAPS_GSA %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))

GAPS_GSA$FactorF <- factor(GAPS_GSA$FactorN,
                            levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                     "LRR sd", "Distance"),
                            labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
                                     "LRR sd", "Distance"))

H1_GSA <- GAPS_GSA %>%
  filter(Factor=="PerfectMatch") %>%
  pull(Gaps) %>%
  median(., na.rm=TRUE)

H2_GSA <- GAPS_GSA %>%
  filter(Factor=="FullSet") %>%
  pull(Gaps) %>%
  median(., na.rm=TRUE)

Plot4_C <- GAPS_GSA %>%
  ggplot(aes(x=FactorF, color=FactorF, fill=FactorF, y=Gaps)) +
  geom_boxplot(position=position_dodge2(preserve="single"), alpha=0.2) +
  geom_hline(aes(yintercept=H1_GSA, color="Perfect Match"), linewidth=1) +
  geom_hline(aes(yintercept=H2_GSA, color="Full Set"), linewidth=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  scale_fill_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                    breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y=expression(bold("LOG"["10"] ~ "[ GAP ]")),
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

print(Plot4_C)

ggsave(plot=Plot4_C,
       filename="FIGURES/Plot4_C.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# SAVE ALL GAPS FILES AS RDATA
save(GAPS_OMNI, GAPS_OEE, GAPS_GSA, file="Table5.RData")

# SAVE THE PANNELED PLOT
ggarrange(ggarrange(Plot4_A + theme(legend.position="none", plot.subtitle=element_blank()),
            align="hv", labels=c("A"), nrow=1, ncol=1, legend="top", common.legend=T),
          ggarrange(Plot4_B + theme(legend.position="none", plot.subtitle=element_blank()),
            Plot4_C + theme(legend.position="none", plot.subtitle=element_blank()),
            align="hv", labels=c("B", "C"), nrow=1, ncol=2),
          align="hv", nrow=2, ncol=1, legend="top", common.legend=T, heights=c(1, 1)) %>%
ggsave(filename="FIGURES/Plot4.png",
       device="png",
       width=8,
       height=10,
       units="in",
       dpi=350,
       bg="white")

# TABULATE GAP SIZES ACROSS TABLES AND CONDITIONS
GAPS_ALL <- bind_rows(
  GAPS_OMNI %>%
    select(c(FactorN, MaxD, Gaps)) %>%
    mutate(Figure_Reference="Figure S2A"),
  GAPS_OEE %>%
    select(c(FactorN, MaxD, Gaps)) %>%
    mutate(Figure_Reference="Figure S2B"),
  GAPS_GSA %>%
    select(c(FactorN, MaxD, Gaps)) %>%
    mutate(Figure_Reference="Figure S2C"))

GAPS_ALL %>%
  group_by(Figure_Reference, FactorN, MaxD) %>%
  summarise(Min = min(Gaps),
            Med = median(Gaps),
            Mean = mean(Gaps),
            Max = max(Gaps),
            StDev = sd(Gaps),
            IQR = IQR(Gaps),
            .groups="keep") %>%
 write_tsv("TABLES/Table5.tsv", col_names=TRUE)
