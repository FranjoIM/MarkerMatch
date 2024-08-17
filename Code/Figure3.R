# Set working directory
W_DIR <- "..."
setwd(W_DIR)

# Load packages
library(tidyverse)
library(ggplot)

# Import array manifests
OMNI_MAN <- read_delim("DATA/OMNI_MAN.txt", delim="\t")
GSA_MAN <- read_delim("DATA/GSA_MAN.txt", delim="\t")
OEE_MAN <- read_delim("DATA/OEE_MAN.txt", delim="\t")

### Graph first-step validation data
# Import the MAT_SELECT information on kept SNPs
FACTORS <- c("Pos", "BAF", "LRRsd", "LRRmean")
MAXDS <- c("10", "50", "100", "500", "1000", "5000", "10000", "50000", 
           "100000", "500000", "1000000", "5000000")

PerfectMatch <- read_delim("DATA/PerfectMatch_MatSelect.csv", delim="\t", col_names=F) %>%
  mutate(PerfectMatch_0=1)

USEDSNPS <- OMNI_MAN %>%
  select(Name) %>%
  mutate(FullSet_0=1) %>%
  left_join(PerfectMatch, by=c("Name"="X1"))

for(i in FACTORS){
  for(j in MAXDS){
    Name <- quo_name(paste0(i, "_", j))
    
    DF_MAT <- read_delim(paste0("DATA/", i, "_", j, "_MatSelect.csv"),
                         delim="\t", col_names=F) %>%
      distinct(.) %>%
      mutate(X2=1) %>%
      rename(!!Name := X2)
    
    USEDSNPS <- USEDSNPS %>%
      left_join(DF_MAT, by=c("Name"="X1"))
  }
}

### NUMBER SNPS PLOTS

NUMBERSNPS <- USEDSNPS %>%
  select(-c(Name)) %>%
  summarize_all(sum, na.rm=TRUE) %>%
  rownames_to_column() %>%
  pivot_longer(!rowname, names_to="col1", values_to="col2") %>% 
  select(-rowname) %>%
  separate(col1, c("Factor", "MaxD")) %>%
  rename(N_SNP=col2) %>%
  mutate(OMNI_Coverage=round(N_SNP/2376441, digits=2),
         GSA_Coverage=round(N_SNP/672590, digits=2),
         MaxD_LOG=log10(as.numeric(MaxD)))

H1 <- NUMBERSNPS %>%
  filter(Factor=="FullSet") %>%
  pull(OMNI_Coverage)

H2 <- NUMBERSNPS %>%
  filter(Factor=="PerfectMatch") %>%
  pull(OMNI_Coverage)

H3 <- NUMBERSNPS %>%
  filter(Factor=="FullSet") %>%
  pull(GSA_Coverage)

H4 <- NUMBERSNPS %>%
  filter(Factor=="PerfectMatch") %>%
  pull(GSA_Coverage)

NUMBERSNPS <- NUMBERSNPS %>%
  mutate(FactorN=case_when(
    Factor=="PerfectMatch" ~ "Perfect Match",
    Factor=="FullSet" ~ "Full Set",
    Factor=="BAF" ~ "BAF",
    Factor=="LRRmean" ~ "LRR mean",
    Factor=="LRRsd" ~ "LRR sd",
    Factor=="Pos" ~ "Distance",
    TRUE ~ NA_character_))


NUMBERSNPS$FactorF <- factor(NUMBERSNPS$FactorN,
  levels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
          "LRR sd", "Distance"),
  labels=c("Full Set", "Perfect Match", "BAF", "LRR mean", 
          "LRR sd", "Distance"))

PLOT_V1_OMNI_COVERAGE <- NUMBERSNPS %>%
  filter(!FactorF %in% c("Full Set", "Perfect Match")) %>%
  ggplot(aes(x=MaxD_LOG, y=OMNI_Coverage, color=FactorF)) +
  geom_hline(aes(yintercept=H2, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H1, color="Full Set"), size=1) +
  geom_line(size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[MAXIMUM MATCHING DISTANCE]")),
       y="COVERAGE",
       subtitle="OMNI Coverage",
       color="MATCHING METHOD") +
  ylim(0, 1) +
  theme_bw() + 
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(face="bold.italic", size=12)) +
  guides(color=guide_legend(title.position="top", nrow=1))

PLOT_V1_GSA_COVERAGE <- NUMBERSNPS %>%
  filter(!FactorF %in% c("Full Set", "Perfect Match")) %>%
  ggplot(aes(x=MaxD_LOG, y=GSA_Coverage, color=FactorF)) +
  geom_hline(aes(yintercept=H4, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=1, color="Full Set"), size=1) +
  geom_line(size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[MAXIMUM MATCHING DISTANCE]")),
       y="COVERAGE",
       subtitle="GSA Coverage",
       color="MATCHING METHOD") +
  ylim(0, 1) +
  theme_bw() + 
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=12, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(face="bold.italic", size=12)) +
  guides(color=guide_legend(title.position="top", nrow=1))
