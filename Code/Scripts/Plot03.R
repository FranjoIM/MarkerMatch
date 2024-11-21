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

# Plot Omni Data
NUMBERSNPS <- read_delim("DATA/UsedSNPs_OMNI.txt", delim="\t", col_names=T) %>%
  select(-c(Name)) %>%
  summarize_all(sum, na.rm=TRUE) %>%
  rownames_to_column() %>%
  pivot_longer(!rowname, names_to="col1", values_to="col2") %>% 
  select(-rowname) %>%
  separate(col1, c("Factor", "D_MAX")) %>%
  rename(N_SNP=col2) %>%
  mutate(OMNI_Coverage=round(N_SNP/nrow(OMNI_MAN), digits=2),
         GSA_Coverage=round(N_SNP/nrow(GSA_MAN), digits=2),
         D_MAX_LOG=log10(as.numeric(D_MAX)))

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
  ggplot(aes(x=D_MAX_LOG, y=OMNI_Coverage, color=FactorF)) +
  geom_hline(aes(yintercept=H2, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=H1, color="Full Set"), size=1) +
  geom_line(size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="COVERAGE",
       subtitle="OMNI Coverage, \nOMNI matched on GSA",
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
  ggplot(aes(x=D_MAX_LOG, y=GSA_Coverage, color=FactorF)) +
  geom_hline(aes(yintercept=H4, color="Perfect Match"), size=1) +
  geom_hline(aes(yintercept=1, color="Full Set"), size=1) +
  geom_line(size=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="COVERAGE",
       subtitle="GSA Coverage, \nOMNI matched on GSA",
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

NUMBERSNPS_OMNI <- NUMBERSNPS

# Plot GSA Data
NUMBERSNPS <- read_delim("DATA/UsedSNPs_GSA.txt", delim="\t", col_names=T) %>%
  select(-c(Name)) %>%
  summarize_all(sum, na.rm=TRUE) %>%
  rownames_to_column() %>%
  pivot_longer(!rowname, names_to="col1", values_to="col2") %>% 
  select(-rowname) %>%
  separate(col1, c("Factor", "D_MAX")) %>%
  rename(N_SNP=col2) %>%
  mutate(OEE_Coverage=round(N_SNP/nrow(OEE_MAN), digits=2),
         GSA_Coverage=round(N_SNP/nrow(GSA_MAN), digits=2),
         D_MAX_LOG=log10(as.numeric(D_MAX)))

I1 <- NUMBERSNPS %>%
  filter(Factor=="FullSet") %>%
  pull(GSA_Coverage)

I2 <- NUMBERSNPS %>%
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

PLOT_V2_GSA_COVERAGE <- NUMBERSNPS %>%
  filter(!FactorF %in% c("Full Set", "Perfect Match")) %>%
  ggplot(aes(x=D_MAX_LOG, y=GSA_Coverage, color=FactorF)) +
  geom_hline(aes(yintercept=I2, color="Perfect Match"), linewidth=1) +
  geom_hline(aes(yintercept=I1, color="Full Set"), linewidth=1) +
  geom_point(position = position_dodge(0.01), size=3, shape=17) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="COVERAGE",
       subtitle="GSA Coverage, \nOEE matched on GSA",
       color="MATCHING METHOD") +
  scale_x_continuous(breaks=4, n.breaks=1, labels=4, limits=c(3.98, 4.02)) +
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

NUMBERSNPS_GSA <- NUMBERSNPS

# Plot OEE Data
NUMBERSNPS <- read_delim("DATA/UsedSNPs_OEE.txt", delim="\t", col_names=T) %>%
  select(-c(Name)) %>%
  summarize_all(sum, na.rm=TRUE) %>%
  rownames_to_column() %>%
  pivot_longer(!rowname, names_to="col1", values_to="col2") %>% 
  select(-rowname) %>%
  separate(col1, c("Factor", "D_MAX")) %>%
  rename(N_SNP=col2) %>%
  mutate(OEE_Coverage=round(N_SNP/nrow(OEE_MAN), digits=2),
         GSA_Coverage=round(N_SNP/nrow(GSA_MAN), digits=2),
         D_MAX_LOG=log10(as.numeric(D_MAX)))

J1 <- NUMBERSNPS %>%
  filter(Factor=="FullSet") %>%
  pull(OEE_Coverage)

J2 <- NUMBERSNPS %>%
  filter(Factor=="PerfectMatch") %>%
  pull(OEE_Coverage)

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

PLOT_V3_OEE_COVERAGE <- NUMBERSNPS %>%
  filter(!FactorF %in% c("Full Set", "Perfect Match")) %>%
  ggplot(aes(x=D_MAX_LOG, y=OEE_Coverage, color=FactorF)) +
  geom_hline(aes(yintercept=J2, color="Perfect Match"), linewidth=1) +
  geom_hline(aes(yintercept=J1, color="Full Set"), linewidth=1) +
  geom_point(position = position_dodge(0.01), size=3, shape=17) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4", "red3", "steelblue3"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match", "Full Set")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="COVERAGE",
       subtitle="OEE Coverage, \nOEE matched on GSA",
       color="MATCHING METHOD") +
  scale_x_continuous(breaks=4, n.breaks=1, labels=4, limits=c(3.98, 4.02)) +
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

NUMBERSNPS_OEE <- NUMBERSNPS

# Save Plots
ggsave(plot=PLOT_V1_OMNI_COVERAGE,
       filename="FIGURES/Plot3_A.png",
       device="png",
       width=10,
       height=7,
       units="in",
       dpi=350,
       bg="white")

ggsave(plot=PLOT_V1_GSA_COVERAGE,
       filename="FIGURES/Plot3_B.png",
       device="png",
       width=10,
       height=7,
       units="in",
       dpi=350,
       bg="white")

ggsave(plot=PLOT_V3_OEE_COVERAGE,
       filename="FIGURES/Plot3_C.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

ggsave(plot=PLOT_V2_GSA_COVERAGE,
       filename="FIGURES/Plot3_D.png",
       device="png",
       width=6,
       height=7,
       units="in",
       dpi=350,
       bg="white")

# Arrange and save plots
Plot3 <- ggarrange(
          ggarrange(PLOT_V1_OMNI_COVERAGE + theme(legend.position="none", plot.subtitle=element_blank()), 
            PLOT_V1_GSA_COVERAGE + theme(legend.position="none", plot.subtitle=element_blank()),
            align="hv", labels=c("A", "B"), nrow=2, ncol=1, legend="top", common.legend=T),        
          ggarrange(PLOT_V3_OEE_COVERAGE + theme(legend.position="none", plot.subtitle=element_blank()),
            PLOT_V2_GSA_COVERAGE + theme(legend.position="none", plot.subtitle=element_blank()),
            align="hv", labels=c("C", "D"), nrow=1, ncol=2),
          align="hv", nrow=2, ncol=1, legend="top", common.legend=T, heights=c(1, 0.5))

ggsave(plot=Plot3,
       filename="FIGURES/Plot3.png",
       device="png",
       width=7,
       height=10,
       units="in",
       dpi=350,
       bg="white")

# OUTPUT TABLE WITH COVERAGE DATA
bind_rows(mutate(NUMBERSNPS_OMNI, Figure_Reference="Figure 3A-B"),
          mutate(NUMBERSNPS_OEE, Figure_Reference="Figure 3C"),
          mutate(NUMBERSNPS_GSA, Figure_Reference="Figure 3D")) %>%
  select(c(Figure_Reference, FactorN, D_MAX, N_SNP, OMNI_Coverage, GSA_Coverage, OEE_Coverage)) %>%
  rename(Factor=FactorN) %>%
  write_tsv("TABLES/Table4.tsv", col_names=TRUE)
