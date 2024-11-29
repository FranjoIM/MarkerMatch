# SET WORKING DIRECTORY, PREPARE THE ENVIRONMENT
setwd("...")
library(tidyverse)

# IMPORT DATA FOR PLOTTING
DF <- read_delim("TABLES/TableS1A.tsv") %>%
  mutate(D_MAX_LOG=log10(D_MAX))

DF$FactorF <- factor(DF$Factor,
                     levels=c("BAF", "LRR mean", "LRR sd", "Distance"),
                     labels=c("BAF", "LRR mean", "LRR sd", "Distance"))

(DF %>%
  rename(R=R_Execution, 
         Python=Python_Execution) %>%
  pivot_longer(c(R, Python), names_to="Language", values_to="Times") %>%
  ggplot(aes(x=D_MAX_LOG, y=Times, linetype=Language, color=FactorF)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("goldenrod1", "slateblue2", "seagreen4", "lightsalmon4"),
                     breaks=c("BAF", "LRR mean", "LRR sd", "Distance", "Perfect Match")) +
  labs(x=expression(bold("LOG"["10"] ~ "[" ~"D"["MAX"] ~ "]")),
       y="EXECUTION TIME",
       linetype="LANGUAGE",
       color="FACTOR") +
  theme_bw() + 
  theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(angle=90, vjust=0.5, hjust=0.5, size=12),
        axis.title.x=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.title.y=element_text(size=15, hjust=0, vjust=0, face="bold", margin=margin(t=0, r=10, b=0, l=0, unit="pt")),
        plot.subtitle=element_text(size=15, hjust=0, vjust=0, face="bold.italic"),
        legend.position="top",
        legend.justification="left",
        legend.title=element_text(size=12, face="bold")) +
  guides(color=guide_legend(title.position="top", nrow=1),
         linetype=guide_legend(title.position="top")) )%>%
  ggsave(filename="FIGURES/Plot2.png",
         device="png",
         width=9,
         height=7,
         units="in",
         dpi=350,
         bg="white")
