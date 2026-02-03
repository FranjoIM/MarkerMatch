##' Summarize_LRR
##'
##' Match genetic markers from a reference and matching manifest based on SNP name,
##' genomic position and/or allele frequency / LRR differences. The function supports
##' several matching methods and writes selected name lists and a summary CSV.
##'
##' @param Input Input file path containing LRR values for each sample (columns) and position (rows) from GenomeStudio.
##' @param Output Output file path to which LRR summaries will be written.
##' @param N Number of samples in the file.
##' @param Max_Rows Number of rows to process at a time. Default 20,000.
##' @return A data frame containing LRR mean and sd values for each probe.
##' @details This function calculates LRR statistics for each probe across all samples.
##' @importFrom rlang .data
##' @importFrom dplyr  mutate pick ends_with rename
##' @importFrom magrittr %>%
##' @importFrom readr write_delim read_delim
##' @importFrom matrixStats rowSds
##' @export

Summarize_LRR <- function(Input=Input, Output=Output, N=N, Max_Rows=20000){

  # COUNT NUMBER OF LINES IN THE LRR FILE
  f <- file(Input, open="rb")

  message("Processing LRR file ", Input)

  nlines <- 0L
  while(length(chunk <- readBin(con=f, what="raw", n=65536L))>0){
    nlines <- nlines + sum(chunk==as.raw(10L))
  }
  close(f); rm(f, chunk)

  # INITIALIZE HOLDING DATA FRAME AND COUNTER VARIABLES
  DF <- data.frame(NULL)
  s <- 1L

  # LOOP THROUGH THE LRR FILES AND CALCULATE THE METRICS
  while(s<nlines){

    # IMPORT FRAGMENT OF THE FILE
    DAT <- read_delim(file=Input, col_names=FALSE, col_types=c(rep("c", 2), rep("d", (N+1))), na=c("", "NA", "NaN"), skip=s, n_max=Max_Rows, lazy = FALSE)

    MOD <- DAT %>%
      mutate(LRR_mean=rowMeans(pick(ends_with(".Log R Ratio")), na.rm=TRUE),
             LRR_sd=rowSds(pick(ends_with(".Log R Ratio")), na.rm=TRUE),
             .keep="unused")

    DF <- rbind(DF, MOD)

    message("Processed rows ", s, " to ", s+Max_Rows, ". ", round(s/nlines*100, digits=2), "% complete.")

    s <- s + Max_Rows
  }

  # FORMAT AND EXPORT FILE
  DF %>% rename(Name=1, Chr=2, Position=3) %>% write_delim(Output, delim="\t", col_names=TRUE, na="NA")

  message("Process complete. LRR summary file written to ", Output)
}

# -----------------------------------------------------------------------------
# End of Summarize_LRR.R
# -----------------------------------------------------------------------------
