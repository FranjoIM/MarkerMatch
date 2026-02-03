##' Marker_Match
##'
##' Match genetic markers from a reference and matching manifest based on SNP name,
##' genomic position and/or allele frequency / LRR differences. The function supports
##' several matching methods and writes selected name lists and a summary CSV.
##'
##' @param Reference A data frame with columns: \code{Name}, \code{Chr}, \code{Position}, \code{BAF}, \code{LRR_mean}, \code{LRR_sd}.
##' @param Matching A data frame with columns: \code{Name}, \code{Chr}, \code{Position}, \code{BAF}, \code{LRR_mean}, \code{LRR_sd}.
##' @param Method A string specifying the matching method. One of "Distance", "BAF_delta", "LRR_mean_delta", or "LRR_sd_delta".
##' @param D_MAX Numeric maximum distance (bp) for matching markers by position.
##' @param OutPath Character path (without extension) used to write output CSV files.
##' @param Chunk Integer chunk/window size used when sliding across genome regions. Default 10,000,000.
##' @param Verbose Logical indicating whether to print detailed progress messages. Default FALSE.
##' @return A data frame containing matched marker pairs and computed distance/delta columns.
##' @details The function attempts exact name matches first, then identical-position matches,
##' and finally nearest matches within `D_MAX` using a sliding window defined by `Chunk`.
##' @importFrom rlang .data
##' @importFrom dplyr filter mutate select group_by ungroup row_number full_join arrange_at cross_join
##' @importFrom magrittr %>%
##' @importFrom stringr str_sort
##' @importFrom readr write_csv
##' @importFrom tidyr drop_na
##' @export

Marker_Match <- function(Reference=Reference, Matching=Matching, Method=Method, D_MAX=D_MAX, OutPath=OutPath, Chunk=10000000, Verbose=FALSE){

  # CHECK FOR COMMON ERRORS IN FUNCTION INPUTS
  TEST_COLNAME <- c("Name","Chr","Position","BAF","LRR_mean","LRR_sd")
  if(sum(intersect(names(Reference), TEST_COLNAME) %in% TEST_COLNAME)!=6){
    stop("Reference is missing one of the following required columns: Name, Chr, Position, BAF, LRR_mean, LRR_sd.")}

  if(sum(intersect(names(Matching), TEST_COLNAME) %in% TEST_COLNAME)!=6){
    stop("Matching is missing one of the following required columns: Name, Chr, Position, BAF, LRR_mean, LRR_sd.")}

  if(!Method %in% c("Distance","BAF_delta","LRR_mean_delta","LRR_sd_delta")){
    stop("Method must be one of the following: Distance, BAF_delta, LRR_mean_delta, or LRR_sd_delta.")}

  if(!is.numeric(D_MAX)){
    stop("D_MAX must be a numeric value.")}

  if(!is.character(OutPath)){
    stop("OutPath must be a string value containing a valid path (without file extensions).")}

  if(!is.numeric(Chunk) | Chunk <= 0){
    stop("Chunk must be a positive numeric value.")}

  if(Chunk < 2*D_MAX){
    stop("Chunk size must be at least twice the D_MAX to ensure proper matching.")}

  if(nrow(Reference)==0){
    stop("Reference manifest is empty.")}

  if(nrow(Matching)==0){
    stop("Matching manifest is empty.")}

  if(nrow(Reference %>% drop_na())==0){
    stop("Reference manifest contains only NA values.")}

  if(nrow(Matching %>% drop_na())==0){
    stop("Matching manifest contains only NA values.")}

  if(Verbose==TRUE){message("Starting MarkerMatch process.\n")}
  start_time <- Sys.time()

  # GET THE LIST OF AVAILABLE CHROMOSOMES
  Chromosomes <- str_sort(unique(Reference$Chr), numeric=TRUE)

  # CREATE HOLDING DATA FRAME
  DF <- data.frame()

  # RUN MARKER MATCH
  for(chr in seq_along(Chromosomes)){
    if(Verbose==TRUE){message("Processing chromosome ", Chromosomes[chr], ".")}

    # FILTER ON CHROMOSOME
    Ref <- Reference %>%
      filter(.data$Chr==Chromosomes[chr])

    Mat <- Matching %>%
      filter(.data$Chr==Chromosomes[chr])

    # PULL PROBES BASED ON IDENTICA NAMES, ANNOTATE DISTANCES, DELTAS, AND NAME IDENTITY
    if(Verbose==TRUE){message("Matching by SNP name -- exact match.")}

    Exact_Match <- full_join(Ref, Mat, by="Name", keep=TRUE, relationship="many-to-many", suffix=c(".Ref", ".Mat")) %>%
      filter(.data$Name.Ref==.data$Name.Mat) %>%
      mutate(Distance=abs(.data$Position.Ref-.data$Position.Mat),
             LRR_sd_delta=abs(.data$LRR_sd.Ref-.data$LRR_sd.Mat),
             LRR_mean_delta=abs(.data$LRR_mean.Ref-.data$LRR_mean.Mat),
             BAF_delta=abs(.data$BAF.Ref-.data$BAF.Mat),
             Position=ifelse(.data$Position.Ref==.data$Position.Mat, "Same", "Different"),
             Name=ifelse(.data$Name.Ref==.data$Name.Mat, "Same", "Different"))

    if(nrow(Exact_Match)>0) {
      # MERGE EXACT MATCH INTO DF
      DF <- rbind(DF, Exact_Match)
    }

    # PURGE MATCHED PROBES
    Mat <- Mat %>%
      filter(!.data$Name %in% DF$Name.Mat)

    Ref <- Ref %>%
      filter(!.data$Name %in% DF$Name.Ref)

    # PULL PROBES WITH SAME POSITION AND CLOSEST METHOD
    if(Verbose==TRUE){message("Matching by identical position and closest ", Method, ".")}

    Position_Match <- full_join(Ref, Mat, by="Position", keep=TRUE, relationship="many-to-many", suffix=c(".Ref", ".Mat")) %>%
      filter(.data$Position.Ref==.data$Position.Mat & .data$Name.Ref!=.data$Name.Mat) %>%
      mutate(Distance=abs(.data$Position.Ref-.data$Position.Mat),
             LRR_sd_delta=abs(.data$LRR_sd.Ref-.data$LRR_sd.Mat),
             LRR_mean_delta=abs(.data$LRR_mean.Ref-.data$LRR_mean.Mat),
             BAF_delta=abs(.data$BAF.Ref-.data$BAF.Mat),
             Position=ifelse(.data$Position.Ref==.data$Position.Mat, "Same", "Different"),
             Name=ifelse(.data$Name.Ref==.data$Name.Mat, "Same", "Different"))

    # PULL BEST MATCHED PROBES FOR IDENTICAL POSITION ON METHOD
    while(nrow(Position_Match) > 0) {
      # SELECT BEST PROBE METHOD-WISE
      Select_PosMatch <- Position_Match %>%
        group_by(.data$Name.Ref) %>%
        arrange_at({{Method}}) %>%
        filter(row_number()==1) %>%
        ungroup()

      # MERGE WITH MASTER DF
      DF <- rbind(DF, Select_PosMatch)

      # PURGE MATCHED
      Mat_Remove <- intersect(Mat$Name, Select_PosMatch$Name.Mat)
      Ref_Remove <- intersect(Ref$Name, Select_PosMatch$Name.Ref)

      # PURGE MATCHED FROM POSITION MATCH
      Position_Match <- Position_Match %>%
        filter(!.data$Name.Mat %in% Mat_Remove) %>%
        filter(!.data$Name.Ref %in% Ref_Remove)
    }

    # PURGE MATCHED FROM INPUTS
    Mat <- Mat %>%
      filter(!.data$Name %in% DF$Name.Mat)

    Ref <- Ref %>%
      filter(!.data$Name %in% DF$Name.Ref)

    # DIFFERENT POSITION DIFFERENT NAMES
    if(Verbose==TRUE){message("Matching SNPs within ", D_MAX, "bp and closest ", Method, ".")}
    St = seq(min(Reference$Position), max(Reference$Position), by=(Chunk-D_MAX))
    En = pmin(St + Chunk - 1, max(Reference$Position))
    Seq <- data.frame(S = St, E = En)

    for(i in 1:nrow(Seq)) {
      if(Verbose==TRUE){message(Chunk, "bp sliding window ", i, " out of ", nrow(Seq), ".")}
      T_REF <- Ref %>%
        filter(.data$Position >= Seq$S[i] & .data$Position <= Seq$E[i])
      T_MAT <- Mat %>%
        filter(.data$Position >= Seq$S[i] & .data$Position <= Seq$E[i])

      Match <- cross_join(T_REF, T_MAT, suffix=c(".Ref", ".Mat")) %>%
        mutate(Distance=abs(.data$Position.Ref-.data$Position.Mat)) %>%
        filter(.data$Distance<=D_MAX) %>%
        mutate(LRR_sd_delta=abs(.data$LRR_sd.Ref-.data$LRR_sd.Mat),
               LRR_mean_delta=abs(.data$LRR_mean.Ref-.data$LRR_mean.Mat),
               BAF_delta=abs(.data$BAF.Ref-.data$BAF.Mat),
               Position=ifelse(.data$Position.Ref==.data$Position.Mat, "Same", "Different"),
               Name=ifelse(.data$Name.Ref==.data$Name.Mat, "Same", "Different"))

      while(nrow(Match) > 0){
        # SELECT BEST PROBE METHOD-WISE
        Select_Match <- Match %>%
          group_by(.data$Name.Ref) %>%
          arrange_at({{Method}}) %>%
          filter(row_number()==1) %>%
          ungroup() %>%
          group_by(.data$Name.Mat) %>%
          arrange_at({{Method}}) %>%
          filter(row_number()==1) %>%
          ungroup()

        # MERGE WITH MASTER DF
        DF <- rbind(DF, Select_Match)

        # PURGE MATCHED PROBES
        Mat_Remove <- intersect(Mat$Name, Select_Match$Name.Mat)
        Ref_Remove <- intersect(Ref$Name, Select_Match$Name.Ref)

        # PURGE MATCHED
        Match <- Match %>%
          filter(!.data$Name.Mat %in% Mat_Remove) %>%
          filter(!.data$Name.Ref %in% Ref_Remove)
      }

      # PURGE MATCHED FROM INPUT
      Mat <- Mat %>%
        filter(!.data$Name %in% DF$Name.Mat)
      Ref <- Ref %>%
        filter(!.data$Name %in% DF$Name.Ref)
    }

    # FINAL CLEANUP
    DF <- DF %>%
      group_by(.data$Name.Ref) %>%
      arrange_at({{Method}}) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      group_by(.data$Name.Mat) %>%
      arrange_at({{Method}}) %>%
      filter(row_number()==1) %>%
      ungroup()
  }

  # SUMMARIZE RESULTS AND SAVE
  if(Verbose==TRUE){message("Saving overall dataframe to ", OutPath, ".csv")}
  DF %>% write_delim(paste0(OutPath, ".tsv"), col_names=TRUE, delim="\t")

  if(Verbose==TRUE){message("Saving Reference SNP name select to ", OutPath, "_RefSelect.csv")}
  DF %>% select(.data$Name.Ref) %>% write_csv(paste0(OutPath, "_RefSelect.txt"), col_names=FALSE, eol="\n")

  if(Verbose==TRUE){message("Saving Matching SNP name select to ", OutPath, "_MatSelect.csv")}
  DF %>% select(.data$Name.Mat) %>% write_csv(paste0(OutPath, "_MatSelect.txt"), col_names=FALSE, eol="\n")
  return(DF)

  end_time <- Sys.time()
  total_time <- round(difftime(end_time, start_time, units='mins'), 4)
  message("MarkerMatch complete in ", total_time," minutes.\n")
}

# -----------------------------------------------------------------------------
# End of Marker_Match.R
# -----------------------------------------------------------------------------
