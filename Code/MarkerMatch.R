MarkerMatch <- function(Reference=Reference, Matching=Matching, Method=Method, D_MAX=D_MAX, OutPath=OutPath){
  require(tidyverse)
  
  # Define errors for function
  TEST_COLNAME <- c("Name","Chr","Position","BAF","LRR_mean","LRR_sd")
  if(sum(intersect(names(Reference), TEST_COLNAME) %in% TEST_COLNAME)!=6){
    stop("Reference is missing one of the following required columns: Name, Chr, Position, BAF, LRR_mean, LRR_sd.")
  }
  
  if(sum(intersect(names(Matching), TEST_COLNAME) %in% TEST_COLNAME)!=6){
    stop("Matching is missing one of the following required columns: Name, Chr, Position, BAF, LRR_mean, LRR_sd.")
  }
  
  if(!Method %in% c("Distance","BAF_delta","LRR_mean_delta","LRR_sd_delta")){
    stop("Method must be one of the following: Distance, BAF_delta, LRR_mean_delta, or LRR_sd_delta.")
  }
  
  if(!is.numeric(D_MAX)){
    stop("D_MAX must be a numeric value.")
  }
  
  if(!is.character(OutPath)){
    stop("OutPath must be a string value containing a valid path (without file extensions).")
  }
  
  # Fetch the list of available chromosomes
  Chromosomes <- str_sort(unique(Reference$Chr), numeric=TRUE)
  
  # Create a holding dataframe
  DF <- data.frame()
  
  # Match markers
  for(chr in 1:length(Chromosomes)){
    message("Processing chromosome ", Chromosomes[chr], ".")
    
    # Filter out chromosome
    Ref <- Reference %>%
      drop_na() %>%
      filter(Chr==Chromosomes[chr])
    Mat <- Matching %>%
      drop_na() %>%
      filter(Chr==Chromosomes[chr])
    
    # Match manifests based on identical position, annotate distances, deltas, and name identity
    PerfectMatch <- full_join(Ref, Mat, by="Position", keep=TRUE, relationship="many-to-many", suffix=c(".Ref", ".Mat")) %>%
      drop_na() %>%
      mutate(Distance=abs(Position.Ref-Position.Mat),
             LRR_sd_delta=abs(LRR_sd.Ref-LRR_sd.Mat),
             LRR_mean_delta=abs(LRR_mean.Ref-LRR_mean.Mat),
             BAF_delta=abs(BAF.Ref-BAF.Mat),
             Position=ifelse(Position.Ref==Position.Mat, "Same", "Different"),
             Name=ifelse(Name.Ref==Name.Mat, "Same", "Different"))
    
    if(nrow(PerfectMatch)>0) {
      # Merge PerfectMatch with master DF
      DF <- rbind(DF, PerfectMatch)
      
      # Purge matched SNPs from Ref and Mat data
      Ref <- Ref %>%
        filter(!Name %in% PerfectMatch$Name.Ref) %>%
        arrange(Position)
      Mat <- Mat %>%
        filter(!Name %in% PerfectMatch$Name.Mat) %>%
        arrange(Position)
    }
    
    for(i in 1:nrow(Ref)){
      # Pull out SNP, filter matches by position distance criteria
      Match <- cross_join(Ref[i,], Mat, suffix=c(".Ref", ".Mat")) %>%
        mutate(Distance=abs(Position.Ref-Position.Mat)) %>%
        filter(Distance<=D_MAX)
      
      # Match manifests based on position, annotate distances, deltas, and name identity
      if(nrow(Match)>0){
        # Match manifests based on position, annotate distances, deltas, and name identity
        Match <- Match %>% 
          mutate(LRR_sd_delta=abs(LRR_sd.Ref-LRR_sd.Mat),
                 LRR_mean_delta=abs(LRR_mean.Ref-LRR_mean.Mat),
                 BAF_delta=abs(BAF.Ref-BAF.Mat),
                 Position=ifelse(Position.Ref==Position.Mat, "Same", "Different"),
                 Name=ifelse(Name.Ref==Name.Mat, "Same", "Different"))%>% 
          arrange_at({{Method}}) %>%
          filter(row_number()==1)
        
        # Merge Match with master DF
        DF <- rbind(DF, Match)
        
        # Purge matched SNPs from Ref and Mat data
        Mat <- Mat %>%
          filter(!Name %in% Match$Name.Mat)
        Ref <- Ref %>%
          filter(!Name %in% Match$Name.Ref)
      }
    }
    message("Completed chromosome ", Chromosomes[chr], ".\n")
  }
  # Summarize results, save, and return dataframe
  message("Saving overall dataframe to ", OutPath, ".csv")
  DF %>% write_csv(paste0(OutPath, ".csv"))
  
  message("Saving Reference SNP name select to ", OutPath, "_RefSelect.csv")
  DF %>% select(Name.Ref) %>% write_csv(paste0(OutPath, "_RefSelect.csv"), col_names=FALSE, eol="\n")
  
  message("Saving Matching SNP name select to ", OutPath, "_MatSelect.csv")
  DF %>% select(Name.Mat) %>% write_csv(paste0(OutPath, "_MatSelect.csv"), col_names=FALSE, eol="\n")
  return(DF)
  
  mesaage("Process complete.\n")
}
