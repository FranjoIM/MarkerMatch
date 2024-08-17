ManifSum <- function(fPath=fPath, oPath=oPath, nSam=nSam){
    # Require packages
    require(tidyverse)
    require(matrixStats)
  
    # Count number of lines in the LRR file
    f <- file(fPath, open="rb")
    nlines <- 0L
    while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
        nlines <- nlines + sum(chunk == as.raw(10L))
    }
    close(f); rm(f, chunk)
    
    # Initialize dataframe and counter variables
    DF <- data.frame()
    s <- 1
    nMax <- 20000
    
    # Loop through the file, calculating the metrics
    while(s < nlines){
        DAT <- read_delim(fPath,
            col_names=FALSE,
            col_types=c(rep("c",2), rep("d",(nSam+3))),
            na = c("", "NA", "NaN"),
            skip=s,
            n_max=nMax)
        
        MOD <- DAT %>%
        mutate(LRR_mean=rowMeans(pick(where(is.numeric), -c(X1, X2, X3)), na.rm=TRUE),
            LRR_sd=rowSds(as.matrix(pick(where(is.numeric), -c(X1, X2, X3))), na.rm=TRUE),
            .keep="unused")
        
        DF <- rbind(DF, MOD)
        
        message("Processed rows ", s," to ", s+nMax, ". ", round(s/nlines*100, digits=2), "% complete.")
        s <- s + nMax
    }
    message("Processed all rows. Saving file", round((s+nMax)/nlines*100, digits=2), "% complete.")

    # Rename variables 
    DF <- DF %>%
        rename(Name=X1, Chr=X2, Position=X3)
    
    # Export file
    DF %>%
        write_csv(oPath)

    message("Process complete.")
}
