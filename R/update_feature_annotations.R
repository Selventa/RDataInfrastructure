###################################################################################################=
## Similar to the "update_annotations.R" file, this file houses functions for 
##  updating platform feature annotation files (GPLs) from their default state 
##  (as downloaded from GEO).
###################################################################################################=


## The "Illumina human-6 v2.0 expression beadchip (extended)" chip.
UpdateFeatureAnnotations_GPL6370 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL6370 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    ## Load biomaRt to be able to translate the nuIDs (Illumina's "nucleotide 
    ##  universal identifiers") into EGID.
    # library(biomaRt)
    # ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
    
    ## Use Lumi to map the nuIDs to EGIDs
    require(lumi)
    require(lumiHumanIDMapping)
    nuID.to.EGID.mapping <- nuID2EntrezID(nuID = NULL, 
                                          lib.mapping = 'lumiHumanIDMapping', 
                                          returnAllInfo = TRUE)
    temp.idx <- match(annotations.df$NuID, rownames(nuID.to.EGID.mapping))
    updated.columns <- cbind(nuID.to.EGID.mapping[temp.idx, c("EntrezID", "Symbol")], 
                             Date=as.character(Sys.Date(), format="%b %d %Y"))
    colnames(updated.columns) <- c("EGID", "GeneSymbol", "EGID_and_GeneSymbol_added_date")
    
    ## Add the EGIDs and the current date onto the annotation data frame.
    updated.annotation.df <- cbind(annotations.df, updated.columns)
    
    ## Write out the updated annotation data frame back to CSV (replacing
    ##  the original annotation file if specified).
    if (replace.original) {
      write.csv(updated.annotation.df, file = paste0(gpl.path, "/", gpl.filename), row.names = FALSE)
    } else {
      write.csv(updated.annotation.df, 
                file = paste0(gpl.path, "/", 
                              gsub(pattern = ".csv$", 
                                   replacement = "_updated.csv", 
                                   gpl.filename, 
                                   ignore.case = TRUE)), 
                row.names = FALSE)
    }
    
    cat("Updated the GPL6370 platform annotation file to include EGIDs and Gene Symbols.\n")
  } else {
    cat("The GPL6370 platform annotation file didn't require updating.\n")
  }
}




