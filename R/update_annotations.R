# This file contains dataset-specific functions for modifying the pData object
# associated with an eset.  The fucntions are named "UpdateAnnotations_GSEXXXXX_GPLYYYYY",
# and are identified by matching each eset by it's name (which is "GSEXXXXX_GPLYYYYY")
# to the appropriate function name.  These functions should not be called directly, 
# rather the eset should be passed to UpdateAnnotations() which will dispatch to the
# correct dataset-specific function herein, and do some simple error checking on
# the eset that is returned.

# Each function must take a single input cur.eset.  The function should access the
# original pData object stored in notes(cur.eset)$original.pData, and use this as
# the basis for making a new pData object which is then placed in pData(cur.eset).
# By always starting from the original pData, then the function can be called
# mulitple times sequentially and produce identical results (which may not
# be the case if it were to start with pData(eset)).  notes(cur.eset)$original.pData
# should generally not be modified (see below).  cur.eset is return from the
# functions.

# if supplementary_file and geo_accession columns exist in original.pData then
# they will be automatically copied to pData(cur.eset) by UpdateAnnotations() if
# they are not put there by the dataset-sepcific function below.

# There may be cases where we want the UpdateAnnotations function to actually
# remove samplse from the eset. No obvious use cases for this come to mind, but
# based on the current design this should be possible without generating an error
# elsewhere.  However, this should be done judiciously.  If this is done, then
# the function must subset pData and notes(cur.eset)$original.pdata to also
# remove those samples.  




get.interesting.annot.cols <- function(annot) {
  ##help! i am a function with no comments!!
  unique.vals <- sapply(annot, function(cur.col) {length(unique(cur.col))})
  
  annot.summary <- c()
  for (col.num in which(unique.vals>1)) {
    cur.col <- annot[[col.num]]
    col.name <- colnames(annot)[col.num]
    if (is.numeric(cur.col) & length(unique(cur.col))>(nrow(annot)/3)) {
      annot.summary[[col.name]] <- summary(cur.col)
    } else if (!is.numeric(cur.col) & length(unique(cur.col))>(nrow(annot)/3)) {
      annot.summary[[col.name]] <- c(as.character(head(cur.col)), "...")
    } else if (is.numeric(cur.col)) {
      annot.summary[[col.name]] <- table(cur.col)
    } else {
      annot.summary[[col.name]] <- table(cur.col)
    }
    
  }
  return(annot.summary)
}


UpdateAnnotations_CSV <- function(cur.eset, annot.csv) {
  if (!file.exists(annot.csv)) {
    return(NULL)
  }
  
  new.annot <- read.csv(annot.csv, stringsAsFactors=F)
  rownames(new.annot) <- new.annot$geo_accession
  
  annot <- notes(cur.eset)$original.pData

  if (!all(rownames(new.annot) %in% rownames(annot))) {
    warning("Annotation csv ", annot.csv, 
            " contains samples missing from the eset.  ",
            "Annotations are not updated for ", notes(cur.eset)$name)
    return(NULL)
  }
    
  if (!all(rownames(annot) %in% rownames(new.annot))) {
    num.missing <- sum(!(rownames(annot) %in% rownames(new.annot))) 
    warning(paste0(num.missing, " samples from ", notes(cur.eset)$name, 
                   " were excluded in the annotation csv ",
                   annot.csv, ".  These samples will be dropped."))
    cur.eset <- cur.eset[,rownames(new.annot)]
    annot <- annot[rownames(new.annot),]
    notes(cur.eset)$original.pData <- annot
    
  }
  
  new.annot <- new.annot[rownames(annot), ,drop=F]
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}

