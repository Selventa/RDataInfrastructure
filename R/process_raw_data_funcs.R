
library(affyio)
library(affy)


ProcessRawGEOData <- function(cur.eset, cache.folder, expt.annot, verbose=T) {
  ##should use getesetname here and/or source load_data_funcs?
  cur.eset.name <- notes(cur.eset)$name
  # get the pData from the eset, so we can keep use the same annotations 
  # for the raw data
  cur.pData <- pData(cur.eset)
  orig.pData <- notes(cur.eset)$original.pData
  
  # get the experiment notes from the eset
  expt.data <- notes(cur.eset)
  
  platform <- expt.data$platform
  
  if (!("supplementary_file" %in% colnames(orig.pData))) {
    if (verbose) {cat("Raw data not available for '", cur.eset.name, "'. ",
                      "Using processed data for this data set.\n", sep="")}
    warning("Raw data not available for '", cur.eset.name, "'. ",
            "Using processed data for this data set.")
    return(cur.eset)
  }

  if (any(as.character(orig.pData$supplementary_file) == "NONE")){
    warning("Some raw data is missing from GEO for ", cur.eset.name, ".
            Files from ", toString(orig.pData$geo_accession[as.character(orig.pData$supplementary_file) == "NONE"]),
            " are at fault. If desired, change the update_annotation function for", cur.eset.name,
            " to remove these samples, and then the raw data can be processed.")
  }
  
  cur.ds.processing.func.name <- paste0("process.data.", cur.eset.name)
  
  # if there is a processing function for this specific dataset, then use that.
  # otherwise look for a processing function for the current platform
  if (exists(cur.ds.processing.func.name)) {
    cur.processing.func.name <- cur.ds.processing.func.name
  } else {
    # get the names of the supp files to see if they are affy cel files
    supp.files <- as.character(orig.pData$supplementary_file)
    
    is.affy.file <- all(grepl("\\.cel$|\\.cel\\.", supp.files, ignore.case = T))
    
    if (is.affy.file) {
      cur.platform.processing.func.name <- "process.data.affy" 
    } else {
      cur.platform.processing.func.name  <- paste0("process.data.", platform)
    }
    
    if (!exists(cur.platform.processing.func.name)) {
      ##this print statement is not entirely accurate?
      if (verbose) {cat("Raw data not available for '", cur.eset.name, "'. ",
                        "Using processed data for this data set.\n", sep="")}
      warning("Raw data processing function not available for '", 
              cur.eset.name, "'. ",
              "Using processed data for this data set.")
      return(cur.eset)
    } else {
      cur.processing.func.name <- cur.platform.processing.func.name
    }
  }
  
  # this is where getGEOSuppFiles should save the files.  So check to see if
  # they exist first to avoid re-downloading something that we already have
  cur.gsm.files <- file.path(cache.folder, orig.pData$geo_accession, basename(as.character(orig.pData$supplementary_file)))
  
  # for each file, see if we have it downloaded already.  If not, download it.
  downloaded.files <- 
    sapply(1:nrow(orig.pData),
           function(i) {
             if (file.exists(cur.gsm.files[i])) {
               return(cur.gsm.files[i])
             } else {
               return(rownames(getGEOSuppFiles(orig.pData$geo_accession[i], 
                                               baseDir=cache.folder)))
             }
           })
  # if it downloaded multiple files per sample, this will be a list.  unlist
  # it to get a vector of downloaded file names
  downloaded.files <- unlist(downloaded.files)
  
  # check to make sure that the expected gsm files are in the downloaded files.
  # If not, may need to update annotation file before running this function to
  # correct supp file names.  Can't just take this list of downloaded files
  # as the files to process because sometimes there are extraneous files (like processed data)...
  if (!all(cur.gsm.files %in% downloaded.files)) {
    #if there is a simple naming problem, we might be able to solve it here
    
    # first see if there there is one downloaded file per gsm number
    downloaded.gsm.files <- sapply(as.character(orig.pData$geo_accession),
                                   function(cur.gsm) {
                                     # find all downloaded files that have the current gsm
                                     matches <- grep(cur.gsm, basename(downloaded.files), ignore.case = T)
                                     if (length(matches)!=1) {
                                       return(NA)  # return NA if there isn't exactly one match
                                     } else {
                                       return(downloaded.files[matches])
                                     }
                                   })
    # if any are NA, then the files don't match up 1-1 with gsm numbers, so
    # can't process raw data
    if (any(is.na(downloaded.gsm.files))) {
      if (verbose) {cat("Raw data not available for '", cur.eset.name, "'. ",
                        "Using processed data for this data set.\n", sep="")}
      warning("For ", cur.eset.name,
              ", downloaded raw data files don't match with files names from ",
              " annotation file, and can't be matched uniquely from GSM IDs.  ",
              "Using processed data.")
      return(cur.eset)
    }
    cur.gsm.files <- downloaded.gsm.files
  }

  # just in case, run normalizePath() to clean up the file paths
  cur.gsm.files <- normalizePath(cur.gsm.files)
  #check if files are larger than 0 kb, which implies that they are corrupted or downloaded improperly
  if (!all(sapply(cur.gsm.files, function(gsm.file) file.info(gsm.file)$size)>0)){
    warning("Some files from GEO downloaded improperly. Files from ",
            toString(orig.pData$geo_accession[sapply(cur.gsm.files, function(gsm.file) file.info(gsm.file)$size)==0]),
            " are at fault.If desired, change the update_annotation function for ", cur.eset.name,
            " to remove these samples, and then the raw data can be processed.")
    return(cur.eset)
  }
  #run original processing function
  #browser()
  cur.eset <- 
    eval(parse(text=paste0(cur.processing.func.name,
                           "(cur.gsm.files",
                           ifelse(is.affy.file, ",expt.annot=expt.annot", ""), 
                           ")")))
  # since we processed files based on the original pData, we want to set the
  # pData for the new eset to be the same as the original pData.
  pData(cur.eset) <- orig.pData
  notes(cur.eset)$original.pData <- orig.pData
  
  # identify the items from expt.data that aren't in the eset, and add them
  expt.data.to.add <- expt.data[!(names(expt.data) %in% 
                                    names(notes(cur.eset)))]
  notes(cur.eset) <- c(notes(cur.eset),
                       expt.data.to.add)
  # overwrite the data.source
  notes(cur.eset)$data.source <- "from_raw"

  
  # set the platform here because the processing function might not know/capture
  # the platform, and we don't want it to change from the platform listed for
  # the pre-processed data
  annotation(cur.eset) <- notes(cur.eset)$platform
  
  return(cur.eset)
  
}



process.data.affy <- function(data.files, expt.annot) {
  
  if (!("brainarray" %in% names(expt.annot))) {
    brainarray <- F
  } else {
    brainarray <- expt.annot["brainarray"]
  }
  
  # get the chip type
  cel.header <- read.celfile.header(data.files[1])

  # get the default cdf name, and clean it up a bit
  cdf.name <- cel.header$cdfName
  cdf.name <- gsub("[^[:alnum:]]","",cel.header$cdfName)
  
  # if it's an exon array, then we need to use brainarray
  if (grepl("HuGene|HuEx|MoGene|MoEx|RaGene|RaEx", cdf.name)) {
    
    brainarray <- T
    
    # in these cases we also need to modify the cdfName in the header slightly 
    # so we can match it with the approrpriate brainarray cdf
    cdf.name <- gsub("v1","",cdf.name)
    cdf.name <- gsub("v2","",cdf.name)
    
  }
    
  if (brainarray) {
    
    # map the cdf name to the brainarray cdf
    cdf.name <- paste(tolower(cdf.name),c("hs","mm","rn"),"entrezgcdf",sep="")
    
    # make sure that the brainarray cdf is installed
    avail.pkg <- installed.packages()[,"Package"]
    cdf.name <- intersect(cdf.name,avail.pkg)
    
    if(length(cdf.name) != 1) {
      warning("Cannot find Brainarray CDF for ",  tmp$cdfName, ".  Using default CDF.")
      cdf.name <- NULL
      brainarray <- F
    }
  } else {
    cdf.name <- NULL
  }
  
  cur.eset <- justRMA(filenames=data.files,
                      background=T,
                      normalize=T,
                      celfile.path=NULL, 
                      cdfname=cdf.name)
  
  # justRMA puts an empty value in notes(cur.eset).  Remove it.
  notes(cur.eset) <- list()
  
  if (brainarray) {
    fData(cur.eset)$EGID <- 
      suppressWarnings(as.numeric(sub("_at", "", rownames(cur.eset))))
  }
  
  notes(cur.eset)$brainarray <- brainarray
  
  return(cur.eset)
}
