library(affyio)
library(affy)

suppressWarnings(setInternet2(use = FALSE))

# Will load from the local file system a list of esets matching the criteria specified
# in expt.annot.
# If from_raw is requested, but only from_processed exists for that GSE.ID, then
# the from_processed data will be returned (with appropriate warnings).  However, 
# if from_raw is requested and a from_raw output exists for a different expt.annot
# configuration, then no eset will be returned for that GSE.ID.  This way we can
# be specific about the data that we process from raw (since there may be a few
# ways of processing it), but if raw data is not available then we can always get
# back the processed data.  If allow.multiple.matches.per.file==F then in cases
# where multiple esets in the same file (ie, same GSE/GPL combination) match 
# the given expt.annot then none will be returned and a warning will be generated.
# This will allow us to identify cases where we're essentially writing duplicate
# or ambiguously annotated esets.
LoadEset <- function(GSE.ID, eset.folder = "~/esets", 
                          expt.annot=list(data.source = "from_raw",
                                       brainarray=T),
                     verbose=T,
                     allow.multiple.matches.per.file=F) {
  
  # if eset.folder is NA, then return an empty list
  if (is.na(eset.folder)) {
    return(invisible(list()))
  }
  
  RData.files <- list.files(path = eset.folder,
                            pattern = "\\.RData$")
  
  # for each GSE.ID given, look for a corresponding .RData file.  If it exists
  # then open it.  Any esets in the file that match the expt.annot are added to 
  # the list.  No warning is given for cases where the eset is not available
  esets  <- list()
  for (cur.GSE.ID in GSE.ID) {
    matching.files <- grep(cur.GSE.ID, RData.files, value=T)
    
    
    # loop through all matches, load the file
    for (cur.file in matching.files) {
      tmp.env <- new.env()
      load(file.path(eset.folder, cur.file), envir = tmp.env)

      
      # first see if we're asking for from_raw and only from_processed exists.  
      # If this is the case, then add the from_processed, throw a warning, and 
      # move to the next file.  If it's not the case, then try to match on all
      # relevant expt.annot fields
      loaded.data.source.types <- sapply(tmp.env$esets,
                                         function(cur.eset) {
                                           notes(cur.eset)$data.source
                                         })      
      loaded.only.from_processed <- (length(loaded.data.source.types)==1 &&
                                       all(loaded.data.source.types=="from_processed"))
      if (expt.annot$data.source=="from_raw" & loaded.only.from_processed) {
        msg <- paste0("Eset processed from raw data is not available for ",
                      cur.GSE.ID, " in ", cur.file, ". Loading from_processed eset.")
        cat(msg, "\n")
        warning(msg)
      }
      
      cur.file.esets <- list()
      # see if any of the esets match the given criteria
      for (cur.tmp.eset in tmp.env$esets) {
        # for those expt.annot fields that are captured in the current eset (ie
        # those that are relevant and were inserted by the data processing function)
        # look for a match
        tmp.expt.annot <- expt.annot[names(expt.annot) %in% names(notes(cur.tmp.eset))]
        expt.annot.matches <- 
          identical(notes(cur.tmp.eset)[names(tmp.expt.annot)],
                    tmp.expt.annot)
        
        if (all(expt.annot.matches)) {
          cur.file.esets <- c(cur.file.esets, list(cur.tmp.eset))
        }
      }
      
      # if allow.multiple.matches.per.file==F and we found more than one match
      # throw a warning and don't load them into the eset that is returned.
      if (!allow.multiple.matches.per.file & length(cur.file.esets)>1) {
        warning(paste0("Multiple matching esets were found in '", cur.file, 
                       "'.  None were loaded because allow.multiple.matches.per.file=F"))
      } else {
        esets <- c(esets, cur.file.esets)
      }
      
      rm(tmp.env)
    }
  }

  names(esets) <- lapply(esets, GetEsetName)
  
  return(invisible(esets))
}


SaveEset <- function(cur.eset, eset.folder) {
  cur.eset.name <- GetEsetName(cur.eset)
  
  save.file <- file.path(eset.folder, paste0(cur.eset.name, ".RData"))
  
  
  if (file.exists(save.file)) {
    load(save.file)
  } else {
    esets <- list()
  }
  
  # see if any of the loaded esets have the same notes.  Ignore original.pData 
  # since this _may_ be modified under some circumstances (like when samples
  # are ommited)
  loaded.eset.to.replace <- 
    sapply(esets,
           function(cur.loaded.eset) {
             cur.loaded.eset.notes <- notes(cur.loaded.eset)
             cur.eset.notes <- notes(cur.eset)
             
             cur.loaded.eset.notes$original.pData <- 
               cur.eset.notes$original.pData <- 
               NULL
             
             identical(cur.loaded.eset.notes,
                       cur.eset.notes)
           })
  if (any(loaded.eset.to.replace)) {
    stopifnot(sum(loaded.eset.to.replace)==1)
    esets[loaded.eset.to.replace] <- cur.eset
  } else {
    esets <- c(esets, cur.eset)
  }
  
  
  save(esets, file=save.file)
  
  return(invisible(NULL))
}


GetEsetName <- function(cur.eset) {
  # get names out of the notes(cur.eset)@$name slot
  return(notes(cur.eset)$name)
}

# This function will either load the esets corresponding to GSE.ID from local files 
# (if they're there already and if overwrite.existing==F), or will go download
# the data from GEO and process it according to expt.annot if desired. The 
# annotation files are automatically updated for any data downloaded from GEO, 
# but not for any data already in local files.

##what other argument options are available for expt.annot?
GetEset <- function(GSE.ID, eset.folder = "S:/Groups/R and D Group Documents/GEO_data/esets/", 
                    overwrite.existing=F,
                    cache.folder = normalizePath("~/../Downloads"),
                     expt.annot=list(data.source = "from_processed"),
                    annot.csv.folder="S:/Groups/R and D Group Documents/GEO_data/sample_annotations/",
                    feature.annotation.path="S:/Groups/R and D Group Documents/GEO_data/feature_annotations/",
                    verbose=T) {
 
  library("GEOquery")	
 
  if (!("data.source" %in% names(expt.annot))) {
    stop("'expt.annot' must contain an element named 'data.source'.")
  }

  if (!(expt.annot$data.source %in% c("from_processed", "from_raw"))) {
    stop("expt.annot$data.source must be 'from_processed' or 'from_raw'.")
  }
  
  # if we're loading processed data, strip away any other annotations because 
  # they won't apply
  if (expt.annot$data.source=="from_processed") {
    expt.annot <- expt.annot["data.source"]
  }
  
  esets <- list()
  
  for (cur.GSE.ID in GSE.ID) {
    if (!overwrite.existing) {
      cur.esets <- LoadEset(cur.GSE.ID, 
                                 eset.folder = eset.folder, 
                                 expt.annot=expt.annot)
      
      # if we managed to load from file we can skip to the next GSE.ID
      if (length(cur.esets)>0) {
        esets <- c(esets, cur.esets)
        if (verbose) {cat("Loaded ", sub("from_", "", expt.annot$data.source), 
                          " data for ", cur.GSE.ID, " from file.\n", sep="")}
        next
      }
    }
    
    
    if (expt.annot$data.source=="from_processed") {
      if (verbose) {cat("Downloading processed data for ", cur.GSE.ID, ".\n", sep="")}
      
      # call function from geoQuery package to get processed data from GEO
      cur.esets <- getGEO(cur.GSE.ID, destdir=cache.folder, GSEMatrix=T, getGPL=F)

      # label each as coming from processed data
      cur.esets <- lapply(cur.esets,
                          function(cur.eset) {
                            notes(cur.eset)$data.source <- "from_processed"
                            return(cur.eset)
                          })

      # cur.esets will be a list of eset objects, one for each platform present 
      # in this geo entry. Rename them based on the platform
      platforms <- sapply(cur.esets,
                          function(x) as.character(annotation(x)))
      names(cur.esets) <- paste0(cur.GSE.ID, "_", platforms)
    
      # store the eset name, and other processing information, in the eset notes
      for (cur.eset.name in names(cur.esets)) {
        cur.eset <- cur.esets[[cur.eset.name]]
        notes(cur.eset)$name <- cur.eset.name
        notes(cur.eset)$GSE.ID <- cur.GSE.ID
        notes(cur.eset)$platform <- as.character(annotation(cur.eset))
        notes(cur.eset)$data.source <- "from_processed"
        notes(cur.eset)$original.pData <- pData(cur.eset)
        
        cur.esets[[cur.eset.name]] <- cur.eset    
      }
      
    } else if (expt.annot$data.source=="from_raw") {
      
      # if we want to process the raw data, first get the pre-processed data (this is
      # where the annotations will come from)
      cur.preprocessed.esets <- GetEset(cur.GSE.ID, 
                                        eset.folder = eset.folder, 
                                        overwrite.existing=overwrite.existing,
                                        cache.folder = cache.folder,
                                        expt.annot = list(data.source="from_processed"),
                                        annot.csv.folder=annot.csv.folder,
                                        feature.annotation.path=feature.annotation.path,
                                        verbose=verbose) 
      if (verbose) {cat("Downloading and processing raw data for ", cur.GSE.ID, ".\n", sep="")}
      
      # now get processed data (will use pData from cur.preprocessed.esets)
      cur.esets <- lapply(cur.preprocessed.esets,
                      ProcessRawGEOData,
                      cache.folder, 
                      expt.annot=expt.annot,
                      verbose=verbose)
      
      # note that the relevant notes fields for each eset are added by  
      # ProcessRawGEOData() (for raw data processing).  
      
    }
    


    

    # we update annotations here so annotations get updated before we save to file
    cur.esets <- lapply(cur.esets,
                    UpdateAnnotations,
                    annot.csv.folder=annot.csv.folder)
    
    # map features to EGID if possible
    cur.esets <- lapply(cur.esets,
                        map.features.to.EGID,
                        feature.annotation.path=feature.annotation.path,
                        cache.folder=cache.folder)
    
    # save to file here instead of outside the loop so that when an eset is loaded
    # from an .RData file we aren't wasting time writing it back (because next()
    # is used to go to the next iteration in the loop)
    if (!is.na(eset.folder)) {
      lapply(cur.esets, 
             SaveEset,
             eset.folder)
    }
    
    esets <- c(esets, cur.esets)
  }
    
    
    
  

  return(invisible(esets))
}




map.features.to.EGID <- function(cur.eset, platform=NULL, feature.annotation.path="~/Datasets/Feature annotation files/", annotation.file=NULL, cache.folder="~/../Downloads") {
  
  # if we already have EGID in the feature annotations, then return immediately
  if ("EGID" %in% colnames(fData(cur.eset))) {
    return(cur.eset)
  }
  
  if (is.null(platform)) {
    platform <- annotation(cur.eset)
  }
  
  # open annotation file from common location
  if (is.null(annotation.file)) {
    annotation.file <- file.path(feature.annotation.path, paste0(platform, ".csv"))
  }
  
  if (!file.exists(annotation.file)) {
    annotation.file <- getGPLFile(cur.GPL=platform, 
                                  feature.annotation.path=feature.annotation.path,
                                  cache.folder=cache.folder)
    if (is.null(annotation.file)) {
      warning("Features are not mapped to EGIDs for ", 
              GetEsetName(cur.eset))
      return(cur.eset)
    }
  }
  
  annot <- read.csv(annotation.file, stringsAsFactors=F)
  
  if (!("EGID" %in% colnames(annot))) {
    warning("EGID column does not exist in feature annotation file ",
         annotation.file, ".  Cannot map features to EGIDs.")
    return(cur.eset)
  }
  
  if (!all(rownames(cur.eset) %in% annot[[1]])) {
    num.in.annot <- sum(rownames(fData(cur.eset)) %in% annot[[1]])
    warning("Only ", num.in.annot, " of ", nrow(cur.eset), 
            " features are found in the feature annotation file ",
            annotation.file, 
            ".  ", nrow(cur.eset) - num.in.annot, 
            " features will be discarded.")
    cur.eset <- cur.eset[rownames(cur.eset) %in% annot[[1]], ]
  }
  
  fData(cur.eset)$EGID <- annot$EGID[match(rownames(cur.eset), annot[[1]])] 
  
  return(cur.eset)
}



getGPLFile <- function(cur.GPL, feature.annotation.path, cache.folder = normalizePath("~/../Downloads")) {
  
  GPL.obj <- tryCatch({
    getGEO(cur.GPL, destdir=cache.folder)
  },
  error=function(err) {
    print(paste("Error: ", err))
    cat("A error occurred while downloading '", cur.GPL, "'.\n", sep="")
    return(NULL)
    
  })
  
  if (is.null(GPL.obj)) return(NULL)

  feat.annot.file <- file.path(feature.annotation.path,
                               paste0(cur.GPL, ".csv"))
  
  updated.date <- Meta(GPL.obj)$last_update_date
  
  feat.annot <- GEOquery::Table(GPL.obj)
  feat.annot$last_updated_date <- updated.date
  
  EGID.col <- grep("EGID|ENTREZ_GENE_ID", colnames(feat.annot))

  if (length(EGID.col)!=1) {
    cat("Unable to determine EGID column from ", cur.GPL, ". Please manually modify ",
        feat.annot.file, " by labeling the EGID column.\n", sep="")
                  
    write.csv(feat.annot, 
              file=feat.annot.file,
              row.names=F)
    
    return(NULL)
  } 
  feat.annot$EGID <- feat.annot[[EGID.col]]
  write.csv(feat.annot, 
            file=feat.annot.file,
            row.names=F)
  
  return(feat.annot.file)  
}


UpdateAnnotations <- function(cur.eset, annot.csv.folder) {

  cur.eset.name <- notes(cur.eset)$name
  
  annot.update.func.name <- paste0("UpdateAnnotations_", cur.eset.name)
  annot.update.csv <- file.path(annot.csv.folder,
                                paste0(cur.eset.name, ".csv"))
  # set the tmp.eset to NULL so we can easily tell if we call a function to update
  # the annotations
  tmp.eset <- NULL
  if (exists(annot.update.func.name)) {
    tmp.eset <- eval(parse(text=paste0(annot.update.func.name,
                                       "(cur.eset)")))
    if (file.exists(annot.update.csv)) {
      warning(paste0("Using function ", annot.update.func.name, "() instead of", 
                     " file ", annot.update.csv, " to update annotations."))
    }
  } else if (file.exists(annot.update.csv)) {
    tmp.eset <- UpdateAnnotations_CSV(cur.eset, annot.update.csv)
  }
  
  if (!is.null(tmp.eset)) {
    # make sure that the number of rows in pData(tmp.eset) and 
    # notes(tmp.eset)$original.pData are the same as the number of columns/samples 
    # in tmp.eset
    if (nrow(pData(tmp.eset)) != nrow(notes(tmp.eset)$original.pData)) {
      warning(paste0(annot.update.func.name, "() did not maintain the same number",
                     " of rows in pData and original.pData.  Annotations were",
                     " not updated"))
    } else if(nrow(pData(tmp.eset)) != ncol(tmp.eset)) {
      warning(paste0(annot.update.func.name, "() did not maintain the same number",
                     " of samples in pData the eset.  Annotations were",
                     " not updated"))
    } else {
      cur.eset <- tmp.eset
    }
    
    
    # if "supplementary_file" in original.pData but not in pData(), add to pData
    # - this is important for processRawGeoData() because it will try to get
    #   sup files from pData, not original.pData, so that if one of the filenames
    #   is wrong in original.pData it can be corrected and the correct name will
    #   be used.
    if (("supplementary_file" %in% names(notes(cur.eset)$original.pData)) &
          !("supplementary_file" %in% names(pData(cur.eset)))) {
      pData(cur.eset)$supplementary_file <- 
        notes(cur.eset)$original.pData$supplementary_file
    }
    
    # if "geo_accession" in original.pData but not in pData(), add to pData
    # - this is important for processRawGeoData() because it will try to get
    #   geo_accession from pData, not original.pData
    if (("geo_accession" %in% names(notes(cur.eset)$original.pData)) &
          !("geo_accession" %in% names(pData(cur.eset)))) {
      pData(cur.eset)$geo_accession <- 
        notes(cur.eset)$original.pData$geo_accession
    }
    

  } else {
    cat("\nNo annotation processing file/function exists for ", cur.eset.name, ".",
        "  Saving phenoData to ", annot.update.csv, "\n",sep="")
    cur.eset.name <- GetEsetName(cur.eset)
    write.csv(pData(cur.eset), 
              file=annot.update.csv,
              row.names=F)
  }
  return(cur.eset)
}



# # generic function for downloading data files from the internet
# downloadData <- function(path.to.data, name, folder=NULL, overwrite=FALSE, md5=NULL) {
#   
#   # TODO
#   # put in update functionality - how to determine file size?
#   # determine whether data is zipped or tar.gz and unpack it
#   # how to associate annotation with this file?
#   
#   if(!is.null(md5) && file.exists(md5)){ # error without shortcircuit
#     md5 <- read.delim2(md5)
#   }
#   
#   # check if file exists - if it does then read it and download data
#   if(file.exists(path.to.data)){
#     print("Assuming input is a text file of URLs to download")
#     # read text file
#     urls <- read.delim2(path.to.data)
#     
#     for(i in 1:length(urls)){
#       fname <- unlist(strsplit(url[i], "/", fixed=TRUE))
#       fname <- fname[length(fname)]
#       
#       # if we should not overwrite file, then check to see if it exists in the current location
#       if(!overwrite){
#         logic <- file.exists(paste0(ifelse(is.null(folder), getwd(), folder), "/", fname))
#         if(logic) stop("Stopping because overwrite=FALSE and the file already exists")
#       }
#       download.file(url[i], destfile=paste0(ifelse(is.null(folder), getwd(), folder), "/", fname))
#       
#       if(!is.null(md5)){
#         stopifnot(md5sum(paste0(ifelse(is.null(folder), getwd(), folder), "/", fname, ))==md5[i])
#       }
#     }
#   } else {
#     print("Assuming input is URL to download")
#     fname <- unlist(strsplit(path.to.data, "/", fixed=TRUE))
#     fname <- fname[length(fname)]
#     
#     # if we should not overwrite file, then check to see if it exists in the current location
#     if(!overwrite){
#       logic <- file.exists(paste0(ifelse(is.null(folder), getwd(), folder), "/", fname))
#       if(logic) stop("Stopping because overwrite=FALSE and the file already exists")
#     }
#     
#     download.file(path.to.data, destfile=paste0(ifelse(is.null(folder), getwd(), folder), "/", fname))
#     
#     if(!is.null(md5)){
#       stopifnot(md5sum(paste0(ifelse(is.null(folder), getwd(), folder), "/", fname, ))==md5)
#     }
#   }
#   return(NULL)
# }

