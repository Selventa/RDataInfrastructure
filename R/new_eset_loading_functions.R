
Download_processed_from_GEO <- function(cur.GSE.ID,
                                        cache.folder=normalizePath("~/../Downloads")) {
  
  # call function from geoQuery package to get processed data from GEO
  cur.esets <- getGEO_retry(cur.GSE.ID, destdir=cache.folder, GSEMatrix=T, getGPL=F)
  
  # label each as coming from processed data
  cur.esets <- lapply(cur.esets,
                      function(cur.eset) {
                        notes(cur.eset)$data.source <- "from_processed"
                        return(cur.eset)
                      })
  
  
  cur.esets <- 
    lapply(cur.esets, 
           function(cur.eset) {
             cur.exprs <- exprs(cur.eset)
             if (max(cur.exprs, na.rm = T)>30) {
               cat("Assuming expression values are not logged.  Logging them now.\n")
               if (any(cur.exprs[!is.na(cur.exprs)] <= 0)) {
                 cat("Intensities less than or equal to zero are replaced with half of the smallest positive intensity.\n")
                 cur.exprs[!is.na(cur.exprs) & cur.exprs <=0] <- min(cur.exprs[!is.na(cur.exprs) & cur.exprs>0])/2
               }
               cur.exprs <- log2(cur.exprs)
               exprs(cur.eset) <- cur.exprs
             }
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
    
    # make sure that data is present in the eset.
    if (length(exprs(cur.eset))==0) {
      notes(cur.eset)$status <- "no data"
    } else {
      notes(cur.eset)$status <- "ok"
    }
    
    cur.esets[[cur.eset.name]] <- cur.eset    
    
  }
  
  return(cur.esets)
}

GetEset_new <- function(GSE.ID, eset.folder = "S:/Groups/R and D Group Documents/GEO_data/esets/", 
                    download.missing=T,
                    overwrite.existing=F,
                    cache.folder = normalizePath("~/../Downloads"),
                    expt.annot=list(data.source = "from_processed"),
                    annot.csv.folder="S:/Groups/R and D Group Documents/GEO_data/sample_annotations/",
                    feature.annotation.path="S:/Groups/R and D Group Documents/GEO_data/feature_annotations/",
                    remove.failed = T,
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
    
    cur.esets <- list()
    
    # if we're not automatically replacing existing data, try loading from_processed
    # data from file
    if (!overwrite.existing) {
      cur.esets <- LoadEset(cur.GSE.ID, 
                            eset.folder = eset.folder, 
                            expt.annot = expt.annot)

      cur.esets <- lapply(cur.esets,
                          UpdateAnnotations,
                          annot.csv.folder=annot.csv.folder)
      
      cur.esets <- lapply(cur.esets,
                          map.features.to.EGID,
                          feature.annotation.path=feature.annotation.path,
                          cache.folder=cache.folder)
      failed.esets <- sapply(cur.eset, function(cur.eset) notes(cur.eset)$status=="ok")
    }
    
    
    
    # we now check to see if we should download from_processed from GEO.  We download
    # it if overwrite existing is TRUE, or if we failed to load it from disk and
    # download.missing is TRUE
    
    if (overwrite.existing || 
        (length(cur.esets)==0 & download.missing) ||
        (any(failed.esets) & download.missing)) {
       
      if (expt.annot$data.source == "from_processed") {
        cur.esets <- Download_processed_from_GEO(cur.GSE.ID,
                                                 cache.folder=cache.folder)
        
      } else if (expt.annot$data.source == "from_raw") {
        # start by getting from_processed.  This is the substrate for the functions
        # that generate the from_raw data
        cur.preprocessed.esets <- GetEset_new(cur.GSE.ID, 
                                              eset.folder = eset.folder, 
                                              download.missing=download.missing,
                                              overwrite.existing=overwrite.existing,
                                              cache.folder = cache.folder,
                                              expt.annot=list(data.source = "from_processed"),
                                              annot.csv.folder=annot.csv.folder,
                                              feature.annotation.path=feature.annotation.path,
                                              removed.failed=F,    # since we only really need the data from the from_processed data, keep all esets even if they are empty
                                              verbose=verbose)
        
        # now get raw data and process it (will use pData from cur.preprocessed.esets)
        cur.esets <- lapply(cur.preprocessed.esets,
                            ProcessRawGEOData,
                            cache.folder, 
                            expt.annot=expt.annot,
                            verbose=verbose)
        # note that the relevant notes fields for each eset are added by  
        # ProcessRawGEOData() (for raw data processing).  
      }

 
      if (length(cur.esets)>0) {
        cur.esets <- lapply(cur.esets,
                            UpdateAnnotations,
                            annot.csv.folder=annot.csv.folder)
        
        cur.esets <- lapply(cur.esets,
                            map.features.to.EGID,
                            feature.annotation.path=feature.annotation.path,
                            cache.folder=cache.folder)
        
        # save the new esets to disk
        lapply(cur.esets, 
               SaveEset,
               eset.folder)      
      }
      
      # now that we've saved the data, swap out the failed from_raw with the from_processed.  note that
      # if the from processed failed as well it will be dealt with below (based on the remove.failed option)
      if (expt.annot$data.source == "from_raw") {
        from.raw.failed <- sapply(esets,
                         function(cur.eset) {
                           notes(cur.eset)$status != "ok"
                         })
        cur.esets[from.raw.failed] <- cur.preprocessed.esets[from.raw.failed]
      }
      
    }

    esets <- c(esets,
               cur.esets)    

  }
  
  if (remove.failed) {
    esets <- esets[sapply(esets,
                          function(cur.eset) {
                            notes(cur.eset)$status!="ok"
                          })]
  }
  
  return(invisible(esets))
}




ProcessRawGEOData_new <- function(cur.eset, cache.folder, expt.annot, verbose=T) {
  
  ##should use getesetname here and/or source load_data_funcs?
  cur.eset.name <- notes(cur.eset)$name
  # get the pData from the eset, so we can keep use the same annotations 
  # for the raw data
  cur.pData <- pData(cur.eset)
  
  # get the experiment notes from the eset
  expt.data <- notes(cur.eset)
  
  platform <- expt.data$platform
  
  # set up an empty eset to return if an error is generated
  eset.to.return = new("ExpressionSet", phenoData=cur.pData,
                   annotation=anotation(cur.eset))
  notes(eset.to.return) <- expt.data
  # replace the expt.annot field to just contain from_raw, so if error is genereated
  # before the actual data processing function then the failed file won't depend on
  # settings in expt.annot that don't apply (e.g., agilent data labeled with brainarray=T)
  notes(eset.to.return)$expt.annot=list(data.source = "from_raw")
  
  if (!("supplementary_file" %in% colnames(cur.pData))) {
    if (verbose) {cat("Raw data not available for '", cur.eset.name, "'. ",
                      "Using processed data for this data set.\n", sep="")}
    warning("Raw data not available for '", cur.eset.name, "'. ",
            "Using processed data for this data set.")
    
    notes(eset.to.return)$status <- "raw data not available"
    return(eset.to.return)
  }
  # get the names of the supp files to see if they are affy cel files
  supp.files <- as.character(cur.pData$supplementary_file)
  
  missing.sup.files <- toupper(supp.files) == "NONE"
  if (all(missing.sup.files)){
    if (verbose) {cat("Raw data not available for '", cur.eset.name, "'. ",
                      "Using processed data for this data set.\n", sep="")}
    warning("Raw data not available for '", cur.eset.name, "'. ",
            "Using processed data for this data set.")
    notes(eset.to.return)$status <- "raw data not available"
    return(eset.to.return)
    
  }
  if (any(missing.sup.files)){
    warning("Some raw data is missing from GEO for ", cur.eset.name, ".
            Specifically, files from ", toString(cur.pData$geo_accession[supp.files == "NONE"]),
            " are at missing. These files will be skipped when processing the data.")
    
    cur.pData <- cur.pData[!missing.sup.files, ,drop=F]
    supp.files <- supp.files[!missing.sup.files]
  }
  
  cur.ds.processing.func.name <- paste0("ProcessData_", cur.eset.name)
  
  # if there is a processing function for this specific dataset, then use that.
  # otherwise look for a processing function for the current platform
  if (exists(cur.ds.processing.func.name)) {
    cur.processing.func.name <- cur.ds.processing.func.name
  } else {
    
    # get the name of the function to process that platform
    cur.platform.processing.func.name  <- paste0("ProcessData_", platform)
    
    # check to see if it's an affy platform
    is.affy.file <- all(grepl("\\.cel$|\\.cel\\.|^none$", supp.files, ignore.case = TRUE))
    
    # if there isn't a specific function for that platform, but it is an affy platform, then
    # use the generic affy processing function
    if (!exists(cur.platform.processing.func.name) & is.affy.file) {
      cur.platform.processing.func.name <- "ProcessData_affy" 
    } 
    
    if (!exists(cur.platform.processing.func.name)) {
      if (verbose) {cat("Raw data processing function not available for '", cur.eset.name, "'. ",
                        "Using processed data for this data set.\n", sep="")}
      warning("Raw data processing function not available for '", 
              cur.eset.name, "'. ",
              "Using processed data for this data set.")
      notes(eset.to.return)$status <- "function to process raw data not available"
      return(eset.to.return)
    } else {
      cur.processing.func.name <- cur.platform.processing.func.name
    }
  }
  
  # this is where getGEOSuppFiles should save the files.  So check to see if
  # they exist first to avoid re-downloading something that we already have
  cur.gsm.files <- file.path(cache.folder, 
                             cur.pData$geo_accession, 
                             basename(supp.files))
  
  # for each file, see if we have it downloaded already.  If not, download it.
  downloaded.files <- 
    sapply(1:nrow(cur.pData),
           function(i) {
             if (file.exists(cur.gsm.files[i]) && file.info(cur.gsm.files[i])$size!=0) {
               return(cur.gsm.files[i])
             } else {
               return(rownames(getGEOSuppFiles_retry(cur.pData$geo_accession[i], 
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
    downloaded.gsm.files <- sapply(as.character(cur.pData$geo_accession),
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
      if (verbose) {cat("For ", cur.eset.name,
                        ", downloaded raw data files don't match with files names from ",
                        " annotation file, and can't be matched uniquely from GSM IDs.  ",
                        "Using processed data.")}
      warning("For ", cur.eset.name,
              ", downloaded raw data files don't match with files names from ",
              " annotation file, and can't be matched uniquely from GSM IDs.  ",
              "Using processed data.")
      notes(eset.to.return)$status <- "downloaded raw data files don't match with files names from annotation file"
      return(eset.to.return)
      
    }
    cur.gsm.files <- downloaded.gsm.files
  }
  
  # just in case, run normalizePath() to clean up the file paths
  cur.gsm.files <- normalizePath(cur.gsm.files)
  #run original processing function
  tmp.eset <- tryCatch({ 
    eval(parse(text=paste0(cur.processing.func.name,
                           "(cur.gsm.files, cur.eset, expt.annot=expt.annot)")))
  }, error = function(err){
    print(paste("Error: ", err))
    if (verbose) {cat("A error occurred while processing '", cur.eset.name, "'. ",
                      "Using processed data for this data set.\n", sep="")
    }
    return(NULL)
  })
  
  if (is.null(tmp.eset)) {
    notes(eset.to.return)$status <- "error generated during processing of raw data"
    return(eset.to.return)
  }
  
  
  # TODO: WHEN SUPP FILES ARE MISSING, ESET IS GENERATED WITHOUT THEM.  THIS TRIGGERS
  # THE FIRST WARNING BELOW AND PROCESSED DATA IS USED INSTEAD - UPDATE THE PDATA???
  
  # make sure that the number of rows in pData(tmp.eset) and 
  # notes(tmp.eset)$original.pData are the same as the number of columns/samples 
  # in tmp.eset
  if (nrow(pData(tmp.eset)) != nrow(notes(tmp.eset)$original.pData)) {
    warning(paste0(cur.processing.func.name, "() did not maintain the same number",
                   " of rows in pData and original.pData. Using processed data", 
                   " for ", cur.eset.name))
    notes(eset.to.return)$status <- "error generated during processing of raw data"
    return(eset.to.return)
  } else if(nrow(pData(tmp.eset)) != ncol(tmp.eset)) {
    warning(paste0(cur.processing.func.name, "() did not maintain the same number",
                   " of samples in pData AND the eset.  Using processed data", 
                   " for ", cur.eset.name))
    notes(eset.to.return)$status <- "error generated during processing of raw data"
    return(eset.to.return)
  } else {
    eset.to.return <- tmp.eset
  }
  
  # identify the items from expt.data that aren't in the eset, and add them
  expt.data.to.add <- expt.data[!(names(expt.data) %in% 
                                    names(notes(eset.to.return)))]
  notes(eset.to.return) <- c(notes(eset.to.return),
                       expt.data.to.add)
  # overwrite the data.source
  notes(eset.to.return)$data.source <- "from_raw"
  
  notes(eset.to.return)$status <- "ok"
  
  # set the platform here because the processing function might not know/capture
  # the platform, and we don't want it to change from the platform listed for
  # the pre-processed data
  annotation(cur.eset) <- notes(cur.eset)$platform
  
  return(cur.eset)
  
}


# MODIFY SAVEESET SO IT WILL REPLACE expt.annot=list(data.source="from_raw") with expt.annot=list(...) (ie, delete more
# general version if a more specific version appears)