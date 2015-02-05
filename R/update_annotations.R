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


UpdateAnnotations_GSE42296_GPL6244 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  # print(get.interesting.annot.cols(annot))
  
  # we can predict 3 things from this dataset:
  #   - response to infliximab in RA patients
  #   - response to infliximab in CD patients
  #   - RA vs CD 
  
  #get the baseline samples
  stopifnot(length(unique(annot$characteristics_ch1.4))==2)
  baseline <- annot$characteristics_ch1.4 == "time: week 0"

  #get the CD patients
  stopifnot(length(unique(annot$characteristics_ch1.5))==2)
  Crohns <- grepl("Crohn", annot$characteristics_ch1.5)
  
  # get the responders
  stopifnot(length(unique(annot$characteristics_ch1.2))==2)
  responder <- annot$characteristics_ch1.2 == "response: R - responder"
  
  # R/NR for CD patients
  new.annot$PRED_Crohns_inflix_response <- ifelse(baseline & Crohns, 
                                          ifelse(responder,
                                                 "R",
                                                 "NR"),
                                          NA)

  # R/NR for RA patients
  new.annot$PRED_RA_inflix_response <- ifelse(baseline & !Crohns, 
                                ifelse(responder,
                                       "R",
                                       "NR"),
                                NA)
  
  # RA vs CD
  new.annot$PRED_disease <- ifelse(!baseline,
                               NA,
                               ifelse(Crohns, 
                                      "CD",
                                      "RA"))

  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}


UpdateAnnotations_GSE17755_GPL1291 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  # get.interesting.annot.cols(annot)
  # table(annot[, c("submission_date", "last_update_date")])  # these columns are equivalent
  
  # we can predict 4 things from this dataset:
  #   - batch (submission_date)
  #   - gender (characteristics_ch1)
  #   - age (characteristics_ch1.1)
  #   - disease (characteristics_ch1.2)
  
  
  # batch date
  new.annot$PRED_submission_date <- annot$submission_date
  
  # gender
  stopifnot(length(unique(annot$characteristics_ch1))==2)
  new.annot$PRED_gender <- sub("gender: ", "", annot$characteristics_ch1)
  
  # age
  new.annot$PRED_age <- as.numeric(sub("age: ", "", annot$characteristics_ch1.1))
  stopifnot(identical(paste0("age: ", new.annot$PRED_age), as.character(annot$characteristics_ch1.1)))
  
  # disease
  new.annot$PRED_age <- sub("disease: ", "", annot$characteristics_ch1.2)

  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}
