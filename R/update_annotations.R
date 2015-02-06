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
# mulitple times sequentially with and produce identical results (which may not
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

UpdateAnnotations_GSE16879_GPL570 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  
  # design column
  new.annot$DESIGN_source_name_ch1 <- annot$source_name_ch1
  
  #contrast between normal and UCR before
  new.annot$CONTRAST_UCRbefore_Normal <- rep(0, times=nrow(new.annot))
  new.annot["Colonic mucosal biopsy from UC responder before first infliximab treatment" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCRbefore_Normal"] <- 1
  new.annot["Colonic mucosal biopsy from control individual" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCRbefore_Normal"] <- -1
  
  #contrast between normal and UCNR before
  new.annot$CONTRAST_UCNRbefore_Normal <- rep(0, times=nrow(new.annot))
  new.annot["Colonic mucosal biopsy from UC non-responder before first infliximab treatment" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCNRbefore_Normal"] <- 1
  new.annot["Colonic mucosal biopsy from control individual" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCNRbefore_Normal"] <- -1
  
  #contrast between normal and CDcR before
  new.annot$CONTRAST_CDcRbefore_Normal <- rep(0, times=nrow(new.annot))
  new.annot["Colonic mucosal biopsy from CDc responder before first infliximab treatment" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDcRbefore_Normal"] <- 1
  new.annot["Colonic mucosal biopsy from control individual" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDcRbefore_Normal"] <- -1
  
  #contrast between normal and CDcNR before
  new.annot$CONTRAST_CDcNRbefore_Normal <- rep(0, times=nrow(new.annot))
  new.annot["Colonic mucosal biopsy from CDc non-responder before first infliximab treatment" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDcNRbefore_Normal"] <- 1
  new.annot["Colonic mucosal biopsy from control individual" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDcNRbefore_Normal"] <- -1
  
  #contrast between ileal normal and CDiR before
  new.annot$CONTRAST_CDiRbefore_ilNormal <- rep(0, times=nrow(new.annot))
  new.annot["Ileal mucosal biopsy from CDi responder before first infliximab treatment" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDiRbefore_ilNormal"] <- 1
  new.annot["Ileal mucosal biopsy from control individual" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDiRbefore_ilNormal"] <- -1
  
  #contrast between ileal normal and CDiNR before
  new.annot$CONTRAST_CDiNRbefore_ilNormal <- rep(0, times=nrow(new.annot))
  new.annot["Ileal mucosal biopsy from CDi non-responder before first infliximab treatment" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDiNRbefore_ilNormal"] <- 1
  new.annot["Ileal mucosal biopsy from control individual" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDiNRbefore_ilNormal"] <- -1
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}

UpdateAnnotations_GSE6731_GPL8300 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  
  # design column
  new.annot$DESIGN_source_name_ch1 <- annot$source_name_ch1
  
  #contrast between normal and CD affected
  new.annot$CONTRAST_CDaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon with Crohn's disease" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDaff_Normal"] <- -1
  
  #contrast between normal and CD unaffected
  new.annot$CONTRAST_CDunaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with unaffected colon with Crohn's disease"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDunaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_CDunaff_Normal"] <- -1
  
  #contrast between normal and UC affected
  new.annot$CONTRAST_UCaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon from adult with UC"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCaff_Normal"] <- -1
  
  #contrast between normal and UC affected
  new.annot$CONTRAST_UCunaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with unaffected colon from adult with UC" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCunaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_UCunaff_Normal"] <- -1
  
  #contrast between normal and INF
  new.annot$CONTRAST_INF_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon from adult with INF" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_INF_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_INF_Normal"] <- -1
  
  #contrast between normal and IC
  new.annot$CONTRAST_IC_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon from adult with IC" == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_IC_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_IC_Normal"] <- -1
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}

UpdateAnnotations_GSE42568_GPL570 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  
  # design column
  new.annot$DESIGN_source_name_ch1 <- sapply(rownames(annot), function(annotrow) 
    paste0(annot[annotrow, "source_name_ch1"], " & ", annot[annotrow, "characteristics_ch1.2"]))
  #new.annot$characteristics_ch1.2 <- annot$characteristics_ch1.2
  
  #contrast between normal and breast cancer
#   new.annot$CONTRAST_BC_Normal <- rep(0, times=nrow(new.annot))
#   new.annot["Breast tissue, cancer" == new.annot$DESIGN_source_name_ch1,
#             "CONTRAST_BC_Normal"] <- 1
#   new.annot["Breast tissue, normal" == new.annot$DESIGN_source_name_ch1,
#             "CONTRAST_BC_Normal"] <- -1
  
  #contrast between normal and ER positive
  new.annot$CONTRAST_ERpos_Normal <- rep(0, times=nrow(new.annot))
  new.annot[("Breast tissue, cancer & er_status: 1"  == new.annot$DESIGN_source_name_ch1), "CONTRAST_ERpos_Normal"] <- 1
  new.annot["Breast tissue, normal & er_status: NA"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_ERpos_Normal"] <- -1
  
  #contrast between normal and ER negative
  new.annot$CONTRAST_ERneg_Normal <- rep(0, times=nrow(new.annot))
  new.annot[("Breast tissue, cancer & er_status: 0"  == new.annot$DESIGN_source_name_ch1), "CONTRAST_ERneg_Normal"] <- 1
  new.annot["Breast tissue, normal & er_status: NA"  == new.annot$DESIGN_source_name_ch1,
            "CONTRAST_ERneg_Normal"] <- -1
  
  #contrast between ER positive and ER negative
  new.annot$CONTRAST_ERPos_ERNeg <- rep(0, times=nrow(new.annot))
  new.annot[("Breast tissue, cancer & er_status: 1"  == new.annot$DESIGN_source_name_ch1), "CONTRAST_ERPos_ERNeg"] <- 1
  new.annot[("Breast tissue, cancer & er_status: 0"  == new.annot$DESIGN_source_name_ch1), "CONTRAST_ERPos_ERNeg"] <- -1

  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}