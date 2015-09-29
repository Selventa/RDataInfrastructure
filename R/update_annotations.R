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

# UpdateAnnotations_GSE16879_GPL570 <- function(cur.eset) {
#   
#   annot <- notes(cur.eset)$original.pData
#   new.annot <- data.frame(row.names=rownames(annot))
#   
#   
#   # design column
#   new.annot$DESIGN_source_name_ch <- annot$source_name_ch1
#   
#   UCR <- (new.annot$DESIGN_source_name_ch == 
#             "Colonic mucosal biopsy from UC responder before first infliximab treatment")
#   UCNR <- (new.annot$DESIGN_source_name_ch == 
#              "Colonic mucosal biopsy from UC non-responder before first infliximab treatment")
#   colon_control <- (new.annot$DESIGN_source_name_ch == 
#                       "Colonic mucosal biopsy from control individual")
#   
#   CDcR <- (new.annot$DESIGN_source_name_ch == 
#              "Colonic mucosal biopsy from CDc responder before first infliximab treatment")
#   CDcNR <- (new.annot$DESIGN_source_name_ch == 
#               "Colonic mucosal biopsy from CDc non-responder before first infliximab treatment")
#   
#   CDiR <- (new.annot$DESIGN_source_name_ch == 
#              "Ileal mucosal biopsy from CDi responder before first infliximab treatment")
#   CDiNR <- (new.annot$DESIGN_source_name_ch == 
#               "Ileal mucosal biopsy from CDi non-responder before first infliximab treatment")
#   
#   ilium_control <- (new.annot$DESIGN_source_name_ch == 
#                       "Ileal mucosal biopsy from control individual")
#   
#   
#   lapply(list(UCR, UCNR, colon_control, CDcR, CDcNR, CDiR, CDiNR, ilium_control),
#          function(x) stopifnot(any(x)))
#   
#   #contrast between normal and UCR before
#   new.annot$CONTRAST_UCRbefore_Normal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_UCRbefore_Normal[UCR] <- 1
#   new.annot$CONTRAST_UCRbefore_Normal[colon_control] <- -1
#   
#   
#   #contrast between normal and UCNR before
#   new.annot$CONTRAST_UCNRbefore_Normal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_UCNRbefore_Normal[UCNR] <- 1
#   new.annot$CONTRAST_UCNRbefore_Normal[colon_control] <- -1
#   
#   #contrast between normal and CDcR before
#   new.annot$CONTRAST_CDcRbefore_Normal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_CDcRbefore_Normal[CDcR] <- 1
#   new.annot$CONTRAST_CDcRbefore_Normal[colon_control] <- -1
#   
#   #contrast between normal and CDcNR before
#   new.annot$CONTRAST_CDcNRbefore_Normal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_CDcNRbefore_Normal[CDcNR] <- 1
#   new.annot$CONTRAST_CDcNRbefore_Normal[colon_control] <- -1
#   
#   #contrast between ileal normal and CDiR before
#   new.annot$CONTRAST_CDiRbefore_ilNormal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_CDiRbefore_ilNormal[CDiR] <- 1
#   new.annot$CONTRAST_CDiRbefore_ilNormal[ilium_control] <- -1
#   
#   #contrast between ileal normal and CDiNR before
#   new.annot$CONTRAST_CDiNRbefore_ilNormal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_CDiNRbefore_ilNormal[CDiNR] <- 1
#   new.annot$CONTRAST_CDiNRbefore_ilNormal[ilium_control] <- -1
#   
# 
#   #contrast between normal and UC (R+NR) before
#   new.annot$CONTRAST_UCbefore_Normal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_UCbefore_Normal[UCR | UCNR] <- 1
#   new.annot$CONTRAST_UCbefore_Normal[colon_control] <- -1
#   
#   #contrast between normal and CDc (R+NR) before
#   new.annot$CONTRAST_CDcbefore_Normal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_CDcbefore_Normal[CDcR | CDcNR] <- 1
#   new.annot$CONTRAST_CDcbefore_Normal[colon_control] <- -1
#   
#   #contrast between normal and CDi (R+NR) before
#   new.annot$CONTRAST_CDibefore_Normal <- rep(0, times=nrow(new.annot))
#   new.annot$CONTRAST_CDibefore_Normal[CDiR | CDiNR] <- 1
#   new.annot$CONTRAST_CDibefore_Normal[ilium_control] <- -1
# 
#   pData(cur.eset) <- new.annot
#   
#   return(cur.eset)
#   
# }

UpdateAnnotations_GSE6731_GPL8300 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  
  # design column
  new.annot$DESIGN_source_name_ch <- annot$source_name_ch1
  
  #contrast between normal and CD affected
  new.annot$CONTRAST_CDaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon with Crohn's disease" == new.annot$DESIGN_source_name_ch,
            "CONTRAST_CDaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult" == new.annot$DESIGN_source_name_ch,
            "CONTRAST_CDaff_Normal"] <- -1
  
  #contrast between normal and CD unaffected
  new.annot$CONTRAST_CDunaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with unaffected colon with Crohn's disease"  == new.annot$DESIGN_source_name_ch,
            "CONTRAST_CDunaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch,
            "CONTRAST_CDunaff_Normal"] <- -1
  
  #contrast between normal and UC affected
  new.annot$CONTRAST_UCaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon from adult with UC"  == new.annot$DESIGN_source_name_ch,
            "CONTRAST_UCaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch,
            "CONTRAST_UCaff_Normal"] <- -1
  
  #contrast between normal and UC affected
  new.annot$CONTRAST_UCunaff_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with unaffected colon from adult with UC" == new.annot$DESIGN_source_name_ch,
            "CONTRAST_UCunaff_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch,
            "CONTRAST_UCunaff_Normal"] <- -1
  
  #contrast between normal and INF
  new.annot$CONTRAST_INF_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon from adult with INF" == new.annot$DESIGN_source_name_ch,
            "CONTRAST_INF_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch,
            "CONTRAST_INF_Normal"] <- -1
  
  #contrast between normal and IC
  new.annot$CONTRAST_IC_Normal <- rep(0, times=nrow(new.annot))
  new.annot["colonoscopic biopsy from adult with affected colon from adult with IC" == new.annot$DESIGN_source_name_ch,
            "CONTRAST_IC_Normal"] <- 1
  new.annot["colonoscopic biopsy from normal adult"  == new.annot$DESIGN_source_name_ch,
            "CONTRAST_IC_Normal"] <- -1
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}

UpdateAnnotations_GSE42568_GPL570 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  
  # design column
  new.annot$DESIGN_1source_name_ch1 <- sapply(rownames(annot), function(annotrow) 
    paste0(annot[annotrow, "source_name_ch1"], " & ", annot[annotrow, "characteristics_ch1.2"]))
  new.annot$DESIGN_2source_name_ch1 <- annot[, "source_name_ch1"]
  
  #contrast between normal and breast cancer
  new.annot$CONTRAST_2BC_Normal <- rep(0, times=nrow(new.annot))
  new.annot["Breast tissue, cancer" == new.annot$DESIGN_2source_name_ch1,
            "CONTRAST_2BC_Normal"] <- 1
  new.annot["Breast tissue, normal" == new.annot$DESIGN_2source_name_ch1,
            "CONTRAST_2BC_Normal"] <- -1
  
  #contrast between normal and ER positive
  new.annot$CONTRAST_1ERpos_Normal <- rep(0, times=nrow(new.annot))
  new.annot[("Breast tissue, cancer & er_status: 1"  == new.annot$DESIGN_1source_name_ch1), "CONTRAST_1ERpos_Normal"] <- 1
  new.annot["Breast tissue, normal & er_status: NA"  == new.annot$DESIGN_1source_name_ch1,
            "CONTRAST_1ERpos_Normal"] <- -1
  
  #contrast between normal and ER negative
  new.annot$CONTRAST_1ERneg_Normal <- rep(0, times=nrow(new.annot))
  new.annot[("Breast tissue, cancer & er_status: 0"  == new.annot$DESIGN_1source_name_ch1), "CONTRAST_1ERneg_Normal"] <- 1
  new.annot["Breast tissue, normal & er_status: NA"  == new.annot$DESIGN_1source_name_ch1,
            "CONTRAST_1ERneg_Normal"] <- -1
  
  #contrast between ER positive and ER negative
  new.annot$CONTRAST_1ERPos_ERNeg <- rep(0, times=nrow(new.annot))
  new.annot[("Breast tissue, cancer & er_status: 1"  == new.annot$DESIGN_1source_name_ch1), "CONTRAST_1ERPos_ERNeg"] <- 1
  new.annot[("Breast tissue, cancer & er_status: 0"  == new.annot$DESIGN_1source_name_ch1), "CONTRAST_1ERPos_ERNeg"] <- -1

  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}

UpdateAnnotations_GSE13911_GPL570 <- function(cur.eset) {
  
  annot <- notes(cur.eset)$original.pData
  new.annot <- data.frame(row.names=rownames(annot))
  
  #contrast between normal and tumor
  new.annot$CONTRAST_Tumor_Normal <- rep(0, times=nrow(new.annot))
  new.annot["TUMOR" == annot$characteristics_ch1.1,
            "CONTRAST_Tumor_Normal"] <- 1
  new.annot["NORMAL" == annot$characteristics_ch1.1,
            "CONTRAST_Tumor_Normal"] <- -1
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
  
}


UpdateAnnotations_GSE10893_GPL1390 <- function(cur.eset) {
  
  ## First, get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  
  ## NOTE! Three samples have some misalignments in the phenotype table: GSM275775, GSM372557, & GSM372558
  ## Fix these annotations manually.
  next.sample <- match("GSM275775", annot$geo_accession)
  annot$source_name_ch2[next.sample] <- annot$characteristics_ch2[next.sample]
  annot$characteristics_ch2[next.sample] <- "age: --"
  
  
  next.sample <- match("GSM372557", annot$geo_accession)
  annot$source_name_ch2[next.sample] <- paste(annot$characteristics_ch2[next.sample], annot$characteristics_ch2.1[next.sample], sep=": ")
  annot$characteristics_ch2[next.sample] <- ""
  annot$characteristics_ch2.1[next.sample] <- ""
  
  
  next.sample <- match("GSM372558", annot$geo_accession)
  to.cols   <- c("characteristics_ch2",   "characteristics_ch2.1", "characteristics_ch2.4", "characteristics_ch2.5", "characteristics_ch2.6", "characteristics_ch2.7", "characteristics_ch2.8", "characteristics_ch2.9", "characteristics_ch2.10", "characteristics_ch2.11")
  from.cols <- c("characteristics_ch2.1", "characteristics_ch2.2", "characteristics_ch2.3", "characteristics_ch2.4", "characteristics_ch2.5", "characteristics_ch2.6", "characteristics_ch2.7", "characteristics_ch2.8", "characteristics_ch2.9",  "characteristics_ch2.10")
  annot[next.sample, to.cols] <- unlist(annot[next.sample, from.cols])
  annot$characteristics_ch2.2[next.sample] <- ""
  annot$characteristics_ch2.3[next.sample] <- ""
  
  
  ## Clean up some of the annotations to make them more computable:
  new.annot$title <- annot$title
  new.annot$geo_accession <- annot$geo_accession
  new.annot$sample_source <- annot$source_name_ch2  ## PROCESS THIS TO ONLTY CONTAIN ONE OF "breast tumor", "normal breast","lymph node met", "autopsy", or "autopsy - skin met". 
  temp.idx <- grep("breast tumor", annot$source_name_ch2, ignore.case = TRUE)
  new.annot$sample_source[temp.idx] <- "Breast Tumor"
  
  temp.idx <- grep("metastasis.*brain", annot$source_name_ch2, ignore.case = TRUE)
  new.annot$sample_source[temp.idx] <- "Brain Metastasis"
  
  ## Manually reset this entry to be a breast tumor (it has "brain metastasis" in the name, but is a primary breast tumor).
  new.annot$sample_source[new.annot$geo_accession == "GSM34432"] <- "Breast Tumor"
  
  temp.idx <- grep("Normal.*Breast", annot$source_name_ch2, ignore.case = TRUE)
  new.annot$sample_source[temp.idx] <- "Normal Breast"
  
  ###########################################################################=
  ## FINISH!!!  This is a rough state, and doesn't indicate which are from autopsy, or are the appropriate  ----
  ##  combinations of Autopsy and metastasis.
  ###########################################################################=
  
  
  #######  FINISH!  #################
  ## Now get the rest of the fields that I want and process them: ER/PR/HER2 statuses, PAM50/Claudin subtype.
  ##  LATER: Get the following fields: Age, node status, grade, size, response (rfs & overall survival).
  ##  ALSO, figure out who got what treatment!...  And was it before or after the gene expression measurment.
  
  new.annot$ER_status <- "-NA-"
  new.annot$ER_status[annot$characteristics_ch2.1 == "er  (1=positive; 0=negative): 1"] <- "pos"
  new.annot$ER_status[annot$characteristics_ch2.1 == "er  (1=positive; 0=negative): 0"] <- "neg"
  
  new.annot$PR_status <- "-NA-"
  new.annot$PR_status[annot$characteristics_ch2.2 == "pgr  (1=positive; 0=negative): 1"] <- "pos"
  new.annot$PR_status[annot$characteristics_ch2.2 == "pgr  (1=positive; 0=negative): 0"] <- "neg"
  
  new.annot$HER2_clin_status <- "-NA-"
  new.annot$HER2_clin_status[annot$characteristics_ch2.3 == "her2 clinical status ihc/fish (1=positive; 0=negative): 1"] <- "pos"
  new.annot$HER2_clin_status[annot$characteristics_ch2.3 == "her2 clinical status ihc/fish (1=positive; 0=negative): 0"] <- "neg"
  
  ## Get the PAM50 intrinsic subtype or claudin-low classification of each tumor.
  new.annot$intrinsic_subtype <- gsub("pam50 predictions plus claudin-low classification (cell line predictor): ", 
                                      "", 
                                      annot$characteristics_ch2.11, 
                                      fixed = TRUE)
  
  
  ##--------------------------------------------------------=
  ## Create the contrasts:
  ##--------------------------------------------------------=
  
  ## Contrast for HER2 breast tumors (by IHC, not by PAM50).
  new.annot$CONTRAST_HER2pos_vs_normal <- 0
  new.annot$CONTRAST_HER2pos_vs_normal[ ((new.annot$HER2_clin_status == "pos") 
                                         & (new.annot$sample_source == "Breast Tumor")) ] <- 1
  new.annot$CONTRAST_HER2pos_vs_normal[ (new.annot$sample_source == "Normal Breast") ]    <- (-1)
  
  ## Contrast for hormone receptor positive (ER+ and/or PR+) vs normal.
  new.annot$CONTRAST_HormoneRecepPos_vs_normal <- 0
  new.annot$CONTRAST_HormoneRecepPos_vs_normal[ ( (new.annot$sample_source == "Breast Tumor")
                                                  & ( ((new.annot$ER_status == "pos") 
                                                       | (new.annot$PR_status == "pos"))
                                                      & (new.annot$HER2_clin_status == "neg")) ) ] <- 1
  new.annot$CONTRAST_HormoneRecepPos_vs_normal[ (new.annot$sample_source == "Normal Breast") ]     <- (-1)
  
  ## Contrast for LuminalA vs normal.
  new.annot$CONTRAST_LumA_vs_normal <- 0
  new.annot$CONTRAST_LumA_vs_normal[ ((new.annot$intrinsic_subtype == "LumA") 
                                      & (new.annot$sample_source == "Breast Tumor")) ] <- 1
  new.annot$CONTRAST_LumA_vs_normal[ (new.annot$sample_source == "Normal Breast") ]    <- (-1)
  
  ## Contrast for LuminalB vs normal.
  new.annot$CONTRAST_LumB_vs_normal <- 0
  new.annot$CONTRAST_LumB_vs_normal[ ((new.annot$intrinsic_subtype == "LumB") 
                                      & (new.annot$sample_source == "Breast Tumor")) ] <- 1
  new.annot$CONTRAST_LumB_vs_normal[ (new.annot$sample_source == "Normal Breast") ]    <- (-1)
  
  ## Contrast for (LuminalA and LuminalB) vs normal.  This should be roughly the 
  ##  same as the hormone-receptor-positives vs normal, but just to be thorough.
  new.annot$CONTRAST_LumA_and_LumB_vs_normal <- 0
  new.annot$CONTRAST_LumA_and_LumB_vs_normal[ ((new.annot$intrinsic_subtype %in% c("LumA", "LumB")) 
                                               & (new.annot$sample_source == "Breast Tumor")) ] <- 1
  new.annot$CONTRAST_LumA_and_LumB_vs_normal[ (new.annot$sample_source == "Normal Breast") ]    <- (-1)
  
  ## Contrast for All breast tumors vs normal.
  new.annot$CONTRAST_AllBreastTumors_vs_normal <- 0
  new.annot$CONTRAST_AllBreastTumors_vs_normal[ (new.annot$sample_source == "Breast Tumor") ]  <- 1
  new.annot$CONTRAST_AllBreastTumors_vs_normal[ (new.annot$sample_source == "Normal Breast") ] <- (-1)
  
  ## Contrast for Triple-negative breast tumors vs normal.
  new.annot$TripleNegative_vs_normal <- 0
  new.annot$TripleNegative_vs_normal[ ((new.annot$sample_source == "Breast Tumor") 
                                       & (new.annot$ER_status == "neg") 
                                       & (new.annot$PR_status == "neg") 
                                       & (new.annot$HER2_clin_status == "neg")) ]  <- 1
  new.annot$TripleNegative_vs_normal[ (new.annot$sample_source == "Normal Breast") ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}
