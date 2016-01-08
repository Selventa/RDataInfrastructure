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
  
  ##--------------------------------------------------------=
  ## Fix misaligned annotations:
  ##--------------------------------------------------------=
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
  
  ##---------=
  ## Below are three different possible contrasts to represent luminal
  ##  tumors, as compared to normal breast.  We're not sure yet which one 
  ##  will best represent the biological differences that characterize 
  ##  luminal breast tumors, so all three are included for possible
  ##  analysis.
  ##---------=
  
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
  
  ##---------=
  ## Below are three different possible contrasts to represent Triple-Negative / Basal
  ##  tumors, as compared to normal breast.  (Similar to the multiple contrasts
  ##  defined for luminal tumors above.)
  ##---------=
  ## Contrast for Triple-negative breast tumors vs normal.
  new.annot$CONTRAST_TripleNegative_vs_normal <- 0
  new.annot$CONTRAST_TripleNegative_vs_normal[ ((new.annot$sample_source == "Breast Tumor") 
                                                & (new.annot$ER_status == "neg") 
                                                & (new.annot$PR_status == "neg") 
                                                & (new.annot$HER2_clin_status == "neg")) ]  <- 1
  new.annot$CONTRAST_TripleNegative_vs_normal[ (new.annot$sample_source == "Normal Breast") ] <- (-1)
  
  ## Contrast for Basal-like breast tumors (as defined by PAM50 intrinsic-subtype) vs normal.
  new.annot$CONTRAST_Basal_vs_normal <- 0
  new.annot$CONTRAST_Basal_vs_normal[ ((new.annot$sample_source == "Breast Tumor") 
                                       & (new.annot$intrinsic_subtype == "Basal"))] <- 1
  new.annot$CONTRAST_Basal_vs_normal[ (new.annot$sample_source == "Normal Breast") ] <- (-1)
  
  ## Contrast for Basal-like breast tumors vs normal, where basal-like 
  ##  consists of either receptor-status-definition (ER-/PR-/HER2-) or
  ##  PAM50 intrinsic subtype.
  new.annot$CONTRAST_TNBC_or_Basal_vs_normal <- 0
  new.annot$CONTRAST_TNBC_or_Basal_vs_normal[ (new.annot$CONTRAST_TripleNegative_vs_normal == 1)
                                              | (new.annot$CONTRAST_Basal_vs_normal == 1) ] <- 1
  new.annot$CONTRAST_TNBC_or_Basal_vs_normal[ (new.annot$sample_source == "Normal Breast") ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE5364_GPL96 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Note - Data annotations appear to be well formatted.  No adjustments necessary!
  
  ##--------------------------------------------------------=
  ## Create the contrast:
  ## -----------=
  ##  - We do not have subtype information for the breast tumors, so we can only
  ##    compare All breast tumors vs breast normal tissue.
  ##  - Note that there are other tumor tissue samples here, but we are only
  ##    interested in Breast tumor-vs-normal for now.
  ##--------------------------------------------------------=
  ## Get the index for the breast tumors, and for the breast normals:
  breast.tumor.idx <- grep(pattern = "Breast tumor", annot$title, ignore.case = TRUE)
  breast.normal.idx <- grep(pattern = "Breast Normal", annot$title, ignore.case = TRUE)
  
  ## Contrast for All breast tumors vs normal.
  new.annot$CONTRAST_AllBreastTumors_vs_normal <- 0
  new.annot$CONTRAST_AllBreastTumors_vs_normal[ breast.tumor.idx ]  <- 1
  new.annot$CONTRAST_AllBreastTumors_vs_normal[ breast.normal.idx ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE14297_GPL6370 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Data annotations appear to be well formatted.  No adjustments necessary!
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_PrimaryCRCTumor_vs_Normal = Compare primary colorectal tumors 
  ##        versus normal colon epithelium.
  ##  - CONTRAST_CRCLiverMets_vs_NormalLiver = Comparison of the liver metastases 
  ##        (metastastases from primary CRC tumors) compared to normal liver tissue.
  ##
  ##-----------=
  ##  * NOTE: Since the liver metastases are paired with the primary tumors, 
  ##    there are more contrasts that we could investigate when not looking
  ##    for positive controls.
  ##---------------------------------------------------------------------------=
  ## Get the index for the different tumors and normals:
  primary.CRC.tumor.idx <- grep(pattern = "Primary Colorectal Cancer", annot$source_name_ch1, ignore.case = TRUE)
  normal.colon.idx <- grep(pattern = "Normal Colon Epithelium", annot$source_name_ch1, ignore.case = TRUE)
  liver.mets.idx <- grep(pattern = "Liver metastasis of Colorectal Cancer", annot$source_name_ch1, ignore.case = TRUE)
  normal.liver.idx <- grep(pattern = "Normal Liver Tissue", annot$source_name_ch1, ignore.case = TRUE)
  
  ## Contrast for Primary Colorectal tumor vs normal colon.
  new.annot$CONTRAST_PrimaryCRCTumor_vs_Normal <- 0
  new.annot$CONTRAST_PrimaryCRCTumor_vs_Normal[ primary.CRC.tumor.idx ]  <- 1
  new.annot$CONTRAST_PrimaryCRCTumor_vs_Normal[ normal.colon.idx ] <- (-1)
  
  ## Contrast for Liver metasteses vs normal liver.
  new.annot$CONTRAST_CRCLiverMets_vs_NormalLiver <- 0
  new.annot$CONTRAST_CRCLiverMets_vs_NormalLiver[ liver.mets.idx ]  <- 1
  new.annot$CONTRAST_CRCLiverMets_vs_NormalLiver[ normal.liver.idx ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE18105_GPL570 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  ## Also keep the titles for this dataset, since they might be useful for reminding us
  ##  which samples are paired, and which are Laser-Capture Microdisection (LCM).
  new.annot$title <- annot$title
  
  ## Data annotations appear to be well formatted.  No adjustments necessary!
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_PairedCRCTumor_vs_Normal = Compare only the 17 colorectal tumors 
  ##        that have adjacent non-cancerous colon samples, to those non-cancerous
  ##        "normals".
  ##  - CONTRAST_AllCRCTumors_vs_Normal = Comparison of the 77 LCM tumor samples + 17 paired tumors
  ##        to the 17 "normal" adjacents.
  ##
  ##-----------=
  ##  * NOTE: These are effectively two different ways of looking at the same 
  ##    contrast (i.e. CRC vs normal).  We need to take this into account 
  ##    when considering results.
  ##---------------------------------------------------------------------------=
  ## Get the index for the different tumors and normals:
  paired.CRC.tumor.idx <- grep(pattern = "cancer, homogenized", annot$title, ignore.case = TRUE)
  paired.normal.colon.idx <- grep(pattern = "normal, homogenized", annot$title, ignore.case = TRUE)
  LCM.CRC.tumor.idx <- grep(pattern = "cancer, LCM", annot$title, ignore.case = TRUE)
  
  ## Contrast for Primary Colorectal tumor vs normal colon.
  new.annot$CONTRAST_PairedCRCTumor_vs_Normal <- 0
  new.annot$CONTRAST_PairedCRCTumor_vs_Normal[ paired.CRC.tumor.idx ]  <- 1
  new.annot$CONTRAST_PairedCRCTumor_vs_Normal[ paired.normal.colon.idx ] <- (-1)
  
  ## Contrast for Liver metasteses vs normal liver.
  new.annot$CONTRAST_AllCRCTumors_vs_Normal <- 0
  new.annot$CONTRAST_AllCRCTumors_vs_Normal[ c(LCM.CRC.tumor.idx, paired.CRC.tumor.idx) ]  <- 1
  new.annot$CONTRAST_AllCRCTumors_vs_Normal[ paired.normal.colon.idx ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}



UpdateAnnotations_GSE5206_GPL570 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  ## Also keep the titles for this dataset, since they might be useful for reminding us
  ##  which samples are paired, and which are Laser-Capture Microdisection (LCM).
  new.annot$title <- annot$title
  
  ## Data annotations appear to be well formatted.  No adjustments necessary!
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_CRCTumor_vs_Normal = the 100 CRC tumors vs 5 normals.
  ##---------------------------------------------------------------------------=
  ## Get the index for the different tumors and normals:
  tumor.idx <- grep(pattern = "ID: T", as.character(annot$title))
  normal.idx <- grep(pattern = "ID: N", as.character(annot$title))
  
  ## Contrast for Colorectal tumor vs normal colon.
  new.annot$CONTRAST_CRCTumor_vs_Normal <- 0
  new.annot$CONTRAST_CRCTumor_vs_Normal[ tumor.idx ]  <- 1
  new.annot$CONTRAST_CRCTumor_vs_Normal[ normal.idx ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE20881_GPL1708 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  ## Also keep the titles for this dataset, since they might be useful for reminding us
  ##  which samples are paired, and which are Laser-Capture Microdisection (LCM).
  new.annot$sample_source <- annot$source_name_ch1
  
  ## Data annotations appear to be well formatted.  No adjustments necessary!
  
  ##---------------------------------------------------------------------------=
  ## Create the contrasts:
  ##-----------=
  ##  - CONTRAST_AscendColon_Crohns_vs_Healthy
  ##  - CONTRAST_DescendColon_Crohns_vs_Healthy
  ##  - CONTRAST_SigmoidColon_Crohns_vs_Healthy
  ##  - CONTRAST_TerminalIleumColon_Crohns_vs_Healthy
  ##---------------------------------------------------------------------------=
  
  ## Contrast for Ascending colon Crohn's vs normal.
  crohns.idx <- grep("ascending colon biopsy from crohns disease subject", as.character(annot$source_name_ch1))
  normal.idx <- grep("ascending colon biopsy from healthy subject", as.character(annot$source_name_ch1))
  new.annot$CONTRAST_AscendColon_Crohns_vs_Healthy <- 0
  new.annot$CONTRAST_AscendColon_Crohns_vs_Healthy[ crohns.idx ]  <- 1
  new.annot$CONTRAST_AscendColon_Crohns_vs_Healthy[ normal.idx ] <- (-1)
  
  ## Contrast for descending colon Crohn's vs normal.
  crohns.idx <- grep("descending colon biopsy from crohns disease subject", as.character(annot$source_name_ch1))
  normal.idx <- grep("descending colon biopsy from healthy subject", as.character(annot$source_name_ch1))
  new.annot$CONTRAST_DescendColon_Crohns_vs_Healthy <- 0
  new.annot$CONTRAST_DescendColon_Crohns_vs_Healthy[ crohns.idx ]  <- 1
  new.annot$CONTRAST_DescendColon_Crohns_vs_Healthy[ normal.idx ] <- (-1)
  
  ## Contrast for sigmoid colon Crohn's vs normal.
  crohns.idx <- grep("sigmoid colon biopsy from crohns disease subject", as.character(annot$source_name_ch1))
  normal.idx <- grep("sigmoid colon biopsy from healthy subject", as.character(annot$source_name_ch1))
  new.annot$CONTRAST_SigmoidColon_Crohns_vs_Healthy <- 0
  new.annot$CONTRAST_SigmoidColon_Crohns_vs_Healthy[ crohns.idx ]  <- 1
  new.annot$CONTRAST_SigmoidColon_Crohns_vs_Healthy[ normal.idx ] <- (-1)
  
  ## Contrast for terminal ileum Crohn's vs normal.
  crohns.idx <- grep("terminal ileum biopsy from crohns disease subject", as.character(annot$source_name_ch1))
  normal.idx <- grep("terminal ileum biopsy from healthy subject", as.character(annot$source_name_ch1))
  new.annot$CONTRAST_TerminalIleumColon_Crohns_vs_Healthy <- 0
  new.annot$CONTRAST_TerminalIleumColon_Crohns_vs_Healthy[ crohns.idx ]  <- 1
  new.annot$CONTRAST_TerminalIleumColon_Crohns_vs_Healthy[ normal.idx ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE24287_GPL6480 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the source of the tissue as another field.
  new.annot$disease_phenotype <- gsub(pattern = "disease phenotype: ", "", annot$characteristics_ch1.1)
  
  
  ##---------------------------------------------------------------------------=
  ## Create the contrasts:
  ##-----------=
  ##  - CONTRAST_unaffIlealCrohns_vs_NoIBD 
  ##  - CONTRAST_unaffUC_vs_NoIBD
  ##-----------=
  ##  * SIDE NOTE: These are 2-chanel chips an the second channel is a comon
  ##    reference RNA sample from a healthy individual's (no-IBD) ileum.
  ##---------------------------------------------------------------------------=
  ## Get the index for the different tumors and normals:
  ileal.Crohns.idx <- grep(pattern = "illeal Crohn's disease", new.annot$disease_phenotype, ignore.case = TRUE)
  UC.idx <- grep(pattern = "ulcerative colitis", new.annot$disease_phenotype, ignore.case = TRUE)
  no.IBD.idx <- grep(pattern = "no inflammatory bowel disease", new.annot$disease_phenotype, ignore.case = TRUE)
  
  new.annot$CONTRAST_unaffIlealCrohns_vs_NoIBD <- 0
  new.annot$CONTRAST_unaffIlealCrohns_vs_NoIBD[ ileal.Crohns.idx ]  <- 1
  new.annot$CONTRAST_unaffIlealCrohns_vs_NoIBD[ no.IBD.idx ] <- (-1)
  
  new.annot$CONTRAST_unaffUC_vs_NoIBD <- 0
  new.annot$CONTRAST_unaffUC_vs_NoIBD[ UC.idx ]  <- 1
  new.annot$CONTRAST_unaffUC_vs_NoIBD[ no.IBD.idx ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}




UpdateAnnotations_GSE14323_GPL571 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the source of the tissue as another field.
  new.annot$sample_source <- annot$source_name_ch1
  
  
  ##---------------------------------------------------------------------------=
  ## Create the contrasts:
  ##-----------=
  ##  - CONTRAST_HCV.cirrhosis.no.HCC_vs_NormalLiver 
  ##  - CONTRAST_HCC_vs_NormalLiver
  ##  - CONTRAST_no.HCC.cirrhosis_vs_HCC.cirrhosis
  ##  - CONTRAST_HCC_vs_adjacent.cirrhotic
  ##  - CONTRAST_all.cirrhotic_vs_NormalLiver
  ##-----------=
  ##  * SIDE NOTE: Probably won't use all of these.  We'll probably only add
  ##    one or two of the contrasts to the targets list.  Though this is worth
  ##    revisiting to determine which diseases to investigate.
  ##---------------------------------------------------------------------------=
  ## Get the index for the different tumors and normals:
  cirrhosis.no.HCC.idx <- which(new.annot$sample_source == "liver tissue with cirrhosis")
  cirrhosis.yes.HCC.idx <- which(new.annot$sample_source == "liver tissue with cirrhosisHCC")
  HCC.idx <- which(new.annot$sample_source == "liver tissue with HCC")
  healthy.liver.idx <- which(new.annot$sample_source == "Normal liver tissue")
  
  new.annot$CONTRAST_HCV.cirrhosis.no.HCC_vs_NormalLiver <- 0
  new.annot$CONTRAST_HCV.cirrhosis.no.HCC_vs_NormalLiver[ cirrhosis.no.HCC.idx ]  <- 1
  new.annot$CONTRAST_HCV.cirrhosis.no.HCC_vs_NormalLiver[ healthy.liver.idx ] <- (-1)
  
  new.annot$CONTRAST_HCC_vs_NormalLiver <- 0
  new.annot$CONTRAST_HCC_vs_NormalLiver[ HCC.idx ]  <- 1
  new.annot$CONTRAST_HCC_vs_NormalLiver[ healthy.liver.idx ] <- (-1)
  
  new.annot$CONTRAST_no.HCC.cirrhosis_vs_HCC.cirrhosis <- 0
  new.annot$CONTRAST_no.HCC.cirrhosis_vs_HCC.cirrhosis[ cirrhosis.yes.HCC.idx ]  <- 1
  new.annot$CONTRAST_no.HCC.cirrhosis_vs_HCC.cirrhosis[ cirrhosis.no.HCC.idx ] <- (-1)
  
  new.annot$CONTRAST_HCC_vs_adjacent.cirrhotic <- 0
  new.annot$CONTRAST_HCC_vs_adjacent.cirrhotic[ HCC.idx ]  <- 1
  new.annot$CONTRAST_HCC_vs_adjacent.cirrhotic[ cirrhosis.yes.HCC.idx ] <- (-1)
  
  new.annot$CONTRAST_all.cirrhotic_vs_NormalLiver <- 0
  new.annot$CONTRAST_all.cirrhotic_vs_NormalLiver[ c(cirrhosis.no.HCC.idx, cirrhosis.yes.HCC.idx) ]  <- 1
  new.annot$CONTRAST_all.cirrhotic_vs_NormalLiver[ healthy.liver.idx  ] <- (-1)
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}



UpdateAnnotations_GSE6764_GPL570 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the tumor descriptors.
  new.annot$sample_description <- annot$characteristics_ch1
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ## For now, let's try grouping all of the HCC's together, and the cirrhosis 
  ## samples together, and compare each to normal liver.
  ##
  ##  - CONTRAST_HCC_vs_NormalLiver.
  ##  - CONTRAST_Cirrhosis_vs_NormalLiver.
  ##-----------=
  ##  * NOTES: This grouping combines samples that perhaps should be considered
  ##    independently.  For the HCC group, all of the different HCC timepoints
  ##    are grouped together, and for cirrhosis, I'm combining both the cirrhotic
  ##    tissues from individuals that have HCC, and those that don't.
  ##     ALSO, for now I'm ignoring the displastic samples, though I will check
  ##    with David whether those should be included in one of the other groups.
  ##---------------------------------------------------------------------------=
  ## Get the index for the different diseases' samples and normals:
  cirrhotic.idx <- c(which(new.annot$sample_description == "cirrhotic liver tissue"),
                     which(new.annot$sample_description == "cirrhotic liver tissue from patients without HCC"))
  HCC.idx <- c(which(new.annot$sample_description == "very early HCC"),
               which(new.annot$sample_description == "early HCC"),
               which(new.annot$sample_description == "advanced HCC"),
               which(new.annot$sample_description == "very advanced HCC"))
  normal.idx <- which(new.annot$sample_description == "normal liver tissue")
  
  ## Contrast for Primary Colorectal tumor vs normal colon.
  new.annot$CONTRAST_HCC_vs_NormalLiver <- 0
  new.annot$CONTRAST_HCC_vs_NormalLiver[ HCC.idx ]  <- 1
  new.annot$CONTRAST_HCC_vs_NormalLiver[ normal.idx ] <- (-1)
  
  ## Contrast for Liver metasteses vs normal liver.
  new.annot$CONTRAST_Cirrhosis_vs_NormalLiver <- 0
  new.annot$CONTRAST_Cirrhosis_vs_NormalLiver[ cirrhotic.idx ]  <- 1
  new.annot$CONTRAST_Cirrhosis_vs_NormalLiver[ normal.idx ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}



UpdateAnnotations_GSE3189_GPL96 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the sample descriptors.
  new.annot$sample_description <- annot$characteristics_ch1
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_Melanoma_vs_NormalSkin.
  ##---------------------------------------------------------------------------=
  ## Get the index for the different diseases' samples and normals:
  melanoma.idx <- which(new.annot$sample_description == "Melanoma")
  normal.idx <- which(new.annot$sample_description == "Normal")
  
  new.annot$CONTRAST_Melanoma_vs_NormalSkin <- 0
  new.annot$CONTRAST_Melanoma_vs_NormalSkin[ melanoma.idx ]  <- 1
  new.annot$CONTRAST_Melanoma_vs_NormalSkin[ normal.idx ] <- (-1)
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}

UpdateAnnotations_GSE19804_GPL570 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the sample descriptors.
  new.annot$sample_description <- annot$characteristics_ch1
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_LungCancer_vs_NormalAdjacent
  ##---------------------------------------------------------------------------=
  ## Get the index for the different diseases' samples and normals:
  lung.cancer.idx <- which(new.annot$sample_description == "tissue: lung cancer")
  normal.idx <- which(new.annot$sample_description == "tissue: paired normal adjacent")
  
  new.annot$CONTRAST_LungCancer_vs_NormalAdjacent <- 0
  new.annot$CONTRAST_LungCancer_vs_NormalAdjacent[ lung.cancer.idx ]  <- 1
  new.annot$CONTRAST_LungCancer_vs_NormalAdjacent[ normal.idx ] <- (-1)
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE32676_GPL570 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the sample description
  new.annot$sample_description <- gsub(pattern = "disease status: ", "", annot$characteristics_ch1.1)
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_PancreaticCancer_vs_NormalPancreas
  ##---------------------------------------------------------------------------=
  ## Get the index for the different diseases' samples and normals:
  tumor.idx <- which(new.annot$sample_description == "pancreatic cancer")
  normal.idx <- which(new.annot$sample_description == "non-malignant pancreas")
  
  new.annot$CONTRAST_PancreaticCancer_vs_NormalPancreas <- 0
  new.annot$CONTRAST_PancreaticCancer_vs_NormalPancreas[ tumor.idx ]  <- 1
  new.annot$CONTRAST_PancreaticCancer_vs_NormalPancreas[ normal.idx ]  <- (-1)
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE16515_GPL570 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the sample description
  new.annot$sample_description <- gsub(pattern = "tissue: ", "", annot$characteristics_ch1)
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_PancreaticCancer_vs_NormalPancreas
  ##---------------------------------------------------------------------------=
  ## Get the index for the different diseases' samples and normals:
  tumor.idx <- which(new.annot$sample_description == "Tumor Tissue in Pancreatic Cancer Sample")
  normal.idx <- which(new.annot$sample_description == "Normal Tissue in Pancreatic Cancer Sample")
  
  new.annot$CONTRAST_PancreaticCancer_vs_NormalPancreas <- 0
  new.annot$CONTRAST_PancreaticCancer_vs_NormalPancreas[ tumor.idx ]  <- 1
  new.annot$CONTRAST_PancreaticCancer_vs_NormalPancreas[ normal.idx ]  <- (-1)
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}


UpdateAnnotations_GSE12021_GPL96 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the sample description
  new.annot$sample_description <- annot$characteristics_ch1.2
  RA.idx <- grep(pattern = "rheumatoid arthritis", new.annot$sample_description)
  new.annot$sample_description[RA.idx] <- "rheumatoid arthritis"
  OA.idx <- grep(pattern = "osteoarthritis", new.annot$sample_description)
  new.annot$sample_description[OA.idx] <- "osteoarthritis"
  normal.idx <- grep(pattern = "normal", new.annot$sample_description)
  new.annot$sample_description[normal.idx] <- "normal control"
  
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_RA_vs_NormalSynovialMembrane
  ##  - CONTRAST_OA_vs_NormalSynovialMembrane
  ##---------------------------------------------------------------------------=
  
  new.annot$CONTRAST_RA_vs_NormalSynovialMembrane <- 0
  new.annot$CONTRAST_RA_vs_NormalSynovialMembrane[ RA.idx ]  <- 1
  new.annot$CONTRAST_RA_vs_NormalSynovialMembrane[ normal.idx ]  <- (-1)
  
  new.annot$CONTRAST_OA_vs_NormalSynovialMembrane <- 0
  new.annot$CONTRAST_OA_vs_NormalSynovialMembrane[ OA.idx ]  <- 1
  new.annot$CONTRAST_OA_vs_NormalSynovialMembrane[ normal.idx ]  <- (-1)
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}



UpdateAnnotations_GSE15573_GPL6102 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Keep the sample description
  new.annot$sample_description <- gsub(pattern = "status: ", replacement = "", annot$characteristics_ch1)
  
  ##---------------------------------------------------------------------------=
  ## Create the contrast:
  ##-----------=
  ##  - CONTRAST_RA_vs_Control_PBMC
  ##---------------------------------------------------------------------------=
  RA.idx <- which(new.annot$sample_description == "Rheumatoid Arthritis Patient")
  normal.idx <- which(new.annot$sample_description == "Control")
  
  new.annot$CONTRAST_RA_vs_Control_PBMC <- 0
  new.annot$CONTRAST_RA_vs_Control_PBMC[ RA.idx ]  <- 1
  new.annot$CONTRAST_RA_vs_Control_PBMC[ normal.idx ]  <- (-1)
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}






UpdateAnnotations_GSE23611_GPL6480 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Note - Data annotations appear to be well formatted.  No adjustments necessary!
  
  ## Get the index for the pre-treatment asthmatic samples, and for the lung normals:
  lung.biopsy.asthmatic.pre.treat <- grep(pattern = "human lung biopsy_asthmatic_pre", annot$source_name_ch1, fixed=T)
  lung.biopsy.normal <- which(annot$source_name_ch1 == "human lung biopsy_healthy")
  
  ## Contrast for asthma vs normal.
  new.annot$CONTRAST_Asthma_biopsy_v_normal <- 0
  new.annot$CONTRAST_Asthma_biopsy_v_normal[ lung.biopsy.asthmatic.pre.treat ]  <- 1
  new.annot$CONTRAST_Asthma_biopsy_v_normal[ lung.biopsy.normal ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}



UpdateAnnotations_GSE4302_GPL570 <- function(cur.eset) {
  ## Get the original pheno-data, and convert all of the factors into strings.
  annot <- notes(cur.eset)$original.pData
  annot <- data.frame(lapply(annot, as.character), stringsAsFactors = FALSE)
  new.annot <- data.frame(row.names=rownames(annot), stringsAsFactors = FALSE)
  ## Keep the geo-accession IDs for each entry just for unique identifiability and consistency.
  new.annot$geo_accession <- annot$geo_accession
  
  ## Note - Data annotations appear to be well formatted.  No adjustments necessary!
  
  ## Get the index for the asthmatic samples, and for the lung normals:
  lung.brushing.asthmatic.pre.treat <- which(annot$characteristics_ch1 == "Asthmatic at baseline")
  lung.brushing.normal <- which(annot$characteristics_ch1 == "Healthy control")
  
  ## Contrast for asthma vs normal.
  new.annot$CONTRAST_Asthma_brushing_v_normal <- 0
  new.annot$CONTRAST_Asthma_brushing_v_normal[ lung.brushing.asthmatic.pre.treat ]  <- 1
  new.annot$CONTRAST_Asthma_brushing_v_normal[ lung.brushing.normal ] <- (-1)
  
  
  pData(cur.eset) <- new.annot
  
  return(cur.eset)
}

