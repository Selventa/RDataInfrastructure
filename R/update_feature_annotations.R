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


## The "Agilent Human 1A Oligo UNC custom Microarrays" chip.
UpdateFeatureAnnotations_GPL1390 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  ## Read in the current GPL1390 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    
    
    library(biomaRt)
    
    ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
    
#     tmp <- biomaRt::listAttributes(ensembl)
#     
#     possible.attributes <- grep("agilent", tmp$name, ignore.case=T, value=T)
#     
#     
#     possible.attributes <- c(possible.attributes)
#     annotations <- list()
#     for (cur.attribute in possible.attributes) {
#       annotations[[cur.attribute]] <- getBM(attributes=c(cur.attribute, "entrezgene"), mart=ensembl)
#       print(cur.attribute)
#       print(sum(annotations.df$SPOT_ID %in% annotations[[cur.attribute]][[cur.attribute]]))
#       
#     }

    tmp.annotations <- getBM(attributes=c("efg_agilent_wholegenome_4x44k_v1", "entrezgene"), mart=ensembl)
    egids.from.agilent <- tmp.annotations$entrezgene[match(annotations.df$SPOT_ID,
                                                           tmp.annotations$efg_agilent_wholegenome_4x44k_v1)]
    sum(!is.na(egids.from.agilent))
    
    tmp.annotations <- getBM(attributes=c("refseq_mrna", "entrezgene"), mart=ensembl)
    egids.from.refseq <- tmp.annotations$entrezgene[match(annotations.df$Refseq.ID,
                                               tmp.annotations$refseq_mrna)]
    sum(!is.na(egids.from.refseq))
    
    
    sum(!is.na(egids.from.refseq) | !is.na(egids.from.agilent))
    
    #look for differences
    diff.mapping <- 
      !is.na(egids.from.refseq) & 
      !is.na(egids.from.agilent) & 
      (egids.from.refseq != egids.from.agilent) 
    
    egids.from.refseq[diff.mapping] <- egids.from.agilent[diff.mapping] <- NA
    
    
    annotations.df$EGID <- ifelse(is.na(egids.from.refseq), 
                                  egids.from.agilent,
                                  egids.from.refseq)
    
    annotations.df$EGID_added_date <- as.character(Sys.Date(), format="%b %d %Y")

    ## Write out the updated annotation data frame back to CSV (replacing
    ##  the original annotation file if specified).
    if (replace.original) {
      write.csv(annotations.df, file = paste0(gpl.path, "/", gpl.filename), row.names = FALSE)
    } else {
      write.csv(annotations.df, 
                file = paste0(gpl.path, "/", 
                              gsub(pattern = ".csv$", 
                                   replacement = "_updated.csv", 
                                   gpl.filename, 
                                   ignore.case = TRUE)), 
                row.names = FALSE)
    }
    
    cat("Updated the GPL1390 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL1390 platform annotation file didn't require updating.\n")
  } 
  return(invisible(annotations.df))
}


##-----------------------------------------------------------------------------=
## GPL4133: "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F" chip.
##------=
## The annotation file doesn't contain an explicitly-named "EGID" column, but
##  does contain the EGIDs in it's "GENE" column, so simply add that column
##  annotation if not already added.
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL4133 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL4133 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    ## The EGID information is already contained in the "GENE" column.
    updated.columns <- cbind(annotations.df$GENE,  
                             Date=as.character(Sys.Date(), format="%b %d %Y"))
    colnames(updated.columns) <- c("EGID", "EGID_added_date")
    
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
    
    cat("Updated the GPL4133 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL4133 platform annotation file didn't require updating.\n")
  }
}



##-----------------------------------------------------------------------------=
## GPL7202: "Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Probe Name version)" chip.
##------=
## The annotation file doesn't contain an explicitly-named "EGID" column, but
##  does contain the EGIDs in it's "GENE" column, so simply add that column
##  annotation if not already added.
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL7202 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL7202 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    ## The EGID information is already contained in the "GENE" column.
    updated.columns <- cbind(annotations.df$GENE,  
                             Date=as.character(Sys.Date(), format="%b %d %Y"))
    colnames(updated.columns) <- c("EGID", "EGID_added_date")
    
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
    
    cat("Updated the GPL7202 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL7202 platform annotation file didn't require updating.\n")
  }
}






##-----------------------------------------------------------------------------=
## GPL13497: "Agilent-026652 Whole Human Genome Microarray 4x44K v2 (Probe Name version)" chip.
##------=
## The annotation file doesn't contain an explicitly-named "EGID" column, but
##  does contain the EGIDs in it's "GENE" column, so simply add that column
##  annotation if not already added.
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL13497 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL13497 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    ## The EGID information is already contained in the "GENE" column.
    updated.columns <- cbind(annotations.df$GENE,  
                             Date=as.character(Sys.Date(), format="%b %d %Y"))
    colnames(updated.columns) <- c("EGID", "EGID_added_date")
    
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
    
    cat("Updated the GPL13497 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL13497 platform annotation file didn't require updating.\n")
  }
}




##-----------------------------------------------------------------------------=
## GPL19943: "Agilent-027516 Human 180K Array (Probe name version)" chip.
##------=
## The annotation file doesn't contain an explicitly-named "EGID" column, but
##  does contain the EGIDs in it's "GENE" column, so simply add that column
##  annotation if not already added.
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL19943 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL19943 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    ## The EGID information is already contained in the "geneID" column.
    updated.columns <- cbind(annotations.df$geneID,  
                             Date=as.character(Sys.Date(), format="%b %d %Y"))
    colnames(updated.columns) <- c("EGID", "EGID_added_date")
    
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
    
    cat("Updated the GPL19943 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL19943 platform annotation file didn't require updating.\n")
  }
}





##-----------------------------------------------------------------------------=
## GPL4134: "Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Feature Number version)" chip.
##------=
## The annotation file doesn't contain an explicitly-named "EGID" column, but
##  does contain the EGIDs in it's "GENE" column, so simply add that column
##  annotation if not already added.
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL4134 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL4134 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    ## The EGID information is already contained in the "GENE" column.
    updated.columns <- cbind(annotations.df$GENE,  
                             Date=as.character(Sys.Date(), format="%b %d %Y"))
    colnames(updated.columns) <- c("EGID", "EGID_added_date")
    
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
    
    cat("Updated the GPL4134 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL4134 platform annotation file didn't require updating.\n")
  }
}




##-----------------------------------------------------------------------------=
## GPL10904: "Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Feature Number version)" chip.
##------=
## The annotation file doesn't contain an "EGID" column, but
##  does contain the gene symbol in it's "ORF" column
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL10904 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL10904 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  

  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {

    gene.names <- annotations.df$ORF
    
    
    library(org.Hs.eg.db)
    
    EGIDs <- select(org.Hs.eg.db, keys=gene.names, columns = "ENTREZID", keytype="SYMBOL")
    symbols.with.two.EGIDs <- EGIDs$SYMBOL[which(duplicated(EGIDs$SYMBOL))]
    # things that map to two EGIDs should be mapped to neither (since we can't disambiguate)
    EGIDs$ENTREZID[EGIDs$SYMBOL %in% symbols.with.two.EGIDs] <- NA
    EGIDs <- unique(EGIDs)
    stopifnot(setequal(EGIDs$SYMBOL, gene.names))
    
    # now for those that didn't map as the true gene symbol, see if they offer a single
    # match as a synonym
    syn.EGIDs <- select(org.Hs.eg.db, keys=gene.names[is.na(EGIDs$ENTREZID)],
                        columns="ENTREZID",
                        keytype="ALIAS")
    
    # get rid of entries that aren't aliases for any EGID
    syn.EGIDs <- syn.EGIDs[!is.na(syn.EGIDs$ENTREZID),]
    
    # get rid of entries for which we already have a mapping (ie, we have the true
    # gene symbol for that EGID)
    syn.EGIDs <- syn.EGIDs[!(syn.EGIDs$ENTREZID %in% EGIDs$ENTREZID),]
    
    # get rid of gene names that are aliases for more than one EGID
    genes.with.dup.EGIDs <- syn.EGIDs$ALIAS[duplicated(syn.EGIDs$ALIAS)]
    syn.EGIDs <- syn.EGIDs[!(syn.EGIDs$ALIAS %in% genes.with.dup.EGIDs),]
    
    # get rid of EGIDs that map to multiple aliases in the list
    EGIDs.with.two.symbols <- syn.EGIDs$ENTREZID[which(duplicated(syn.EGIDs$ENTREZID))]
    syn.EGIDs <- syn.EGIDs[!(syn.EGIDs$ENTREZID %in% EGIDs.with.two.symbols),]
    
    # make sure we've got a 1-to-1 mapping    
    stopifnot(!any(duplicated(syn.EGIDs$ALIAS)))
    stopifnot(!any(duplicated(syn.EGIDs$ENTREZID)))
    
    # put the remaining EGIDs into the EGIDs object
    EGIDs$ENTREZID[match(syn.EGIDs$ALIAS, EGIDs$SYMBOL)] <- syn.EGIDs$ENTREZID
    
    
    ## Add the EGIDs and the current date onto the annotation data frame.
    updated.annotation.df <- data.frame(annotations.df,
                                         EGID = as.numeric(EGIDs$ENTREZID[match(annotations.df$ORF, EGIDs$SYMBOL)]),
                                         EGID_added_date = as.character(Sys.Date(), format="%b %d %Y"),
                                        stringsAsFactors=F)
    
    
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
    
    cat("Updated the GPL10904 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL10904 platform annotation file didn't require updating.\n")
  }
}




##-----------------------------------------------------------------------------=
## GPL18943: "	NimbleGen Human Gene Expression 12x135K Array [100718_HG18_opt_expr]" chip.
##------=
## The annotation file doesn't contain an explicitly-named "EGID" column, but
##  does contain the EGIDs in it's "GENE" column, so simply add that column
##  annotation if not already added.
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL18943 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL18943 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    ## The EGID information is already contained in the "DESCRIPTION" column, but
    ## we need to strsplit to get it out.
    
    desc.split <- strsplit(annotations.df$DESCRIPTION, ";")
    gene.id.text <- lapply(desc.split,
                          function(cur.vec) {
                            grep("gene_id", cur.vec, value=T)
                            })
    table(sapply(gene.id.text, length))
    # ok, gene_id appears zero or one times per probe...
    
    
    gene.id <- sapply(gene.id.text,
                      function(x) {
                        if (length(x)==0) {
                          return(NA_integer_)
                        } else {
                          stopifnot(length(x)==1)
                          cur.id <- sub(" gene_id ", "", x, fixed=T)
                          cur.id <- gsub("'", "", cur.id, fixed=T)
                          if (is.na(as.numeric(cur.id))) {
                            browser()
                            print(cur.id)
                          }
                          return(as.numeric(cur.id))
                        }
                      })
    
    
    ## Add the EGIDs and the current date onto the annotation data frame.
    updated.annotation.df <- cbind(annotations.df, 
                                   data.frame(EGID=gene.id,
                                              EGID_added_date=as.character(Sys.Date(), format="%b %d %Y")))
    
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
    
    cat("Updated the GPL18943 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL18943 platform annotation file didn't require updating.\n")
  }
}




##-----------------------------------------------------------------------------=
## GPL19462: "	Illumina MouseRef-8 v2.0 expression beadchip (ILMN_SYMBOL)" chip.
##------=
## The annotation file doesn't contain an "EGID" column, but
##  does contain the gene symbol in it's "ORF" column
##-----------------------------------------------------------------------------=
UpdateFeatureAnnotations_GPL19462 <- function( gpl.filename, gpl.path = getwd(), replace.original=TRUE ) {
  
  ## Read in the current GPL19462 platform annotation CSV file.  
  annotations.df <- read.csv(file = paste0(gpl.path, "/", gpl.filename), stringsAsFactors = FALSE)
  
  
  ## If there is no EGID (Entrez gene ID) column, then we will add one.
  if (!("EGID" %in% colnames(annotations.df))) {
    
    gene.names <- annotations.df$ORF
    
    
    library(org.Mm.eg.db)
    
    EGIDs <- select(org.Mm.eg.db, keys=gene.names, columns = "ENTREZID", keytype="SYMBOL")
    symbols.with.two.EGIDs <- EGIDs$SYMBOL[which(duplicated(EGIDs$SYMBOL))]
    # things that map to two EGIDs should be mapped to neither (since we can't disambiguate)
    EGIDs$ENTREZID[EGIDs$SYMBOL %in% symbols.with.two.EGIDs] <- NA
    EGIDs <- unique(EGIDs)
    stopifnot(setequal(EGIDs$SYMBOL, gene.names))
    
    # now for those that didn't map as the true gene symbol, see if they offer a single
    # match as a synonym
    syn.EGIDs <- select(org.Mm.eg.db, keys=gene.names[is.na(EGIDs$ENTREZID)],
                        columns="ENTREZID",
                        keytype="ALIAS")
    
    # get rid of entries that aren't aliases for any EGID
    syn.EGIDs <- syn.EGIDs[!is.na(syn.EGIDs$ENTREZID),]
    
    # get rid of entries for which we already have a mapping (ie, we have the true
    # gene symbol for that EGID)
    syn.EGIDs <- syn.EGIDs[!(syn.EGIDs$ENTREZID %in% EGIDs$ENTREZID),]
    
    # get rid of gene names that are aliases for more than one EGID
    genes.with.dup.EGIDs <- syn.EGIDs$ALIAS[duplicated(syn.EGIDs$ALIAS)]
    syn.EGIDs <- syn.EGIDs[!(syn.EGIDs$ALIAS %in% genes.with.dup.EGIDs),]
    
    # get rid of EGIDs that map to multiple aliases in the list
    EGIDs.with.two.symbols <- syn.EGIDs$ENTREZID[which(duplicated(syn.EGIDs$ENTREZID))]
    syn.EGIDs <- syn.EGIDs[!(syn.EGIDs$ENTREZID %in% EGIDs.with.two.symbols),]
    
    # make sure we've got a 1-to-1 mapping    
    stopifnot(!any(duplicated(syn.EGIDs$ALIAS)))
    stopifnot(!any(duplicated(syn.EGIDs$ENTREZID)))
    
    # put the remaining EGIDs into the EGIDs object
    EGIDs$ENTREZID[match(syn.EGIDs$ALIAS, EGIDs$SYMBOL)] <- syn.EGIDs$ENTREZID
    
    
    ## Add the EGIDs and the current date onto the annotation data frame.
    updated.annotation.df <- data.frame(annotations.df,
                                        EGID = as.numeric(EGIDs$ENTREZID[match(annotations.df$ORF, EGIDs$SYMBOL)]),
                                        EGID_added_date = as.character(Sys.Date(), format="%b %d %Y"),
                                        stringsAsFactors=F)
    
    
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
    
    cat("Updated the GPL19462 platform annotation file to include EGIDs.\n")
  } else {
    cat("The GPL19462 platform annotation file didn't require updating.\n")
  }
}


