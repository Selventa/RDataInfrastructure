library(limma)

GEOPipelineLimma <- function(relev.eset){
  #This function will run the appropriate limma functions on an expression set list generated through the GEO pipeline.
  #typeof(relev.eset) == expressionSet. This should be the eset of interest. 
  #If using LimmaGeneric, read the following:
  #Atleast one column in pData should start with DESIGN_ preceding the column name. This indicates that this column holds the
  #information that is the reason for these contrasts (disease type, tissue type, etc). If there
  #are multiple design columns, they should start with DESIGN_1, DESIGN_2, etc.
  #One or more columns in pData should start with CONTRAST_.These should have 1's indicating which rows are included 
  #in the numerator of the contrast and -1's indicating which rows are included in the denominator of that contrast. 
  #There should only be one contrast marked per row. If there are multiple design columns, the
  #contrasts should be annotated so that they can be matched with the relevant design column.
  #Example: design columns - DESIGN_1disease, DESIGN_2gender
  #contrast columns -- CONTRAST_1crohns_vs_normal, CONTRAST_1ibd_vs_normal, CONTRAST_1UC_vs_normal,
  #CONTRAST_2male_vs_female
  
  #check if there is a specific limma function for this dataset
  relev.eset.limma.func <- paste0("Limma_", notes(relev.eset)$name)
  if (exists(relev.eset.limma.func)) {
    limma.func.name <- relev.eset.limma.func
  }else{limma.func.name <- "LimmaGeneric"}
  #if there is a specific function, run it; if not, run generic
  eFits <- eval(parse(text=paste0(limma.func.name,"(relev.eset)")))
  return(eFits)
}

LimmaGeneric <- function(relev.eset){
  #typeof(relev.eset) == expressionSet. This is the eset of interest.
  
  #locate design columns
  design.cols <- grep("^DESIGN_", colnames(pData(relev.eset)), value=TRUE)
  #locate columns which delineate contrasts
  contrast.cols <- grep("^CONTRAST_", colnames(pData(relev.eset)), value=TRUE)
  #find the numbers associated with the designs and contrasts
  num.design.cols <- unique(na.omit(as.numeric(unlist(strsplit(design.cols, "[^[:digit:]]")))))
  num.contrast.cols <- unique(na.omit(as.numeric(unlist(strsplit(contrast.cols, "[^[:digit:]]")))))
  #if the same numbers (order irrespective) are not found in both the design columns and contrast
  #columns, end the analysis (updateannotations error)
  stopifnot(all(num.design.cols %in% num.contrast.cols) &
              length(num.design.cols) == length(num.contrast.cols))
  #if there is no numeric annotation (implying that there is only one design column), change the
  #num.design.cols to be a vector of 1 so that the for loop below still applies
  if (length(num.design.cols) == 0 & length(num.contrast.cols) == 0){
    num.design.cols = c(1)
  }
  
  
  eFit.res <- list()
  for (comparison in num.design.cols){
    #create model matrix and linear fit for this design column
    #if there are multiple design columns, then select only the one currently being considered in 
    #the loop
    if (length(num.contrast.cols) > 1){
      design.col <- grep(paste0("DESIGN_", comparison), design.cols, value=TRUE)
      relev.contrast.cols <- grep(paste0("CONTRAST_", comparison), contrast.cols, value = TRUE)
    } else{
      design.col <- design.cols
      relev.contrast.cols <- contrast.cols
    }
    stopifnot(length(design.col)==1)
    design <- eval(parse(text = paste0("model.matrix(~0+", design.col,", pData(relev.eset))", sep="")))
    #remove odd column naming convention (adds design.col to the beginning of the types in design.col)
    colnames(design) <-  sub(design.col, "", colnames(design))
    fit <- lmFit(relev.eset, design)
    eFit.list <- list()
    
  
  #run each contrast using RunContrasts; save eBayes fit to a list
    for (contrast in 1:length(relev.contrast.cols)){
      contrast.eFit <- tryCatch({
        RunContrasts(relev.contrast.cols[contrast], design.col, design, relev.eset, fit)
      }, error = function(err){
        print(paste("Error:", err))
        print(paste("An error occurred. Contrast ", relev.contrast.cols[contrast], "will not be calculated."))
        return(err)
      })
      eFit.list[[contrast]] <- contrast.eFit
      names(eFit.list)[contrast] <- relev.contrast.cols[contrast]
    }
  
    #add results from each comparison to the final results list
    eFit.res[[comparison]] <- eFit.list
    names(eFit.res)[[comparison]] <- design.col
  
  }
  return(eFit.res)
}

RunContrasts <- function(contrast.column, design.column, design.mat, relev.eset, linfit){
  #typeof(contrast.column) == string. This is the name of the column that holds the designation for which samples should
  #be compared in the contrast.
  #typeof(design.column) == string. This is the name of the column that holds the relevant design information (see docstring
  #for GEOPipelineLimma for more information).
  #typeof(design.mat) == matrix. This is the model/design matrix for the linear fit model for the eset of interest.
  #typeof(relev.eset) == expressionSet. This is the eset of interest.
  
  contrast.matrix <- matrix(0, nrow=ncol(design.mat), ncol=1, 
                            dimnames=list(colnames(design.mat), "Contrast Value"))
  contrast.numerator <- pData(relev.eset)[pData(relev.eset)[,contrast.column]==1, design.column]
  #ensures that there is no user error in selecting the samples for contrast; they should all be of the same type in the
  #design column
  contrast.denominator <- pData(relev.eset)[pData(relev.eset)[,contrast.column]== -1, design.column]
  if(length(unique(contrast.numerator)) != 1){
    stop("Error: More or less than one design type among contrast numerator samples")
  }
  if(length(unique(contrast.denominator)) != 1){
    stop("Error: More or less than one design type among contrast denominator samples")
  }
  contrast.matrix[unique(contrast.numerator), 1] <- 1
  contrast.matrix[unique(contrast.denominator), 1] <- -1
  tmp.fit <- contrasts.fit(linfit, contrast.matrix)
  tmp.eFit <- eBayes(tmp.fit)
  return(tmp.eFit)
}
