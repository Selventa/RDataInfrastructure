library(limma)

GEOPipelineLimma <- function(from.raw.num = 1, eset.list = esets){
  #This function will run the appropriate limma functions on an expression set list generated through the GEO pipeline.
  #typeof(from.raw.num) == integer/double. This gives the index in the eset.list where the expressionset of interest is 
  #located.
  #typeof(eset.list) == list. This should be a list of eSet items. The eset of interest should have a specially modified
  #pData attribute. This modification should be made through an UpdateAnnotations script.
  #One column in pData should start with DESIGN_ preceding the column name. This indicates that this column holds the
  #information that is the reason for these contrasts (disease type, tissue type, etc).
  #One or more columns in pData should start with CONTRAST_.These should have 1's indicating which rows are included 
  #in the numerator of the contrast and -1's indicating which rows are included in the denominator of that contrast. 
  #There should only be one contrast marked per row.
  
  #locate the design column
  design.col <- colnames(pData(eset.list[[from.raw.num]]))[grepl("DESIGN_", colnames(pData(eset.list[[from.raw.num]])))]
  stopifnot(length(design.col)==1)
  #locate columns which delineate contrasts
  contrast.cols <- colnames(pData(eset.list[[from.raw.num]]))[grepl("CONTRAST_", colnames(pData(eset.list[[from.raw.num]])))]
  #create model matrix and linear fit for this eset
  design <- eval(parse(text = paste0("model.matrix(~0+", design.col,", pData(eset.list[[", from.raw.num, "]]))", sep="")))
  fit <- lmFit(eset.list[[from.raw.num]], design)
  eFits <- list()

  #run each contrast and save eBayes fit and top differential probesets to a list
  for (contrast in 1:length(contrast.cols)){
    contrast.matrix <- matrix(0, nrow=ncol(design), ncol=1, dimnames=list(colnames(design), "Contrast Value"))
    contrast.numerator <- pData(eset.list[[from.raw.num]])[pData(eset.list[[from.raw.num]])[,contrast.cols[contrast]]==1, design.col]
    stopifnot(length(unique(contrast.numerator)) == 1)
    contrast.denominator <- pData(eset.list[[from.raw.num]])[pData(eset.list[[from.raw.num]])[,contrast.cols[contrast]]== -1, design.col]
    stopifnot(length(unique(contrast.denominator)) == 1)
    contrast.matrix[unique(contrast.numerator), 1] <- 1
    contrast.matrix[unique(contrast.denominator), 1] <- -1
    tmp.fit <- contrasts.fit(fit, contrast.matrix)
    tmp.eFit <- eBayes(tmp.fit)
    eFits[[contrast]] <- tmp.eFit
    names(eFits)[contrast] <- contrast.cols[contrast]
  }
  return(eFits)
}
eFits <- GEOPipelineLimma()
