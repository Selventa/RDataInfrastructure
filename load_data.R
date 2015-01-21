
library("GEOquery")

source("~/GitHub/RDataInfrastructure/process_raw_data_funcs.R")
source("~/GitHub/RDataInfrastructure/load_data_funcs.R")
source("~/GitHub/RDataInfrastructure/update_annotations.R")


#PROBLEMS

# prob: GSE58795 uses custom affy exon array, so needs custom cdf
# soln: install custom cdf with filename to match pattern in process.data.affy()


# to resolve:  
GSE.IDs <- c("GSE42296", "GSE17755", "GSE33377")
# GSE.IDs <- ""
esets <- GetEset(GSE.IDs, 
                 eset.folder = "~/esets", 
                 overwrite.existing=T,
                 cache.folder = normalizePath("~/../Downloads"),
                 expt.annot=c(data.source = "from_raw", brainarray=T),
                 verbose=T)

# esets <- lapply(esets,
#                 map.features.to.EGID)
# 
# ds <- lapply(esets,
#              make.ds.from.eset)


load("~/esets/GSE17755_GPL1291.RData")

x <- lapply(esets, function(cur.eset) experimentData(cur.eset))
