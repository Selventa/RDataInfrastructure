
UpdateAnnotations_GSE42296_GPL6244 <- function(cur.eset) {
  
#   annot <- pData(cur.eset)
#   
#   baseline <- annot$characterisics_ch1.4 == "time: week 0"
#   Crohns <- grepl("Crohn", annot$characterisics_ch1.5)
#   annot$PRED_RA_inflix_response <- ifelse(baseline & !Crohns, 
#                                 ifelse(grepl("non-responder", annot$characterics_ch1.2, fixed = T),
#                                        "NR",
#                                        "R"),
#                                 NA)
#   pData(cur.eset) <- annot
  return(cur.eset)
  
}