WriteSummaryToFile <- function(fit.info, ATLANTIS_Summary, output.dir){
  
  library(caTools)
  
  if(is.character(fit.info) == TRUE){
    load(fit.info)
  } 
  if(is.character(ATLANTIS_Summary) == TRUE){
    load(ATLANTIS_Summary)
  }
  
  #ATLANTIS Summary output file name
  SummaryFileName <- file.path(output.dir, paste(paste("ATLANTIS_Summary",gsub(" ","_",ATLANTIS_Summary$time),ATLANTIS_Summary$targetID,sep="_"),".txt",sep=""))
  
  #Write to ATLANTIS Summary output file
  cat(c("ATLANTIS Summary",paste("Date and Time of Analysis: ",ATLANTIS_Summary$time),"",paste("Analysis ID: ",ATLANTIS_Summary$analysisID),"",paste("The target matrix file, ",ATLANTIS_Summary$targetMat.file,", contains ",ATLANTIS_Summary$targetMatFeatures," feature(s). The feature ",ATLANTIS_Summary$targetID, " was analyzed during this run of ATLANTIS.\n"),paste("The predictive matrix file, ",ATLANTIS_Summary$predMat.file,", contains ",ATLANTIS_Summary$predMatFeatures, " feature(s).\n")),file = SummaryFileName,sep='\n',append=FALSE)
  
  #finish Writing to ATLANTIS Summary Output File  
  summaryMatrix <- matrix(c(length(ATLANTIS_Summary$targetMatCellLines),length(ATLANTIS_Summary$targetMatCellLinesANALYZED),length(ATLANTIS_Summary$predMatCellLines),length(ATLANTIS_Summary$predMatCellLinesANALYZED)),nrow=2,ncol=2,byrow=TRUE,dimnames=list(c("Target Matrix","Predictive Matrix"),c("Number of Cell Lines Inputted","Number of Cell lines Analyzed")))
  
  capture.output(print(summaryMatrix,print.gap=2),file=SummaryFileName,append=TRUE)
  
  predProfile <- fit.info$OOB$prediction
    
  metric <- fit.info$OOB$metric
  value <- fit.info$OOB$quality
  modelType <- fit.info$modelType
  
  if (is.null(ATLANTIS_Summary[["pvalue"]]) == FALSE){
    cat(c(paste("\nMode: ",modelType),paste("Metric: ",metric),paste("Value: ",value),paste("P-value: ",ATLANTIS_Summary$pvalue)),file=SummaryFileName,append=TRUE)
  } else{
    cat(c(paste("\nMode: ",modelType),paste("Metric: ",metric),paste("Value: ",value)),file=SummaryFileName,append=TRUE)
  }
    
  
  cat(c("\nCell lines that were not analyzed: \n"),file=SummaryFileName,append=TRUE)
  write(ATLANTIS_Summary$missingTargetCellLines,SummaryFileName,ncolumns=1,append=TRUE,sep='\n')
  write("\n\nCell lines that were analyzed: \n",file=SummaryFileName,append=TRUE)
  write(ATLANTIS_Summary$targetMatCellLinesANALYZED,SummaryFileName,ncolumns=1,append=TRUE,sep='\n') 
  
}