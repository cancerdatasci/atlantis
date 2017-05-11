#' @useDynLib atlantis

write.oos.predictions <- function(targetID, model, newdata, features, target.predictOnly, oos.fit.table.file) {
      oos.predictions <- predict(model, newdata=newdata)
      
      fit.table <- data.frame(target = target.predictOnly, pred_target=oos.predictions, newdata[,features, drop=F], row.names=rownames(newdata), check.names=F)
      colnames(fit.table)[1] <- targetID
      colnames(fit.table)[2] <- paste(targetID, 'pred', sep='_')

      write.csv(fit.table, file=oos.fit.table.file)
}

subset.target.and.features <- function(targetMat, predMat, predictOnly, targetID, predFeatureRegExps) {

  targetMat.full <- targetMat

  # remove any predictOnly datapoints from the targetMat
  targetMat <- targetMat.full[setdiff(rownames(targetMat.full), predictOnly),,drop=F]

  #check if targetID is index or column name, if index: make numeric and find corresponding column name
  if (grepl("^[[:digit:]]+$",targetID)){
    targetID <- as.numeric(targetID)
  }
  
  if (is.numeric(targetID)) {
    targetID <- colnames(targetMat)[targetID]
  }

  # coerce to dataframe so that we get NAs for missing cell lines
  predictOnly <- intersect(predictOnly, rownames(predMat))
  target.predictOnly <- as.matrix(as.data.frame(targetMat.full[, targetID,drop=F])[predictOnly,])

  targetVec <- targetMat[,targetID]
  names(targetVec) <- rownames(targetMat)
  
  if (is.null(names(targetVec))) {
    stop("targetVec does not have names")
  }

  #check if all cell lines in target mat are in pred mat
  #remove cell lines from target mat that are not in pred mat
  missingCellLines <- setdiff(names(targetVec),rownames(predMat))
  if (length(missingCellLines) != 0){
    print("removing cell lines from target matrix")
    keepCellLines <- setdiff(names(targetVec),missingCellLines)
    targetVec <- targetVec[keepCellLines]
  } else {
    keepCellLines <- names(targetVec)
  }
  
  # remove targetID if it appears in predMat
  if (targetID %in% colnames(predMat)) {
    cat(paste("Removing predictive feature", targetID, "from predMat as it is the target feature\n"))
    predMat <- predMat[, colnames(predMat) != targetID]
  }

  grep.with.warning <- function(pattern, predMatCols) {
    targetFeatures <- grep(pattern, predMatCols)
    if (length(targetFeatures) == 0) {
       print(length(predMatCols))
       if(length(predMatCols) > 10) {
          predMatCols <- c(predMatCols[1:10], "...")
       }
       cat(paste('WARNING: grep(', pattern,', ', do.call(paste, as.list(predMatCols)),") returned no results\n", sep=''))
    }

    targetFeatures
  }
  
  # filter predictive features to only those whose names are specified in predFeatureNamesToUse
  targetFeatures <- c()
  for (pattern in predFeatureRegExps) {
      matches <- grep.with.warning(pattern, colnames(predMat))
      #cat("for", pattern, "found", matches, "\n")
      targetFeatures <- c(targetFeatures, matches)
  }
  targetFeatures <- unique(targetFeatures)
  predMat <- predMat[, targetFeatures, drop=F]
  cat(paste("Using ", length(targetFeatures), "predictive features using patterns:"), paste(predFeatureRegExps, collapse=" "), "\n")

  predOnlyMat <- predMat[predictOnly, ,drop=F]
  if(ncol(predMat) > 0) {
    predMat <- predMat[names(targetVec), ,drop=F] # limit predMat to only cell lines of targetVec
  }

  predMat <- predMat[!is.na(targetVec),,drop=F]
  targetVec <- targetVec[!is.na(targetVec)]

  if(ncol(predMat) > 0) {
    # illegal column might have been created by removing the cell lines in previous line. remove them.
    predMat <- removeNAcols(predMat)
  }

  if(ncol(predMat) > 0) {
    sparse.columns <- columns.w.rare.values(predMat)
    predMat <- predMat[,!sparse.columns,drop=F]
  }

  if(ncol(predMat) > 0) {
    entire.row.is.na <- apply(predMat, 1, function(x) { all(is.na(x)) })
    predMat <- predMat[!entire.row.is.na,,drop=F]
    targetVec <- targetVec[!entire.row.is.na]
  }

  return(list(predMat=predMat, predOnlyMat=predOnlyMat, targetVec=targetVec, missingCellLines=missingCellLines, target.predictOnly=target.predictOnly))
}

columns.w.rare.values <- function(m) {
  apply(m, 2, function(x) { 
    max.freq <- max(table(x))
    max.freq >= sum(!is.na(x))-1
  })
}

#' Main entry point for invoking ATLANTIS
#' 
#' Fit a model to the provided data via ATLANTIS
#' @param analysis.ID A label to be used in the output files.  (For example "20130801")
#' @param targetID Either the number or the name of the feature in targetMat to use as response variable
#' @param predictOnly the names of the datapoints which should be only predicted, but not used for training
#' @param makePlot If TRUE, will generate a PDF report summarizing model
#' @param makeSummary If TRUE, will write out a file with summary of what ATLANTIS was run on
#' @param fitControlSettings Used to select a canned configuration for party.  (Use when biased sampling is desired)
#' @param predMat.file Path to Rdata file which contains features to use
#' @param targetMat.file Path to Rdata file which contains target values to predict
#' @param anno.file Path to annotation file containing information about features needed to generate plots
#' @param report.summary Any additional text to include in plot (to describe where the data came from)
#' @export
runATLANTIS <- function(
  analysis.ID, # e.g. "SAN_by_CN.MUT_20130306", "20130801"
  targetID, # either the number or the name of the feature in targetMat to use as response variable
  predictOnly=NULL, # the datapoints which should be only predicted, but not used for training
  makePlot = FALSE,
  makeSummary = FALSE,
  data.dir = paste('./', analysis.ID, '/data/', sep=''),
  output.dir = paste('./', analysis.ID, '/analysis/', sep=''),
  fitControlSettings = "default", 
  mode = "regression",
  additionalFeatures=NULL,
  predFeatureNamesToUse=NULL, # a list of the only genes whose predictive 
                                # features should be used
  predFeatureRegExps=NULL, # features which match this regexp will also be included
  predMat.file=file.path(data.dir, paste(analysis.ID, '_PredMat.AnnTable.Rdata', sep='')),
  targetMat.file=file.path(data.dir, paste(analysis.ID, '_targetMat.Rdata', sep='')),
  targetMat=NULL,
  anno.file=file.path(data.dir, 'annTable.txt'),
  kMinNumSamples=50,
  save.model=T,
  save.params=T,
  output.prefix=targetID,
  report.summary=c(),
  save.featureData=T
  ) {

  library(foreach)
  
  #get time of analysis
  timeOfRun <- Sys.time()
  
  load(predMat.file)
  
  if(is.null(targetMat)) {
    load(targetMat.file)
  }

  #ATLANTIS Summary Object
  ATLANTIS_Summary <- list()
  ATLANTIS_Summary$file.name <- file.path(output.dir,paste(output.prefix,'_ATLANTIS_Summary.Rdata',sep=''))
  ATLANTIS_Summary$analysisID <- analysis.ID
  ATLANTIS_Summary$targetID <- targetID
  ATLANTIS_Summary$time <- timeOfRun
  ATLANTIS_Summary$targetMat.file <- basename(targetMat.file)
  ATLANTIS_Summary$predMat.file <- basename(predMat.file)
  ATLANTIS_Summary$targetMatCellLines <- rownames(targetMat)
  ATLANTIS_Summary$predMatCellLines <- rownames(predMat)
  ATLANTIS_Summary$targetMatFeatures <- ncol(targetMat)
  ATLANTIS_Summary$predMatFeatures <- ncol(predMat)
  ATLANTIS_Summary$additionalFeatures <- additionalFeatures
  ATLANTIS_Summary$predFeatureNamesToUse <- predFeatureNamesToUse   
  ATLANTIS_Summary$predFeatureRegExps <- predFeatureRegExps

  # Translate gene names into regexps
  if(is.null(predFeatureRegExps) && is.null(predFeatureNamesToUse) ) {
    predFeatureRegExps <- ".*"
  } else {
    predFeatureRegExps <- c(predFeatureRegExps, sapply(predFeatureNamesToUse, function(featureName) { paste('(^|_|:)', featureName, '(:|_| |$)', sep='') } ))
  }

  subset <- subset.target.and.features(targetMat, predMat, predictOnly, targetID, predFeatureRegExps)
  targetVec <- subset$targetVec
  predMat <- subset$predMat
  predOnlyMat <- subset$predOnlyMat
  missingCellLines <- subset$missingCellLines
  target.predictOnly <- subset$target.predictOnly
  rm(targetMat)
  gc()

  stopifnot(all(!is.na(targetVec)))
  
  run.params <- list()
  run.params$cell.lines <- rownames(predMat)
  run.params$features <- colnames(predMat)
  run.params$targetID <- targetID
  
  run.params$file.name <- file.path(output.dir, 
                                    paste(analysis.ID, '_', targetID,
                                          '_run.params.Rdata', sep=''))
  if(save.params) {
    save(run.params, file=run.params$file.name)  
  }

  ATLANTIS_Summary$targetMatCellLinesANALYZED <- names(targetVec)
  ATLANTIS_Summary$missingTargetCellLines <- missingCellLines
  ATLANTIS_Summary$predMatCellLinesANALYZED <- rownames(predMat)  
  
  if(ncol(predMat) == 0) {
    model <- "no-features"
    fit.info <- list(failure.reason=model, targetVec=targetVec)
  } else if(! has.sufficent.sensitive.lines(targetVec, fitControlSettings)) {
    model <- "too-few-sensitive"
    fit.info <- list(failure.reason=model, targetVec=targetVec)  
  } else {
    stopifnot(all(dim(predMat) > 0))

    res <- fitOptimizedCForest(targetMat = as.matrix(targetVec),
                                 predMat = predMat, fitControlSettings = fitControlSettings)

    model <- res$model
    fit.info <- res$fitInfo
    if((is.character(model))) {
      fit.info$failure.reason = model
    }
  }
  
  fit.file <- file.path(output.dir, paste(output.prefix, '_fit.info_cforest.Rdata', sep=''))
  fit.info$report.summary <- c(
      fit.info$report.summary, 
      report.summary, 
      paste("# of inital features:", ncol(predMat)))

  if(save.featureData && !is.null(fit.info$predData)) {  
      save(fit.info, file=fit.file)

      fit.table.file <- file.path(output.dir, paste(output.prefix, '_fit.table.csv', sep=''))
      fit.table <- generatePredictorsTable(fit.info, ATLANTIS_Summary$targetID)
      write.csv(fit.table, file=fit.table.file)
  } else {
      # to save space, drop the feature data
      fit.info$predData <- NULL
      save(fit.info, file=fit.file)
      fit.table.file <- NULL
  }

  if((is.character(model))) {
    list(fit.successful=F,
         preprocess.successful=T,
         message=model,
         fit.file=fit.file,
         fit.table.file=fit.table.file)
  } else {
    if(save.model) { 
      model.file <- file.path(output.dir, paste(output.prefix, '_model_cforest.Rdata', sep=''))
      save(model, file=model.file)
    }

    if(length(predictOnly) > 0) {
      newdata<-as.data.frame(predOnlyMat)
      write.oos.predictions(targetID, model, newdata, names(fit.info$varImp), target.predictOnly, file.path(output.dir, paste(output.prefix, '_oos_pred.csv', sep='')))
    }
    
    ATLANTIS_Summary$metric <- fit.info$OOB$metric
    ATLANTIS_Summary$value <- fit.info$OOB$quality
    
    if(makeSummary || makePlot) {
      save(ATLANTIS_Summary, file=ATLANTIS_Summary$file.name)
    }
    
    if(makePlot){
      PlotATLANTISresults(fit.info, ATLANTIS_Summary$file.name, data.dir, output.dir, anno.file, output.prefix=output.prefix)
    }

    if(makeSummary) {
      WriteSummaryToFile(fit.info, ATLANTIS_Summary$file.name, output.dir)
    }
    
    list(fit.successful=T,
         preprocess.successful=T,
         message="Success",
         fit.file=fit.file,
         fit.table.file=fit.table.file)
  }  
}
