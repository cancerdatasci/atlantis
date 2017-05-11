#' Generate ATLANTIS plots from intermediate files
#'
#' @param fi fit.info object (or file name)
#' @param ATLANTIS_Summary the structure itself or file to load it from
#' @param code.dir unused -- to be dropped 
#' @param output.dir if output.dir is NULL then the plot is written to the current device, otherwise a PDF is written to the output.dir
#' @param annTable.file File containing feature annotations (the gene, location, etc)
#' @param output.prefix a prefix to prepend to all output files
#' @export
PlotATLANTISresults <- function(fi, ATLANTIS_Summary, code.dir, output.dir, annTable.file, output.prefix = NULL){
  
  if(is.character(fi) == TRUE){
    load(fi)
  }

  if(is.character(ATLANTIS_Summary) == TRUE){
    load(ATLANTIS_Summary)
  }

  #read in annTable file stored in libdir
  
  annTable <- read.table(annTable.file,header=TRUE,comment.char="",sep="\t", fill=FALSE,stringsAsFactors = FALSE, check.names = FALSE,row.names=1)
  
  #get analysisID, targetID, targetVec, predMat from ATLANTIS_Summary
  targetVec <- fi$targetVec
  predMat <- fi$predData
  targetID <- ATLANTIS_Summary$targetID
  analysisID <- ATLANTIS_Summary$analysisID
  additionalFeatures <- ATLANTIS_Summary$additionalFeatures

  if(is.null(output.prefix)) {
    output.prefix = targetID
  }
  
  # create plots
  if(!is.null(output.dir)) {
    plots.file <- paste(output.prefix, '_plots.pdf', sep='') 
    pdf(file.path(output.dir, plots.file), w=11*1.5, h=8.5*1.5)
  }
 
  report.summary <- fi$report.summary

  for (nominalScale in c(F, T)) {
    for (orderByPred in c(F,T)) {
      tmp <- plots.plotModel(fi, targetVec, targetID, 
                         predMat, annTable,
                         maxNumPreds=20,
                         report.summary,
                         orderByPred=orderByPred,
                         additionalFeatures=additionalFeatures, nominalScale=nominalScale)
    } 
  }

  if(!is.null(output.dir)) {
    dev.off()
  }
}

# draw text "msg" scaling to fit in the specified bounding box
drawBoundedText <- function(msg, left, bottom, right, top) {
  target.width <- right - left;
  target.height <- top - bottom;
  scale <- min(target.width/strwidth(msg), target.height/strheight(msg))
  text(left, bottom, msg, cex=scale, adj=c(0, 0))
}

# generates a plot with a summary of predictive variables used in gbm model
# fi - fit.info of model
# targetVec - target vector (from targetMat)
# targetID - the ID of the target variable (its name in targetMat)
# maxNumPreds - number of features to plot
plots.plotModel <- function(fi, targetVec, targetID, predMat, annTable, maxNumPreds, report.summary, ...) {
  library(caTools) # for colAUC
  
  # top 20 features, ordered by cluster & contribution of their clusters
  # head(used.features.inf[order(used.features.inf[1:20,"clusterRank"]), ], 20)
  
  # can move the following to model creation GBM too
  predProfile <- fi$OOB$prediction
  
  model_quality <- list()
  model_quality$value <- fi$OOB$quality
  model_quality$metric <- fi$OOB$metric
  
  model.features.ordered <- fi$varImp
  
# convert targetVec to numeric binary vector if is factor
  if (fi$modelType == "Classification") {
    targetVec <- ifelse(as.numeric(targetVec) == 1, 0, 1)  
  }

  plots.predictorsBars(targetVec, predMat, model.features.ordered, 
                       targetID, model_quality, annTable=annTable,
                       report.summary=report.summary,
                       predProfile=predProfile, 
                       maxNumPreds=maxNumPreds, ...)
}




# model.features - relative influence for predictors.  a named vector
# predProfile - predictions of the targetVec, to be plotted on-top
# model_quality - a list. $metric is a string ("R^2" or "AUC"), $value holds the value
plots.predictorsBars <- function(depProfile, predMat, model.features, gsID, model_quality,
                                 annTable, report.summary,
                                 additionalFeatures=NULL, maxNumPreds=10, predProfile=NULL,
                                 orderByPred=FALSE, nominalScale=F) {

  # subset influence table to include only top predictive features (whose influence sum to 50)
  ## order.by.inf <- order(model.features[, "rel.inf"], decreasing=T)
  ## model.features <- model.features[order.by.inf, ]
  
  # how many features are used by the model?
  numModelFeatures <- length(model.features)
  
  #lastPredToPlot <- min(max(which(cumsum(model.features[,"rel.inf"]) < 50)) + 1, maxNumPreds)
  lastPredToPlot <- min(numModelFeatures, maxNumPreds)
  
  
  # we're going to plot only these predictors
  features.to.plot <- model.features[1:lastPredToPlot]
  stopifnot(all(names(features.to.plot) %in% colnames(predMat)))
  feat.data <- predMat[, names(features.to.plot), drop=F]
  
  # plot barplots for all predictors, target GS, influence weights too
  num.feats <- length(features.to.plot)
  if (orderByPred) {
    samp.ord <- order(predProfile, depProfile)  # the order of the samples in the plots
  } else {
    samp.ord <- order(depProfile, predProfile)  # order by true then prediction
  }
  layout(matrix(c(num.feats+4, rep(1, num.feats+2), 2, num.feats+5, 3:(num.feats+3)), ncol=2), 
         width=c(1, 5),
         heights=c(3, 0.5, rep(1, num.feats), 1.5))
  
  # plot rel.inf values
  op <- par(mar = c(3, 2, 0, 0))
  
  featCorToTarget <- cor(feat.data, 
                         as.numeric(depProfile), 
                         use="pairwise.complete.obs", 
                         method='spearman')
  cols <- rev(c("goldenrod", "deeppink3")[as.numeric(featCorToTarget > 0)+1])
  
  barplot(rev(features.to.plot), horiz=T, col=cols, names.arg="")
  par(op)
  
  ######
  # plot response variable
  ######
  # figure out chr loci for response varibale
  responseHUGO <- annTable[gsID, "HUGOsymbol"]
  responseChrLoci <- annTable[gsID, "ChrLoci"]
  gsLabel <- paste(gsID, "\n", responseChrLoci, sep='')
  
  cols <- ifelse(depProfile[samp.ord] < -2.0, "red", "grey")
  op <- par(mar = c(0, 2, 1, 13)) # c(bottom, left, top, right) 
  mp <- barplot(depProfile[samp.ord], names.arg="", axes=T, col=cols, border=NA)
  points(mp, predProfile[samp.ord], pch='+', col="palegreen4", type='l', lwd=2)
  mtext(gsLabel, side=4, las=1, cex=0.9)
  par(op)
  
  # plot barplots for predictive features
  uniqChrLoci <- unique(annTable[names(features.to.plot), "ChrLoci"], border=NA)
  chrBands <- sub("\\..+$", "", uniqChrLoci)
  names(chrBands) <- uniqChrLoci
  
  colMap <- plots.colorChrLoc(chrBands[annTable[names(features.to.plot), "ChrLoci"]])
  for (i in 1:num.feats) {
    featID <- names(features.to.plot)[i]
    plotdata <- feat.data[samp.ord, i]
    if (annTable[featID, "VarType"] == "CN") {
      # for CN features
      # color amps in red, dels in blue
      cols <- ifelse(plotdata > 0.5, "red", "grey")
      cols[plotdata < -0.5] <- "blue"
      ylim=c(-1,1)
    }
    if (annTable[featID, "VarType"] %in% c("MUT", "Mmis", "Mbad", "Mcos", "Mall", "RsMmis", "RsMnon")) {
      # color mutations in purple
      f <- feat.data[samp.ord, i]
      stopifnot(all(f %in% c(0,1,2) | is.na(f)))
      cols <- c("grey", "sienna4", "darkorange")[f+1]
      ylim=c(0,2)      
    }
    if (annTable[featID, "VarType"] %in% c("ACH", "GS", "SS")) {
      # color high dependency in red
      cols <- ifelse(feat.data[samp.ord, i] > 0.44, "orangered3", "grey")
      cols[plotdata < -0.44] <- "seagreen3"
      ylim=c(-2,2)      
    }
    if (annTable[featID, "VarType"] %in% c("Exp", "GSE", "CTD", "mirE")) {
      # transform data to robust Z
      z.plotdata <- (plotdata - median(plotdata, na.rm=T)) / mad(plotdata, na.rm=T)
      if(!nominalScale) {
        plotdata <- z.plotdata
      }
      # color high expression
      cols <- ifelse(z.plotdata > 2, "orangered3", "grey")
      cols[z.plotdata < -2] <- "seagreen3"
      ylim=c(-3,3)      
    }
    if (annTable[featID, "VarType"] == "SI") {
      if (length(unique(feat.data[,i])) <= 3)  { # binary + NA
        cols <- ifelse(feat.data[samp.ord, i] > 0, "orange3", "grey")
        ylim=c(0,1)      
      } else {
        ylim <- NULL
        cols <- NULL
      }
    }
    if (annTable[featID, "VarType"] == "San") {
      # transform data to robust Z
      plotdata <- (plotdata - median(plotdata, na.rm=T)) / mad(plotdata, na.rm=T)      
      # color high expression
      cols <- ifelse(plotdata > 2, "orangered3", "grey")
      cols[plotdata < -2] <- "seagreen3"
      ylim=c(-3,3)      
#       # color high sensitivity in red
#       cols <- ifelse(feat.data[samp.ord, i] < 0.3, "orangered3", "grey")
#       ylim=c(0,1)      
    }

feature.barplot <- function(plotdata, axes.range, cols, ylim) {
  #par(oma=rep(1,4))
  #par(oma=c(1,0,1,0))
  orig.mar <- par()$mar
  par(mar=c(1,0,1,0)+orig.mar) 
  par(mgp=c(0,-2,-3))

  r.range <- round( c(min(plotdata,na.rm=T), max(plotdata,na.rm=T)), 0)
  if(is.null(ylim)) {
    ylim <- c( min(c(plotdata, r.range[1]), na.rm=T), max(c(plotdata, r.range[2]), na.rm=T) )
  }

  if(any(is.infinite(ylim))) {
     # hack -- ylim is infinite in some case. Need to investigate where how this happens
     ylim <- c(-1e10, 1e10)
     warning("ylim was infinite so skipping this bar")
  } else {

    barplot(plotdata, names.arg="", axes=F, col=cols, ylim=ylim, border=NA)
    if(!is.null(axes.range)) {
      axes.range <-round( c(min(plotdata,na.rm=T), max(plotdata,na.rm=T)), 0)
  #    axis(ifelse((i %% 2) == 1,2,4), axes.range, las=2)
      axis(2, axes.range, las=2)
    }
  }

  par(mar=orig.mar)
}
    
    op <- par(mar = c(0, 2, 0, 13)) # c(bottom, left, top, right)
    axes.range <- NULL
#    if(nominalScale) { 
      ylim <- NULL 
      axes.range = pretty(plotdata, n=1)
#    }
    feature.barplot(plotdata, axes.range, cols, ylim)

    if (annTable[featID, "VarType"] %in% c("GSE", "SI")) {
      first.line.text <- substr(featID,1,15)
      second.line.text <- paste(substr(featID,16,32), "\n", substr(featID,33,49), sep='')
    } else {
      first.line.text <- substr(featID,1,25)
      second.line.text <- annTable[featID, "ChrLoci"]    
    }
    
    feature.label <- paste(first.line.text, "\n", second.line.text, sep='')
    mtext(feature.label, side=4, las=1, col=colMap[chrBands[annTable[featID, "ChrLoci"]]], cex=0.9) 
    par(op)
  }
  
# plot text on top left corner
  op <- par(mar = rep(0, 4))
  plot.new()
#  text(0.85, 0.5, pos=4, paste(report.summary, collapse="\n"))
  par(op)
  
  # model info in the top left corner
  op <- par(mar = rep(0, 4))
  plot.new()
  
  drawBoundedText(paste("CV", model_quality$metric, "=", format(model_quality$value, digits=2)), 0, 0.9, 1, 1)
  drawBoundedText(paste(report.summary, collapse="\n"), 0, 0, 1, 0.85)
  par(op)
  
  list(feat.data=feat.data, samp.ord=samp.ord)
}


# given a vector of strings containing chromosomal locations, assign
# a color to each one so that chromosomal loci that are represented more
# than once get colored in the same color.
plots.colorChrLoc <- function(chr.loci, def.color="black") {
  library(RColorBrewer)
  
  freqs <- table(chr.loci)
  num.unique <- length(freqs)
  cols <- rep(def.color, num.unique)
  names(cols) <- names(freqs)
  
  repeats <- names(freqs)[freqs >= 2]
  num.repeats <- length(repeats)
  cols[repeats] <- brewer.pal(8, "Set1")[1:length(repeats)]
  
  cols
}

