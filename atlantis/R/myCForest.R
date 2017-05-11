#' @import atlantisPartyMod

compress.via.pack <- function(tree) { 
  fs <- .Call("R_packNode", tree, PACKAGE="atlantis")
  cat("Calling pack\n");
#  if(exists("mem_used")){
#    print(mem_used());
#  }
  list(fs)
}

expand_count <- 0
expand.via.unpack <- function(fs) {
  tree <- .Call("R_unpackNode", fs[[1]], PACKAGE="atlantis")
  tree
}

myRF <- atlantisPartyMod:::RandomForest

# altered version of RandomForest's dpp()
myRF@dpp  <- function (predMat, targetVec, subset = NULL, na.action = NULL, 
                       xtrafo = ptrafo, ytrafo = ptrafo, scores = NULL, ...) 
{
  dat <- ModelEnvMatrix(designMatrix = predMat, responseMatrix = targetVec,
                        subset = subset, ...)
  inp <- initVariableFrame(predMat, trafo = xtrafo, 
                           scores = scores)
  response <- as.data.frame(targetVec)
  if (any(is.na(response))) 
    stop("missing values in response variable not allowed")
  resp <- initVariableFrame(response, trafo = ytrafo, response = TRUE, 
                            scores = scores)
  RET <- new("LearningSampleFormula", inputs = inp, responses = resp, 
             weights = rep(1, inp@nobs), nobs = inp@nobs, ninputs = inp@ninputs, 
             menv = dat)
  RET
}


# altered version of cforest()
myCforest <- function (predMat, targetVec, subset = NULL, weights = NULL, 
                       controls = cforest_unbiased(), xtrafo = ptrafo, ytrafo = ptrafo, 
                       scores = NULL) 
{
  ls <- dpp(myRF, predMat, targetVec, subset, xtrafo = xtrafo, 
            ytrafo = ytrafo, scores = scores)
  fitmem <- ctree_memory(ls, TRUE)
  fit(myRF, ls, controls = controls, weights = weights, 
      fitmem = fitmem)
}


sample.from.classes <- function(n, sample.class, split.frac=0.8, replace=F, balanced=T) {
  tp <- which(sample.class)
  fn <- which(!sample.class)
  tp.n <- as.integer(split.frac * length(tp))
  if(length(sample.class) - tp.n <= 1) {
    # if our fraction would have us holding out one or less positive examples per tree, then
    # explictly hold one out at a time
    tp.to.exclude <- sample(tp, 1)
    tp.to.use <- setdiff(tp, tp.to.exclude) 
  } else {
    tp.to.use <- sample(tp, tp.n)
  }

  fn.to.use <- if(length(fn) == 1) {
    rep(fn, n-tp.n)
  } else {
    sample(fn, n-tp.n, replace=replace)
  }

  v <- rep(0, length(sample.class))
  if(balanced) {
    for(i in tp.to.use) {
      v[i] <- v[i] + length(fn.to.use)
    }
    
    for(i in fn.to.use) {
      v[i] <- v[i] + length(tp.to.use)
    }
  } else {
    for(i in tp.to.use) {
      v[i] <- v[i] + 1
    }
    
    for(i in fn.to.use) {
      v[i] <- v[i] + 1
    }
  }

  v
}

sample.weights <- function(sample.count, total.frac, split.frac, classes, replace=F, balanced=T) {
  n <- as.integer(length(classes)*total.frac)
  w <- foreach(i=seq(sample.count), .combine=rbind) %do% {
    sample.from.classes(n, classes, replace=replace, balanced=balanced)
  }
  t(w)
}

create.balanced.weights <- function(pos.samples, max.weight=NULL, mask=rep(T, length(pos.samples))) {
    stopifnot(is.logical(pos.samples))
    stopifnot(is.logical(mask))

    pos.count <- sum(pos.samples & mask)
    neg.count <- sum(!pos.samples & mask)

    if(pos.count == 0 || neg.count == 0) {
      return(ifelse(mask, 1, 0))
    }

    if(is.null(max.weight)) {
      if(pos.count < 10) {
          pos.count <- 10
      }
    } else {
      a <- (1 - max.weight * pos.count)
      if(a > 0) {
          max.neg.count <- round(max.weight * neg.count / a * pos.count)
          if(neg.count > max.neg.count) { 
              neg.count <- as.integer(max.neg.count)
          } 
      }
    }
    
    stopifnot(neg.count > 0)

    v <- ifelse(pos.samples, neg.count, pos.count)
    v[!mask] <- 0
    v
}

alt.biased.sampling.weights <- function(n, pos.samples, prob=0.8, max.weight=NULL) {
    stopifnot(is.logical(pos.samples))
    selected.pos <- which(runif(length(pos.samples)) < prob & pos.samples)
    if(n - length(selected.pos) < 0){
      selected.neg <- sample(which(!pos.samples), 1, replace=F)
    } else if(n - length(selected.pos) > length(which(!pos.samples))){
      selected.neg <- sample(which(!pos.samples), length(which(!pos.samples)), replace=F)
    } else {
      selected.neg <- sample(which(!pos.samples), n-length(selected.pos), replace=F)
    }
    
    mask <- rep(F, length(pos.samples))
    mask[c(selected.pos, selected.neg)] <- T
    
    create.balanced.weights(pos.samples, mask=mask, max.weight=max.weight)
    }

sample.weights.w.no.min <- function(ntrees, sample.fraction, prob, pos.samples, max.weight=NULL) {

    foreach(i=seq(ntrees), .combine=cbind) %do% {
        alt.biased.sampling.weights(round(length(pos.samples) * sample.fraction), pos.samples, prob=prob, max.weight=max.weight)
    }
}

has.sufficent.sensitive.lines <- function(targetVec, fitControlSettings) {
  if(length(grep("biased", fitControlSettings) > 0)) {
    if(grepl("viability", fitControlSettings)){
      sum(is.sensitive(targetVec, viability = TRUE), na.rm=T) >= 2 
    } else if(grepl("viability-auc", fitControlSettings)){
      sum(is.sensitive(targetVec, viability = TRUE, auc = TRUE), na.rm=T) >= 2 
    } else {
      sum(is.sensitive(targetVec), na.rm=T) >= 2 
    }
  } else {
    TRUE
  }
}

is.sensitive <- function(targetVec, viability=FALSE, auc=FALSE) {
  if(viability){
    if(auc){
      return(targetVec < 0.6)
    } else {
      return(targetVec < 40)
    }
  } else {
  	return(targetVec < -2)
  }
}

# builds a cforest (as defined in the party package)
# targetMat - a matrix with one column of the target vector
# preadMat - predictor matrix with predictors in columns
#   targetMat and predMat must have identical rownames (samples)
# ntree - number of trees to build
buildCForest <- function(targetMat, predMat, ntree, fitControlSettings) {
  
  cat("Building a cforest\n")

  compress.fn <- identity
  expand.fn <- identity
 
  if(ncol(predMat) > 5000) {
    cat("Lots of features (",ncol(predMat),") enabling tree compression\n");
    compress.fn <- compress.via.pack
    expand.fn <- expand.via.unpack
  }

  case.weights <- NULL
  per.sample.weights <- NULL
  if (fitControlSettings == "default") {
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
                      trace=T,
                      compress=compress.fn,
                      expand=expand.fn)
  } else if(fitControlSettings == "biased-1-deep") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T, stump=T,
                      compress=compress.fn,
                      expand=expand.fn)
  } else if(fitControlSettings == "biased-varOnce") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T, varOnce=T,
                      compress=compress.fn,
                      expand=expand.fn)
  } else if(fitControlSettings == "biased-varOnce-z") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T, varOnce=T,
                      compress=compress.fn,
                      expand=expand.fn)
    targetMat <- ifelse(lt.neg.2, targetMat, 0)
  } else if(fitControlSettings == "biased-2-deep") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T, maxdepth=2,
                      compress=compress.fn,
                      expand=expand.fn)
  } else if(fitControlSettings == "biased-1-deep-z") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T, stump=T,
                      compress=compress.fn,
                      expand=expand.fn)
    targetMat <- ifelse(lt.neg.2, targetMat, 0)
  } else if(fitControlSettings == "biased-2-deep-z") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T, maxdepth=2,
                      compress=compress.fn,
                      expand=expand.fn)
    targetMat <- ifelse(lt.neg.2, targetMat, 0)
  } else if(fitControlSettings == "biased-sampling-wmin4") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T,
                      compress=compress.fn,
                      expand=expand.fn)
  } else if(fitControlSettings == "biased-sampling-wmin4-viability") {
    lt.neg.2 <- is.sensitive(targetMat, viability=TRUE)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T,
                      compress=compress.fn,
                      expand=expand.fn)
  } else if(fitControlSettings == "biased-sampling-wmin4-viability-auc") {
    lt.neg.2 <- is.sensitive(targetMat, viability=TRUE, auc=TRUE)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
      ntree=ntree, 
      mtry=0, # equiv. to ncol(predMat)
      savesplitstats=F,
      minsplit = 2*3*max.weight, minbucket = 3*max.weight,
      trace=T,
      compress=compress.fn,
      expand=expand.fn)
    
  } else if(fitControlSettings == "biased-zero-wmin4") {
    lt.neg.2 <- is.sensitive(targetMat)
    case.weights <- sample.weights.w.no.min(ntree, 0.67, 0.8, lt.neg.2, max.weight=0.05)
    per.sample.weights <- create.balanced.weights(lt.neg.2, max.weight=0.05)
    max.weight <- max(case.weights)
    controls <- cforest_unbiased(
                      ntree=ntree, 
                      mtry=0, # equiv. to ncol(predMat)
                      savesplitstats=F,
			minsplit = 2*3*max.weight, minbucket = 3*max.weight,
                      trace=T,
                      compress=compress.fn,
                      expand=expand.fn)
    targetMat <- ifelse(lt.neg.2, targetMat, 0)
  } else {
    stop(paste("Unknown value for fitControlSettings:", fitControlSettings))
  }

  mycf <- myCforest(predMat = predMat, 
                    weights = case.weights,
                    targetVec = targetMat, 
                    controls = controls)
  
      # predict relies on access to "input" in the menv, however, that only exists for formulaModelEnvs or something -- so patch it to use get designMatrix instead which
      # appears reasonable for MatrixModelEnvs
      orig.get <- environment(mycf@predict_response)$object@menv@get
      new.get <- function(which, data = NULL, ...) { if(which == "input") { orig.get("designMatrix", data=data, ...) } else { orig.get(which, data=data, ...) } }
      environment(mycf@predict_response)$object@menv@get <- new.get
      mycf@data@get <- new.get

  OOBinfo <- getOOBInfo(mycf, per.sample.weights)
  cat(paste("OOB", OOBinfo$metric, "is:", signif(OOBinfo$quality, 3), "\n"))
  
  list(model=mycf, OOBinfo=OOBinfo)
}

collapse.weights <- function(weights) { 
  if (is.null(weights)) {
    return(NULL)
  }

  per.sample <- apply(weights, 1, mean)
  for(i in seq(ncol(weights))) {
    stopifnot( all(weights[i, ] == per.sample | weights[i, ] == 0) )
  }
  
  per.sample
}

# returns the relative predictor importance for a given myCForest.
# If quantileCutOff is not NULL (and between 0 and 1) it includes only predictors with 
# importance values greater than the absolute
# value of the minimal (negative) importance of any predictor.
getVarImp <- function(mycf, quantileCutOff=NULL, plotHist=FALSE) {
  cat("Computing predictor importance...")

  if (mycf@responses@is_nominal) {
    
    num_permutations <- 10
    num_features <- ncol(mycf@data@get("designMatrix"))
    
    keep <- rep(TRUE,num_features)
    vi_sums <- rep(0,num_features)
    
    for (i in 1:num_permutations){
      vi <- varimpAUC(mycf, nperm=1)
      vi_sums <- vi_sums + vi
      vi_positive <- vi > 0
      keep <- keep & vi_positive
    }
    
    vi_avg <- vi_sums / num_permutations
    vi_avg <- vi_avg[keep]
    return(vi_avg[order(vi_avg, decreasing=T)] / sum(vi_avg))
    
  } else {
    vi <- varimp(mycf, nperm=10)
  }

  cat("Done.\n")
  if (plotHist) {
    hist(vi[vi!=0], 20)
  }
  
  cat(paste("There are", sum(vi != 0), "predictors with non-zero importance."))
  cat(paste("Of them,", sum(vi<0), "have negative importance.\n"))
  
  if ((!is.null(quantileCutOff)) && (sum(vi<=0) > 0)) { # if there are predictors with negative (or zero) importance, use them to set cut-off for positive importance
    cutoff <- abs(quantile(vi[vi<=0], quantileCutOff))
    cat(paste("Selecting only variables with importance >", cutoff, "\n."))
    vi <- vi[vi > cutoff]
  }
  cat(paste("Returning importance for", 
            length(vi), 
            "predictors.\n"))
  
  vi[order(vi, decreasing=T)] / sum(vi)
}

get.tree.depths <- function(tree, depth=1) {
  if(tree[[4]]) {
    depth
  } else {
    c(get.tree.depths(tree[[8]], depth+1), get.tree.depths(tree[[9]], depth+1))
  }
}

get.ensemble.depths <- function(model) { 
  table(do.call(c, lapply(model@ensemble, function(x) { get.tree.depths(model@expand(x)) })))
}


# fits a cforest model that predicts targetMat[,targetID] using predMat
# features in columns (for bost targetMat and predMat)
# targetMat and predMat have to have the same rownames
# targetID specifieis which column of targetMat to use as the targetVec
# returns the cforest model
fitOptimizedCForest <- function(targetMat, predMat, targetID=1, ntrees=200, showPlots=FALSE, fitControlSettings = "default") {
  stopifnot(all(rownames(targetMat)==rownames(predMat)))
  
#  .Call("R_print_pack_alloc_summary")
#  print(gc());
  
  targetVec <- targetMat[, targetID, drop=F]
  res <- buildCForest(targetVec, predMat, ntrees, fitControlSettings)
  mycf <- res$model
  OOBinfo <- res$OOBinfo
  
  if(showPlots) {
    if (mycf@responses@is_nominal) {
      plot(sapply(treeresponse(mycf), I)[2,], targetVec)
    } else {
      plot(predict(mycf, OOB=T), targetVec)
    }
  }
  
  vi  <- getVarImp(mycf, 0.01, plotHist = showPlots)
  fitInfo <- extractFitInfo(mycf, OOBinfo, vi)

  print(vi)
  print(get.ensemble.depths(mycf))

#  .Call("R_print_pack_alloc_summary")
#  print(gc());  

  if( !is.na(OOBinfo$weighted.cor) && OOBinfo$weighted.cor < 0 ) {
     fitInfo$OOB$quality <- 0

     return (list(model="inv-corr", fitInfo=fitInfo))
  }
  
  if(length(vi) == 0) {
    return (list(model="no-vars", fitInfo=fitInfo))
  }

  res <- buildCForest(targetVec, predMat[, names(vi), drop=F], 500, fitControlSettings)
  mycf1 <- res$model
    
  if(showPlots) {
    if (mycf1@responses@is_nominal) {
      plot(sapply(treeresponse(mycf1), I)[2,], targetVec)
    } else {
      plot(predict(mycf1, OOB=T), targetVec)
    }
  }
    
  system.time( vi1  <- getVarImp(mycf1, quantileCutOff = 0.0, plotHist = showPlots) )
  print(vi1)
  print(get.ensemble.depths(mycf1))
  rm(mycf1)
  
  fitInfo$varImp <- vi1
  if(length(vi1) == 0) {
    return (list(model="no-vars", fitInfo=fitInfo))
  } 

  res <- buildCForest(targetVec, predMat[, names(vi1), drop=F], 500, fitControlSettings)
  mycf2 <- res$model

  if(showPlots) {
    if (mycf2@responses@is_nominal) {
      plot(sapply(treeresponse(mycf2), I)[2,], targetVec)
    } else {
      plot(predict(mycf2, OOB=T), targetVec)
    }
  }

  system.time( vi2  <- getVarImp(mycf2, quantileCutOff = 0.0, plotHist = showPlots) )
  print(vi2)
  print(get.ensemble.depths(mycf2))
  
  if(length(vi2) > 30) {
    vi2 <- head(vi2, 30)
  }

  fitInfo$varImp <- vi2
  if(length(vi2) == 0) {
    return (list(model="no-vars", fitInfo=fitInfo))
  }

  return (list(model=mycf2, fitInfo=fitInfo, first.model=mycf))
}

# When OOBinfo is passed we will use it instead of computing it here. This is used to pass
# quality metrics from the first model built in buildCForest() as later models' are contaminated' R2 is
# contaminated (due to feature selection using all samples)
extractFitInfo <- function(model, OOBinfo, vi) {
  
  i <- list()
  i$modelType <- ifelse(model@responses@is_nominal, "Classification", "Regression")
  i$nobs <- model@responses@nobs
  i$targetVec <- model@responses@variables[[1]] # save as vector
  i$sampleNames <- rownames(model@responses@variables)
  i$OOB <- OOBinfo
  i$varImp <- vi
  i$predictors <- names(i$varImp)
  i$predData <- model@data@get("designMatrix")[, i$predictors, drop=F]
  i$ntrees <- length(model@ensemble)
  i$report.summary <- c(sprintf("# of samples: %d", i$nobs), sprintf("nTrees: %d", i$ntrees), sprintf("selectedFeatures: %d", length(i$predictors)))
  i
}

# returns OOB classification predictions
predictOOBProbs <- function(mycf, OOB=T) {
  p <- treeresponse(mycf)
  preds <- sapply(p, I)[2,]
  
  preds
}

weighted.cor <- function(x, y, w){
  stopifnot(length(w) == length(x))
  w <- w/sum(w, na.rm=T)
  w.mean.x <- sum(x*w, na.rm=T)/sum(w, na.rm=T)
  w.mean.y <- sum(y*w, na.rm=T)/sum(w, na.rm=T)
  v <- sum(w) - (sum(w**2, na.rm=T)/sum(w, na.rm=T))
  w.cov <- function(x,y,w) { sum(w*(x-w.mean.x)*(y-w.mean.y), na.rm=T)/v }
  cor.r <- w.cov(x,y,w)/((w.cov(x,x,w)*w.cov(y,y,w))**0.5)

  cor.r
}

# returns prediction, quality, metric used
# handles both regression and classification cases
# used by extractFitInfo()
getOOBInfo <- function(model, weights=NULL) {
  OOB <- list()
  is_classification <- model@responses@is_nominal
  targetVec <- model@responses@variables[[1]] #get("responseMatrix", env=model@data@env)[,1]
  OOB$nfeatures <- ncol(model@data@get("designMatrix"))
  
  if (is_classification) {
    stopifnot(is.null(weights))
    library(AUC) # auc, roc
    OOB$prediction <- predictOOBProbs(model)
    OOB$quality <- auc(roc(OOB$prediction, targetVec))
    OOB$metric <- "AUC"
    OOB$weights <- weights
  } else {
    OOB$prediction <- as.numeric(predict(model, OOB=T))
    if(is.null(weights)) {
       weights <- rep(1, length(targetVec))
    }
    OOB$metric <- "cor R2"
    OOB$quality <- weighted.cor(OOB$prediction, targetVec, rep(1, length(targetVec))) ** 2
    OOB$weights <- weights
    OOB$weighted.cor <- weighted.cor(OOB$prediction, targetVec, weights)
    OOB$weighted.cor.R2 <- OOB$weighted.cor * OOB$weighted.cor
    if(!all(weights == 1)) {
      OOB$metric <- "wcor R2"
      OOB$quality <- OOB$weighted.cor.R2
    }
  }
  
  OOB
}

# generates a data frame with cell lines in rows, where 
# first column is the target feature,
# second is predicted target feature,
# the rest contain all the predictive features used in model.
# Params:
# fit.info - a fit.info object
# targetID - target feature ID 
generatePredictorsTable <- function(fit.info, targetID) {
  fit.table <- data.frame(target=fit.info$targetVec, 
                          pred_target=fit.info$OOB$prediction, 
                          fit.info$predData[,names(fit.info$varImp),drop=F], 
                          row.names=fit.info$sampleNames, 
                          check.names=F)
  colnames(fit.table)[1] <- targetID
  colnames(fit.table)[2] <- paste(targetID, 'pred', sep='_')
  
  fit.table
}
