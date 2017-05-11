library(gridExtra)
library(foreach)
library(atlantis)

createBulkAtlantisFlock <- function(analysis.description, run_mode, predMat.file, targetMat.file, anno.file, featureSetToUse, kMinNumSamples=50, seed=300, batch.size=1, predFeaturesPerTarget=NULL, fitControlSettings="default", config.summary=NULL, predictOnly=NULL, permutation.count=NULL, only.permutations=F, 
  code.dir=".",
  test.fraction=NA,
  predFeatureRegExps=NULL) {
  # analysis.description: a description which will be saved in the run directory
  # run_mode: ignored
  # mustHave: a vector of dataset ids which must have a value present to be included in the prediction matrix
  # niceToHave: a vector of dataset ids to include in the prediction matrix, but these may have missing values in the prediction matrix
  # featureSetToUse: one of "all", "perTarget", "perTargetAndSelf" or "single" determining how to subset the prediction matrix when
  #                  for a given target.
  # predFeaturesPerTarget: if featureSetToUse = "perTarget" or "perTargetAndSelf" then this should be a dataframe with two columns: c("target", "partner") which will be used to look up which genes to use as predicitive features based on a given target

  if (featureSetToUse == "perTarget" || featureSetToUse == "perTargetAndSelf") {
	  stopifnot(!is.null(predFeaturesPerTarget) && is.data.frame(predFeaturesPerTarget) && all(colnames(predFeaturesPerTarget) == c("target", "partner")))
  } else if(featureSetToUse == "perGroup") {
     stopifnot(!is.null(predFeaturesPerTarget) && is.data.frame(predFeaturesPerTarget) && all(colnames(predFeaturesPerTarget) == c("group", "target", "group.name")))
  }

  # flock_run_dir is the name of the directory of where the results will be written
  # the analysis.ID is the trailing dir name
  analysis.ID <- basename(flock_run_dir)
  output.dir <- paste(flock_run_dir, '/results/', sep='')
  dir.create(output.dir)

  # cache the functional related genes if we'll need them
  predFeaturesPerTarget.file <- NULL
  if(featureSetToUse == "perTarget" || featureSetToUse == "perTargetAndSelf" || featureSetToUse == "perGroup") {
    stopifnot(!is.null(predFeaturesPerTarget))

    predFeaturesPerTarget.file <- paste(output.dir, '/featuresPerTarget.Rdata', sep="")
    save(predFeaturesPerTarget, file=predFeaturesPerTarget.file)
  }

  # load the target matrix to determine the targets
  load(targetMat.file);
  stopifnot(is.matrix(targetMat) || is.data.frame(targetMat))
  targets <- colnames(targetMat);
  if(is.null(permutation.count)) {
    permutation.count <- length(targets)
  }

  if(featureSetToUse == "perGroup") {
    # then we want to create a target per member in our predFeaturesPerTarget table
    with.valid.targets <- predFeaturesPerTarget[paste("GS_",predFeaturesPerTarget$target,sep='') %in% targets,,drop=F]
    targets <- lapply(seq(nrow( with.valid.targets)), function(i) { 
	x<-with.valid.targets[i,]; 
        list(ID=paste("GS_",as.character(x$target),sep=''), 
             gene=as.character(x$target), 
             group=as.character(x$group), 
             output.prefix=paste(x$group,'-',x$target,sep=''),
             group.name=as.character(x$group.name)) 
        } )
    stopifnot(! ("GS_NA" %in% sapply(targets, function(t){t$ID})))
  } else {
    # filter out targets which have a single output value
    print("before filtering out singles")
    print(length(targets))
    is.bad.target <- sapply(seq(dim(targetMat)[[2]]), function(i) { t<-targetMat[,i]; t<-t[!is.na(t)]; length(unique(t))==1 } )
    targets <- targets[!is.bad.target]
    targets <- lapply(targets, function(t) { list(ID=t, output.prefix=t) })
    print("after filtering out targets with a single target value")
    print(length(targets))
    stopifnot(grep("GS_NA", targets)!=0)
  }

  common <- list( 
    fitControlSettings=fitControlSettings,
    targetMat.file=targetMat.file,
    predMat.file=predMat.file,
    anno.file=anno.file,
    output.dir=output.dir,
    featureSetToUse=featureSetToUse,
    kMinNumSamples=kMinNumSamples,
    predFeaturesPerTarget.file=predFeaturesPerTarget.file,
    seed=seed,
    quality.dist.file=paste(output.dir, "null_quality_distribution.txt", sep=""),
    config.summary=config.summary,
    predictOnly=predictOnly,
    permutation.count=permutation.count,
    only.permutations=only.permutations,
    batch.size=batch.size,
    predFeatureRegExps=predFeatureRegExps
    )

  sources <- c("createBulkAtlantisFlock.R")

  common$sources <- sources

  # add in runs for each permuted version used for the null
  permutations <- sample(targets, permutation.count, replace=T)
  for(i in seq_along(permutations)) { 
    permutations[[i]]$permute.seed <- i + seed
  }

  if(common$only.permutations) {
    tasks <- permutations
  } else {
    tasks <- c(targets, permutations)
  }

  tasks <- split(tasks, as.integer(seq(length(tasks))/common$batch.size))  

  if(!is.na(test.fraction)) {
    task.count <- as.integer(ceiling(test.fraction * length(tasks)))
    cat("Only running", task.count, "tasks out of", length(tasks), "because test.fraction was set to", test.fraction, "\n")
    tasks <- tasks [ seq(task.count) ]
  }

  flock.run(tasks, task_function='per_target', gather_function='gather', sources=common$sources, flock_common_state=common)
}

prev.proc.time <- NULL

time.checkpoint <- function(label) {
#  now <- proc.time()
#  if(!is.null(prev.proc.time)) {
#    cat(sprintf("%s: %f\n", label, (now["elapsed"] - prev.proc.time["elapsed"])))
#  } else {
#    cat(sprintf("%s: first\n", label))
#  }
#  prev.proc.time <<- now
}

# given a gene symbol, returns a list of functionaly related genes based on dataframe which was persisted
get.related.genes <- function(predFeaturesPerTarget, targetToLookup, target, include.self) {
  matches <- as.character(predFeaturesPerTarget$target) == targetToLookup
  partners <- as.character(predFeaturesPerTarget[matches, "partner"])
  if(include.self) {
    ret <- unique(c(target, partners))
  } else {
    ret <- setdiff(unique(partners), target)
  }

  return(ret)
#  sapply(ret, getFeatureNameFromFeatureID)
}

get.gene.group <- function(predFeaturesPerTarget, groupToLookup, target, include.self) {
  partners <- as.character(predFeaturesPerTarget[predFeaturesPerTarget$group == groupToLookup, "target"])
  if(include.self) {
    ret <- unique(c(target, partners))
  } else {
    ret <- setdiff(unique(partners), target)
  }
  return(ret)
#  sapply(ret, getFeatureNameFromFeatureID)
}

determine.feature.genes <- function(target, featureSetToUse, predFeaturesPerTarget.file) {
  if(!is.null(predFeaturesPerTarget.file)) {
    load(predFeaturesPerTarget.file)
  }

  if (featureSetToUse == "all") {
    genes <- NULL
  } else if (featureSetToUse == "perTarget" || featureSetToUse == "perTargetAndSelf") {
    targetName <- getFeatureNameFromFeatureID(target$ID)
    stopifnot(!is.null(targetName))
    genes <- get.related.genes(predFeaturesPerTarget, targetName, targetName, featureSetToUse == "perTargetAndSelf")
    cat("found", length(genes), "gene names to be used as features for", targetName, "\n")
  } else if (featureSetToUse == "perGroup") {
    stopifnot(all(colnames(predFeaturesPerTarget) == c("group", "target", "group.name")))
    group <- target$group
    stopifnot(! is.null(group))
    stopifnot(! is.null(target$gene))
    genes <- get.gene.group(predFeaturesPerTarget, group, target$gene, FALSE)
  } else {
    stopifnot(featureSetToUse ==  "single")
    targetName <- getFeatureNameFromFeatureID(target$ID)
    genes <- targetName
    stopifnot(!is.null(genes))
  }

  return(genes)
}

per_target <- function(per.task.state, common.state, output.file, run) {
  library(foreach)

  flock_run_dir <- run$dir

  analysis.dir <- basename(flock_run_dir);
  analysis.ID <- analysis.dir;
  output.dir <- flock_common_state$output.dir
  predMat.file <- flock_common_state$predMat.file
  targetMat.file <- flock_common_state$targetMat.file
  anno.file <- flock_common_state$anno
  kMinNumSamples <- flock_common_state$kMinNumSamples
  predFeaturesPerTarget.file <- flock_common_state$predFeaturesPerTarget.file
  config.summary <- flock_common_state$config.summary
  predFeatureRegExps <- flock_common_state$predFeatureRegExps

  fitControlSettings = flock_common_state$fitControlSettings

  # the following figures out what are the 100 most similar genes to target gene and uses only them
  # for prediction.
  ####
  targetMat.file <- flock_common_state$targetMat.file
  load(targetMat.file)

  featureSetToUse <- flock_common_state$featureSetToUse

  fit.single.target <- function(target) {
    # allow target to be specified as either a string or a list of (ID, group)
    if(is.character(target)) {
      target <- list(ID=target, targetUsedToLookupRelated=target, permute.seed=NULL)
    }

    targetID <- target$ID
    permute.seed <- target$permute.seed

    print(targetID)
    cat(sprintf("Generating model for %s\n", targetID))
    set.seed(flock_common_state$seed)

    if (is.numeric(targetID)) {
      targetID <- colnames(targetMat)[targetID]
    }

    stopifnot(length(grep("NO_CURRENT", targetID)) == 0)

    permuteRows <- function(seed, targetMat) {
      saved <- .Random.seed
      set.seed(seed)
      permutedTargetMat <- targetMat
      for(i in seq(ncol(targetMat))) { 
        permutedTargetMat[,i] = sample(targetMat[,i]) 
      }
      .Random.seed <- saved
      permutedTargetMat
    }

    if(featureSetToUse == "all") {
      predFeatureRegExps <- ".*"
      genes <- NULL
    } else {
      genes <- determine.feature.genes(target, featureSetToUse, predFeaturesPerTarget.file)
    }

    output.prefix <- target$output.prefix
    group.name <- target$group.name
    if(is.null(group.name)) {
      group.name <- NA
    } else {
      config.summary <- c(config.summary, paste("group:", group.name))
    }

    # if this is for the null distribution
    if(!is.null(permute.seed)) {
      permutedTargetMat <- permuteRows(permute.seed, targetMat)

      res <- runATLANTIS(
        analysis.ID=analysis.ID,
        targetID=targetID,
        output.prefix=paste(output.prefix, "-NULL-",permute.seed, sep=''),
        makePlot=FALSE,
        output.dir = paste(dirname(output.dir), '/temp', sep=''),
        fitControlSettings = fitControlSettings,
        additionalFeatures=NULL,
        predFeatureNamesToUse=genes,
        predFeatureRegExps=predFeatureRegExps,
        predMat.file = predMat.file,
        targetMat = permutedTargetMat,
        anno.file=anno.file,
        kMinNumSamples=kMinNumSamples,
        save.params=F,
        save.model=F,
        save.featureData=F)
    } else {
      res <- runATLANTIS(
        analysis.ID=analysis.ID,
        targetID=targetID, 
        makePlot=TRUE,
        output.dir = output.dir,
        output.prefix=output.prefix,
        fitControlSettings = fitControlSettings,
        additionalFeatures=NULL,
        predFeatureNamesToUse=genes,
        predFeatureRegExps=predFeatureRegExps,
        predMat.file = predMat.file,
        targetMat = targetMat,
        anno.file=anno.file,
        kMinNumSamples=kMinNumSamples,
        save.params=F,
        save.model=F,
        report.summary=config.summary,
        predictOnly=flock_common_state$predictOnly
        )
    }

    list(targetID=targetID,
         result=res,
         permute.seed=permute.seed,
         output.prefix=output.prefix,
         group.name=group.name)

  }

  task.results <- lapply(flock_per_task_state, fit.single.target)

  save(task.results, file=output.file)
}

    zscore <- function(x) {
     (x-mean(x, na.rm=T))/sd(x, na.rm=T)
    }


gather <- function(per.task.state, common.state, output.file, run) {
  task.results <- do.call(c, lapply(flock_per_task_state, function(job.details) { 
    e <- new.env()
    load(job.details$flock_output_file, envir=e)
    e$task.results
  }))

  is.permuted <- sapply(task.results, function(x) { !is.null(x$permute.seed) } )

  null.dist.details <- foreach(task.result=task.results[is.permuted], .combine=rbind) %do% {
    load(task.result$result$fit.file)
    if(!is.null(fit.info$OOB$weighted.cor)) {
      data.frame(quality=fit.info$OOB$quality, weighted.cor=fit.info$OOB$weighted.cor, weighted.cor.R2=fit.info$OOB$weighted.cor.R2, with.more.weight=sum(min(fit.info$OOB$weights) != fit.info$OOB$weights), sample.count=length(fit.info$targetVec), failure.reason=NA, targetID=task.result$targetID, nfeatures=fit.info$OOB$nfeatures)
    } else {
      data.frame(quality=NA, weighted.cor=NA, weighted.cor.R2=NA, with.more.weight=NA, sample.count=NA, failure.reason=fit.info$failure.reason, targetID=task.result$targetID, nfeatures=NA)
    }
  }

  null.distribution.values <- na.omit(null.dist.details$quality[!is.na(null.dist.details$quality)])
  pval.from.null.distribution <- function(x) { (sum(null.distribution.values >= x)+1)/(length(null.distribution.values)+1) }

  write.table(null.distribution.values, file=flock_common_state$quality.dist.file, col.names=F, row.names=F)
  write.csv(null.dist.details, file=paste(flock_common_state$quality.dist.file,"-details.csv",sep=''),row.names=F)

  # summarize the runs into a table 

  fits <- lapply(task.results[!is.permuted], function(job) { 
    if(!is.null(job$result) && !is.null(job$result$fit.file)) {
        load(job$result$fit.file)
        # reduce memory by dropping unused fields
        top.feature <- names(which.max(fit.info$varImp))
        if(!is.null(top.feature)) {
          top.feature.values <- fit.info$predData[,top.feature]
          wt.cor.top.feature <- atlantis:::weighted.cor(top.feature.values, fit.info$targetVec, fit.info$OOB$weights)
          target.top.pred.corr <- cor(fit.info$targetVec, top.feature.values, use='pairwise.complete')
          lines.with.more.weight <- which(min(fit.info$OOB$weights) < fit.info$OOB$weights)
          sensitive.lines.zmean.top.feature <- mean(zscore(top.feature.values)[lines.with.more.weight], na.rm=T)
        } else {
          target.top.pred.corr <- NA
          wt.cor.top.feature <- NA
          sensitive.lines.zmean.top.feature <- NA
        }
        fit.info$top.feature <- top.feature
        fit.info$target.top.pred.corr <- target.top.pred.corr
        fit.info$wt.cor.top.feature <- wt.cor.top.feature
        fit.info$sensitive.lines.zmean.top.feature <- sensitive.lines.zmean.top.feature
        fit.info$predData <- NULL
        fit.info$predictors <- NULL
        fit.info
    } else {
      NULL
    }
  })

  taskResults <- task.results[!is.permuted]

  varsPerModel <- foreach(taskResult=taskResults, fit.info=fits, .combine=rbind) %do% {
    if(length(fit.info$varImp) == 0) {
      NULL
    } else {
      d <- data.frame(variable=names(fit.info$varImp), varImportance=fit.info$varImp, targetID=taskResult$targetID, group.name=taskResult$group.name)
      rownames(d) <- NULL
      d
    }
  }

  coerce.to.na <- function(x) { ifelse(is.null(x), NA, x) }

  qualPerModel <- foreach(taskResult=taskResults, fit.info=fits, .combine=rbind) %do% {
    targetID <- taskResult$targetID
    Rsquared <- coerce.to.na(fit.info$OOB$quality)
    if(!is.null(fit.info$failure.reason)) {
      Rsquared <- NA
    }
    target.top.pred.corr <- fit.info$target.top.pred.corr
    top.feature <- fit.info$top.feature
    wt.Rsquared <- fit.info$OOB$weighted.cor.R2
    Pvalue <- pval.from.null.distribution(Rsquared)
    n <- names(fit.info$varImp[top.feature])
    group.name <- taskResult$group.name
    if(is.null(group.name)) {
      group.name <- NA
    }
    output.prefix <- taskResult$output.prefix
    target.min <- min(fit.info$targetVec, na.rm=T)
    sample.count <- length(fit.info$targetVec)
    with.more.weight <- sum(min(fit.info$OOB$weights) != fit.info$OOB$weights)
    nfeatures=fit.info$OOB$nfeatures
    if(is.null(nfeatures)) { 
      nfeatures <- NA
    }
    print("nfeatures")
    print(nfeatures)

    #stopifnot(!is.null(Rsquared))
    #stopifnot(!is.null(targetID))

    data.frame(Rsquared=Rsquared,
               targetID=targetID, 
               topFeature=coerce.to.na(n), 
               Pvalue=coerce.to.na(Pvalue), 
               topFeatureCor=target.top.pred.corr, 
               targetMin=target.min,
               group.name=group.name,
               output.prefix=output.prefix,
               wt.Rsquared=coerce.to.na(wt.Rsquared),
               sample.count=sample.count,
               with.more.weight=with.more.weight,
               failure.reason=coerce.to.na(fit.info$failure.reason),
               nfeatures=nfeatures,
               sensitive.lines.zmean.top.feature=fit.info$sensitive.lines.zmean.top.feature,
               wt.cor.top.feature=fit.info$wt.cor.top.feature)
  }

  qualPerModel$percentile <- rank(qualPerModel$Rsquared,na.last=F)/nrow(qualPerModel)
  qualPerModel$fdr <- p.adjust(qualPerModel$Pvalue, method="fdr")

  cat("writing summary\n")
  write.csv(varsPerModel, file=paste(flock_run_dir,'/results/varsPerModel.csv',sep=''), row.names=F)
  write.csv(qualPerModel, file=paste(flock_run_dir,'/results/qualPerModel.csv',sep=''), row.names=F)

  annTable.file <- flock_common_state$anno.file

  # make report summary plot with:
  # 1. BRAF, KRAS, CTNNB1, ESR1
  # 2. random 4 out of the top 20 models of that bulk run in terms of R^2
  # 3. random 4 models with p-value > 0.1 

  if(length(qualPerModel) > 4) {

    sample.at.most <- function(v, n) {
       sample(v, min(length(v), n))
    }

    randTop4 <- sample.at.most(qualPerModel$targetID[rank(-qualPerModel$Rsquared) <= 20], 4)
    randBad4 <- sample.at.most(qualPerModel$targetID[qualPerModel$Pvalue > 0.1], 4)

    add.atlantis.plots <- function(targetIds) {
      for(targetID in targetIds) {
        e <- new.env()
        fi <- sprintf("%s/results/%s_fit.info_cforest.Rdata", flock_run_dir, targetID)
        sum.file <- sprintf("%s/results/%s_ATLANTIS_Summary.Rdata", flock_run_dir, targetID)
        if(file.exists(fi) && file.exists(sum.file)) {
          load(fi, envir=e)
          PlotATLANTISresults(e$fit.info, sum.file,
            NA, NULL, annTable.file)
        } else {
          # write something to file to indicate target not present
        }
      }
    }

    pdf(file.path(flock_run_dir, 'results/summary.pdf'), w=11*1.5, h=8.5*1.5)
    pardefault <- par()
    grid.table(head(qualPerModel[order(qualPerModel$Rsquared, decreasing=T),],20))
    par(pardefault)
    add.atlantis.plots(c("GS_BRAF", "GS_KRAS", "GS_CTNNB1", "GS_ESR1"))
    par(pardefault)
    add.atlantis.plots(randTop4)
    par(pardefault)
    add.atlantis.plots(randBad4)
    par(pardefault)
    if(length(null.distribution.values) > 0) {
      hist(null.distribution.values)
    }
    par(pardefault)
    dev.off()
  }
}

# given a feature ID (a compound string that includes a feature type and a feature name),
# returns the feature name
getFeatureNameFromFeatureID <- function(ID) {
  res <- strsplit(ID, "[_:]")[[1]]
  if (length(res) < 2)
    NULL
  else
    res[[2]]
}
