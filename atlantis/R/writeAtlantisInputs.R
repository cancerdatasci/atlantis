# Construct a matrix and an annotation table for a set of matrices.
# "sources" is a list of lists, where each list has: dat, prefix, vartype, must.have.cellline
constructMergedMatrix <- function(sources, loci.by.symbol, limitToCellLines=NULL) {
  stopifnot(sapply(sources, function(s) is.matrix(s$dat)))


  if(is.null(limitToCellLines)) {
    # if the list of cell lines was not explicitly provided, derive it based on the cell lines in the 
    # data matrices and the must.have.cellline flag
    
    must.have.cellline <- sapply(sources, function(s) {s$must.have.cellline} )
    # if there is at least one dataset which must have the cell lines, then take the intersection of all those
    # that have must.have.cellline == TRUE
    if(any(must.have.cellline)) {
      limitToCellLines <- foreach(source=sources[must.have.cellline], .combine=intersect) %do% {
        row.names(source$dat)
      }
    } else {
      # if there was no dataset which had must.have.cellline, then we take all cell lines
      limitToCellLines <- foreach(source=sources, .combine=union) %do% {
        row.names(source$dat)
      }
    }
  }

  # one way or another, this should be non-null by the time we reach this point
  stopifnot(!is.null(limitToCellLines))

  # generate the annotation table
  annotable <- foreach(source=sources, .combine=rbind) %do% {
    symbols <- colnames(source$dat)

    # SI feature names are not genes, so for this vartype, na out the gene symbols
    if(source$vartype == "SI") {
      symbols[] <- NA
    } else if(source$vartype == "Exp") {
      symbols <- sub(" \\(ENSG.*\\)$", "", symbols)
    }

    loci <- loci.by.symbol[symbols]

    stopifnot(source$vartype %in% c("CN","MUT","ACH","Exp","SI","SAN","GSE", "GS", "SS"))

    data <- data.frame(HUGOsymbol=symbols, VarType=source$vartype, ChrLoci=loci,row.names=paste(source$prefix, colnames(source$dat), sep="_"))
  }

  # generate the full matrix by normalizing the rows to be the cell lines chosen and concatenating the columns together
  cmat <- foreach (source=sources, .combine=cbind) %do% {
    sdat <- source$dat

    missing.celllines <- setdiff(limitToCellLines, row.names(col))
    if(length(missing.celllines)) {
      columns <- colnames(sdat)
      sdat <- rbind(sdat, matrix(NA, length(missing.celllines), length(columns), dimnames=list(missing.celllines, columns)))
    }

    data <- sdat[limitToCellLines,,drop=F]

    # generate full feature name: [feature type]_[feature name]
    colnames(data) <- paste(source$prefix, colnames(data), sep="_")

    data
  }

  list(matrix=cmat, annotable=annotable)
}

constructPredictionMatrix <- function(sources, loci.by.symbol, limitToCellLines=NULL) {
  ret <- constructMergedMatrix(sources, loci.by.symbol, limitToCellLines)
  
  ret$matrix <- removeNAcols(ret$matrix)
  sparse.columns <- columns.w.rare.values(ret$matrix)
  ret$matrix <- ret$matrix[,!sparse.columns,drop=F]
  
  ret
}

#' Construct a matrix and an annotation table for a set of matrices for use with ATLANTIS
#'
#' @param output.dir The directory to write the two matrices as Rdata files
#' @param pred.sources The definitions of the predictive features to use (See below)
#' @param target.source The definitions of the target matrix to use (See below)
#' @param loci.by.symbol A character vector where the elements are named by gene symbols
#'
#' @section Specifying data sources:
#' Both pred.sources and target.source are defined as a list of definitions where each definition has 
#' the following fields: dat, prefix, vartype, must.have.cellline
#' 
#' The 'dat' field is the numeric matrix of values.  'prefix' is a short code to prepend to each feature name from that matrix.  (For example "CN" for copy number to differentiate the feature from the gene expression for the same gene.)
#' 'vartype' must be one of 'MUT', 'Exp', 'CN' and controls how the values will be plotted.   Lastly, must.have.cellline controls whether NAs should be used to fill in features that are missing a 
#' cell line from this dataset, or whether the cell line should be dropped from the others if its missing.
#' 
#' Upon completion, returns a list with the following fields which point to the generated files: target.file, pred.file, anno.file
#' @export
writeAtlantisInputs <- function(output.dir,
                           pred.sources, 
                           target.sources,
                           loci.by.symbol=NULL,
                           analysisID="myAnalysis") {
    

  library(foreach)

  if(is.null(loci.by.symbol)) { 
    gene.names <- foreach(source=pred.sources, .combine=union) %do% colnames(source$dat)
    loci.by.symbol <- rep("NA", length(gene.names))
    names(loci.by.symbol) <- gene.names
  }

  pred <- constructPredictionMatrix(pred.sources, loci.by.symbol);
  target <- constructMergedMatrix(target.sources, loci.by.symbol);
  
  #Convert numeric vectors to factors for classification if all targets are 0 or 1.
  if (all(apply(target$matrix,2,function(x) { all(x %in% c(0,1,NA)) }))){
    cat("Target matrix is binary. \n")
    cat("Converting target matrix to factor data.frame for classification. \n")
    target$matrix <- data.frame(target$matrix)
    target$matrix <- data.frame(lapply(target$matrix,factor),row.names=rownames(target$matrix))
  }
  
  targetDims <- dim(target$matrix)
  predDims <- dim(pred$matrix)
  cat(paste("Target Matrix:", targetDims[2], "targets for", targetDims[1], "samples\n"))
  cat(paste("Prediction Matrix:", predDims[2], "features for", predDims[1], "samples\n"))
  
  # concatenate annotables
  anno <- pred$annotable
  anno <- rbind(anno, target$annotable[setdiff(row.names(target$annotable), row.names(anno)),])
  
  used.lines <- rownames(pred$matrix)
  
  shared.lines <- intersect(used.lines, rownames(target$matrix))

  target$matrix <- target$matrix[rownames(target$matrix) %in% shared.lines, ,drop=F] # limit to only lines used in predMat
  targetDims <- dim(target$matrix)
  cat(paste("After matching sample names: ", targetDims[2], "targets for ", targetDims[1], "samples\n"))
  stopifnot(is.matrix(target$matrix) || is.data.frame(target$matrix))
  
  predMat.file <- file.path(output.dir, paste(analysisID, '_PredMatrix.Rdata', sep=''))
  predMat <- pred$matrix
  stopifnot(all(dim(predMat) > 0))
  save(predMat, file=predMat.file)
  
  targetMat.file <- file.path(output.dir, paste(analysisID, '_targetMatrix.Rdata', sep=''))
  targetMat <- target$matrix 
  stopifnot(all(dim(targetMat) > 0))
  save(targetMat, file=targetMat.file)

  anno.file <- file.path(output.dir, paste(analysisID, '_anno.txt', sep=''))
  write.table(anno, file=anno.file, sep='\t', quote=FALSE)
  
  ret <- list(target.file=targetMat.file, pred.file=predMat.file, anno.file=anno.file)

  #if (generateTestMat) {
  #  testMat.file <- file.path(output.dir, paste(analysisID, '_testMatrix.txt', sep=''))
  #  write.table(testMat, file=testMat.file, sep='\t', quote=FALSE)
  #  ret <- c(ret, test.file=testMat.file)
  #}

  ret
}


#a <- matrix(c(1,2,3,4), 2, 2, dimnames=list(c("a","b"), c("c","d")))
#loci <- c("loci1", "loci2")
#names(loci) <- c("c","d")
#constructMatrix(list(list(dat=a, prefix="a", vartype="SI", must.have.cellline=F), list(dat=b, prefix="b", vartype="EXP", must.have.cellline=F)), loci)

# returns the input mat without all columns of  mat that are composed of solely NA/NaNs
removeNAcols <- function(mat) {
  colIDsToRemove <- which(simplify2array(apply(mat,2,function(x) all(is.na(x)))))
  
  cat(paste("removing", length(colIDsToRemove), "columns that are all NAs\n"))
  if (length(colIDsToRemove) > 0)
    mat[,-colIDsToRemove,drop=F]
  else
    mat
}


