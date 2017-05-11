
# library("mboost")
if (FALSE) {
cvrisk2 <- function(object, folds = cv(object@weights), 
    mincriterion = seq(from = 0.6, to = 0.95, by = 0.05),
    papply = if (require("multicore")) mclapply else lapply,
    ...){

    weights <- object@weights

    if (any(weights == 0))
        warning("zero weights")
    if (is.null(folds)) {
        folds <- rmultinom(25, length(weights), weights/sum(weights))
    } else {
        stopifnot(is.matrix(folds) && nrow(folds) == length(weights))
    }
    fitfct <- object@update
    oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))

    dummyfct <- function(weights) {
        mod <- fitfct(weights = weights)
        sapply(mincriterion, function(m) {
            p <- predict(mod, mincriterion = m)
            if (is.factor(p)) 
                err <- p != mod@responses@variables[[1]]
            if (is.numeric(p))
                err <- (p - mod@responses@variables[[1]])^2
            mean(err[weights == 0])
        })
    }

    OOBweights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
    OOBweights[folds > 0] <- 0
    oobrisk <- papply(1:ncol(folds),
        function(i) dummyfct(weights = folds[, i]))

    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk / colSums(OOBweights)
    colnames(oobrisk) <- mincriterion
    rownames(oobrisk) <- 1:nrow(oobrisk)
    oobrisk
}

set.seed(290875)
     
### regression
airq <- subset(airquality, !is.na(Ozone))
airct <- ctree(Ozone ~ ., data = airq, 
               controls = ctree_control(maxsurrogate = 3, mincriterion = 0.5))
### bootstrap
cvm <- cvrisk2(airct, mincriterion = 50:95/100)
boxplot(cvm)

### 10-fold CV
cvm <- cvrisk2(airct, mincriterion = 50:95/100, cv(airct@weights, type = "kfold"))
boxplot(cvm)

}

