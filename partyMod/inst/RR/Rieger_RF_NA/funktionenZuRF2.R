##############################################################################
#                                                                            #
#                           R-Code zur HiWi-Stelle                           #
#                                                                            #
#                                  Teil 3:                                   #
#                     Funktionen zu Random Forests                           #
#                                                                            #
##############################################################################

mylog <- function(x) {
    x[x < 1e-5] <- 1e-5                # nötig wegen log(0) = -Inf
    log(x)
}

#****************************************************************************#

# Funktion zum Berechnen der Binomial-Loglikelihood
# für der Daten generierenden Prozess 1 *RF1*:
# erzeugt - basierend auf einer Liste von gleich großen Datensätzen - einen Ran-
# dom Forest, sagt die Werte für einen Testdatensatz vorher und berechnet die
# mittlere Binomial-Loglikelihood

RF1 <- function(mvfun = NULL, test, nIter = 500,
    lernMV = FALSE, testMV = FALSE, imp = FALSE, na.omit = FALSE, ...) {
    # Eingabe
    # mvfun: eine der *delete*-Funktionen; nur nötig, falls auch fehlende Werte
    #   eingestreut werden
    # test: Test-Datensatz ohne fehlende Werte als Liste (mit einem Element)
    # nIter: Anzahl der zu erzeugenden Lerndatensätze
    # lernMV, testMV: logischer Wert: TRUE, falls im Lern- bzw. Testdatensatz
    #   fehlende Werte eingestreut werden sollen
    # mvlern, mvtest: Vektor mit den Anzahlen an fehlenden Werten (für MCAR)
    #   bzw. Matrix mit Anzahl fehlender Werte, Streichvariable, Beurteilungs-
    #   variable (pro Variable mit fehlenden Werten eine Zeile); bei fehlender
    #   Angabe wird der Default aus den *delete*-Funktionen verwendet (d.h.
    #   keine fehlenden Werte)
    # mögliche weitere Eingaben für dgpfun:
    #   n, sigma, evtl. coef (für *dgp1*)
    # Ausgabe
    # loglik: Vektor mit der Binomial-Log-Likelihood pro erzeugtem Lerndatensatz

    dat <- dgp1(niter = nIter, n = 200, ...)

    # fehlende Werte in den Lern-Datensatz einstreuen:
    if (lernMV) {
        dat <- mvfun(dat, ...)
        cat("* NA in die Lern-Datensaetze: erledigt \n")
    }

    # fehlende Werte in den Test-Datensatz einstreuen:
    if (testMV) {
        test <- mvfun(test, ...)
        cat("* NA in den Test-Datensatz: erledigt \n")
    }
    test <- test[[1]]

    # Imputation des Test-Datensatzes:
    if (imp & testMV) {
        stopifnot(require("impute"))
            # aus: "Missing value estimation methods for DNA microarrays" von
            # Troyanskaya et al., 2001
        test1 <- test[, -1]
        ###save(test1, file = "test1.Rda")
        test2  <- try(as.data.frame(impute.knn(as.matrix(test1))))
        if (inherits(test2, "try-error")) {
            save(test2, file = "errortest.Rda")
            stop("Impute.knn hat Fehler geschmissen.")
        }
        test[, -1] <- test2
    }

    stopifnot(require("party"))

    loglik <- vector("numeric")
    niter <- length(dat)

    for (i in 1:niter) {
        loglik[i] <- NA

        # Imputation der Lern-Datensätze:
        if (imp & lernMV) {
            dat[[i]][-1]  <- as.data.frame(impute.knn(as.matrix(dat[[i]][-1])))

#            dat[[i]][-1]  <- try( as.data.frame(impute.knn(as.matrix(dat[[i]][-1]))) )
#            if (inherits(dat[[i]][-1], "try-error")) next();
        ######
        #print(summary(dat[[i]]))
        #dati <- dat[[i]]
        #save(dati, file = paste("dat", i, ".Rda", sep=""))
        ######
        }

        if (na.omit) dat[[i]] <- dat[[i]][complete.cases(dat[[i]]), ]

        # Berechnung des Random Forests:
        rs <- .Random.seed
        rf <- try(cforest(y ~ ., data = dat[[i]],
            control = cforest_control(maxsurrogate=3, ntree=50, minsplit=30)))
        if (inherits(rf, "try-error")) {
            #save(rf, rs, file = "error.Rda")
            #save(test, file = "errortest.Rda")
            #save(dat, file = "errordat.Rda")
            print("Random Forest hat Fehler geschmissen.")
            loglik[i] <- NA ### einfach NA als loglik einfuegen
            next()          ### naechstes i
        }

        print(paste("RF", i))

        # Vorhersage für den Test-Datensatz:
        p <- treeresponse(rf, newdata = test)
        P <- matrix(unlist(p), byrow = TRUE, ncol = 2)
        p <- P[, 2]
        yvec <- as.numeric(test$y) - 1
        loglik[i] <- mean(yvec * mylog(p) + (1 - yvec) * mylog(1 - p))
    }
    return(loglik)
}

# Test:
#test <- dgp1(niter=1, n=5000)
#dat <- RF1(deleteMCAR, test, dgpfun.niter=10, lernMV=TRUE, testMV=TRUE,
#           mvlern=c(0.2*200, 0.2*200, 0.1*200, 0, 0),
#           mvtest=c(0.2*5000, 0.2*5000, 0.1*5000, 0, 0), imp=TRUE)
#boxplot(dat, ylim=c(-15000, 0))
#rm(test, dat)

#****************************************************************************#

# Funktion zum Berechnen des Mean Squared Errors
# für der Daten generierenden Prozess 2 *RF2*:
# erzeugt - basierend auf einer Liste von gleich großen Datensätzen - einen Ran-
# dom Forest, sagt die Werte für einen Testdatensatz vorher und berechnet den
# Mean Squared Error

RF2 <- function(mvfun = NULL, test, nIter = 500,
    lernMV = FALSE, testMV = FALSE, imp = FALSE, na.omit = FALSE, ...) {
    # Eingabe
    # mvfun: eine der *delete*-Funktionen; nur nötig, falls auch fehlende Werte
    #   eingestreut werden
    # test: Test-Datensatz ohne fehlende Werte als Liste (mit einem Element)
    # dgpfun.niter: Anzahl der zu erzeugenden Lerndatensätze
    # lernMV, testMV: logischer Wert: TRUE, falls im Lern- bzw. Testdatensatz
    #   fehlende WErte eingestreut werden sollen
    # mvlern, mvtest: Vektor mit den Anzahlen an fehlenden Werten (für MCAR)
    #   bzw. Matrix mit Anzahl fehlender Werte, Streichvariable, Beurteilungs-
    #   variable (pro Variable mit fehlenden Werten eine Zeile); bei fehlender
    #   Angabe wird der Default aus den *delete*-Funktionen verwendet (d.h.
    #   keine fehlenden Werte)
    # mögliche weitere Eingaben für dgpfun:
    #   n, sigma
    # Ausgabe
    # MSE: Vektor mit dem Mean Squared Error pro erzeugtem Lerndatensatz

    dat <- dgp2(niter = nIter, n = 200, ...)

    # fehlende Werte in die Lern-Datensaetze einstreuen:
    if (lernMV) {
        dat <- mvfun(dat, ...)
        cat("* NA in die Lern-Datensaetze: erledigt \n")
    }

    # fehlende Werte in den Test-Datensatz einstreuen:
    if (testMV) {
        test <- mvfun(test, ...)
        cat("* NA in den Test-Datensatz: erledigt \n")
    }
    test <- test[[1]]

    # Imputation des Test-Datensatzes:
    if (imp & testMV) {
        stopifnot(require("impute"))
            # aus: "Missing value estimation methods for DNA microarrays" von
            # Troyanskaya et al., 2001
        test[, -1]  <- as.data.frame(impute.knn(as.matrix(test[, -1])))
    }

    stopifnot(require("party"))

    MSE <- vector("numeric")
    niter <- length(dat)

    for (i in 1:nIter) {
        MSE[i] <- NA

        # Imputation der Lern-Datensätze:
        if (imp & lernMV) {
            dat[[i]][-1]  <- as.data.frame(impute.knn(as.matrix(dat[[i]][-1])))
        }

        if (na.omit) dat[[i]] <- dat[[i]][complete.cases(dat[[i]]), ]

        # Berechnung des Random Forests:
        rs <- .Random.seed
        rf <- try(cforest(y ~ ., data = dat[[i]],
            control = cforest_control(maxsurrogate=3, ntree=50, minsplit=30)))
        if (inherits(rf, "try-error")) {
            #save(rf, rs, file = "error.Rda")
            #save(test, file = "errortest.Rda")
            #save(dat, file = "errordat.Rda")
            print("Random Forest hat Fehler geschmissen.")
            MSE[i] <- NA ### einfach NA als loglik einfuegen
            next()          ### naechstes i

        }

        print(paste("RF", i))

        # Vorhersage für den Test-Datensatz:
        f <- predict(rf, newdata=test)
        MSE[i] <- mean((test$y - f)^2)
    }
    return(MSE)
}

# Test:
#test <- dgp2(niter=1, n=5000)
#dat <- RF2(deleteMCAR, test, dgpfun.niter=10, lernMV=TRUE, testMV=TRUE,
#           mvlern=c(0.2*200, 0.2*200, 0.1*200, 0, 0),
#           mvtest=c(0.2*5000, 0.2*5000, 0.1*5000, 0, 0), imp=TRUE)
#boxplot(dat, ylim=c(0, 11))
#rm(test, dat)

##############################################################################