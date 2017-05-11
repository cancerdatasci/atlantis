##############################################################################
#                                                                            #
#                           R-Code zur HiWi-Stelle                           #
#                                                                            #
#                                  Teil 2c:                                  #
#                     neue Funktionen zu Missing Values                      #
#                          FÜR REGRESSIONS-BÄUME!!!                          #
#                                                                            #
##############################################################################

# Funktion zum Erzeugen von fehlenden Werten per Logit-Modellierung *deleteLOG*:
# erzeugt in einer Liste von Datensätzen pro Zufallsvariable beliebig viele
# fehlende Werte

deleteLOG <- function(datsimul, mv = c(1, 0, 1, 1, 0),
    Coeff = list( c(-3, rep(1, times=5)), NULL, c(-3, rep(1, times=5)),
    c(-3, rep(1, times=5)) ), ...) {

    # Eingabe
    # datsimul: Liste von gleich großen Datensätzen (z.B. aus einer *dgp*-
    #   Funktion)
    # mv: Vektor als Indikator für fehlende Werte (missing values); 1 - feh-
    #   lende Werte einstreuen, 0 - Variable wird voll beobachtet (ohne MVs)
    # Coeff: Liste, die pro Eintrag einen Vektor zur zugehörigen Variable enthält,
    #   der die Berechnung der Logits steuert
    # Ausgabe
    # Liste von Datensätzen mit zufällig fehlenden Werten

    niter <- length(datsimul)          # Anzahl erzeugte Datensätze
    nvar <- length(mv)                 # Anzahl der Variablen x1 - xp
    n <- nrow(datsimul[[1]])           # Anzahl Beob eines erzeugten Datensatzes

    for (i in 1:niter) {
        x <- datsimul[[i]][, -1]       # speichert nur die Variablen ab
        logit.dep <- list()

        for (j in 1:nvar) {
            if(mv[j] > 0) {
                # Logit-Berechnung nur in den Variablen, die im Vektor mv angegeben sind

                coeff <- Coeff[[j]]

                ld <- coeff[1]         # Intercept
                for (k in 1:nvar) {
                    if (k != j) {
                        ld <- (coeff[k+1] * x[, k]) + ld  # [k+1] wegen Intercept
                        # ld ist Vektor
                        # z.B. P(X_1 missing) = logit( \sum_{k = 2}^n \beta_{1,k} * X_k) ~~ 0.2
                        # P(X_3 missing) = logit( \sum_{k = 1,2, 4,5} \beta_{3,k} * X_k) ~~ 0.1
                    }
                }
                logit.dep[[j]] <- ld   # logit.dep ist Liste von Vektoren
            }
        }

        for (j in 1:nvar) {
            if(mv[j] > 0) {
                # streichen nur in den Variablen, die im Vektor mv angegeben sind

                pNA <- as.vector( exp(logit.dep[[j]])/(1 + exp(logit.dep[[j]])) )
                # Wahrscheinlichkeit, in der j-ten Variable einen fehlenden Wert zu haben

                kriegtNA <- rbinom(n, 1, pNA)
                # n Bernoulli-Experimente, d.h. entweder 0 oder 1 als Ergebnis

                cat("Anteil NA in Variable", j, "Datensatz", i, ": ", sum(kriegtNA, na.rm=TRUE)/n, "\n")

                kriegtNA[kriegtNA == 1] <- NA
                kriegtNA <- 1 - kriegtNA
                # Damit bei * kriegtNA die ursprünglichen Werte erhalten bleiben,
                # wird der Vektor "umgedreht".
            }
        }
        datsimul[[i]][, -1] <- x
        # speichert Variablen mit fehlenden Werten zurück
    }
    datsimul
}

# Test:
#setwd("Z:/Eigene Dateien/HiWi/HiWi-Code")
#setwd("C:/HiWi-Code")
#source("funktionenZumDGP.R")
#setwd("C:/HiWi-Code/zusatzNA")
#source("funktionenZuMV2regres.R")
#d <- dgp2(niter=10, n=1000)
#dNA <- deleteLOG(d)


#****************************************************************************#

# Funktion zum Erzeugen von zufällig fehlenden Werten in Abhängigkeit vom
# Response *deleteDEPy*:
# erzeugt in einer Liste von Datensätzen pro Zufallsvariable _zufällig_ viele
# fehlende Werte in Abhängigkeit vom jeweiligen Response-Wert

deleteDEPy <- function(datsimul, mv = c(1, 0, 1, 1, 0), ...) {
    # Eingabe
    # datsimul: Liste von gleich großen Datensätzen (z.B. aus einer *dgp*-
    #   Funktion)
    # Funktioniert nur bei Klassifikation!
    # mv: Vektor als Indikator für fehlende Werte (missing values); 1 - feh-
    #   lende Werte einstreuen, 0 - Variable wird voll beobachtet (ohne MVs)
    # Ausgabe
    # Liste von Datensätzen mit zufällig fehlenden Werten nach MAR2

    niter <- length(datsimul)          # Anzahl erzeugte Datensätze
    nvar <- length(mv)                 # Anzahl der Variablen x1 - xp
    n <- nrow(datsimul[[1]])           # Anzahl Beob eines erzeugten Datensatzes
    mv <- mv * n

    for (i in 1:niter) {
        x <- datsimul[[i]][, -1]       # speichert nur die Variablen ab
        y <- datsimul[[i]][, 1]        # speichert nur den Response ab

        for (j in 1:nvar) {
            if(mv[j] > 0) {
                # streichen nur in den Variablen, die im Vektor mv angegeben sind

                p <- rep(0.1, times = n)
                # diese Zahl gibt Wahrscheinlichkeit für "1" an, also ist 0.1
                # die Wahrscheinlichkeit für kriegtNA == 1, also für NA, da es
                # ja später umgedreht wird.

                if (is.factor(y)) {
                    p[y == 2] <- 0.3
                    # diese Zahl gibt Wahrscheinlichkeit für "1" an, also ist
                    # 0.3 die Wahrscheinlichkeit für kriegtNA == 1, also für NA,
                    # da es ja später umgedreht wird.
                }

                if (!is.factor(y)) {
                    p[y < 13] <- 0.4
                    # diese Zahl gibt Wahrscheinlichkeit für "1" an, also ist
                    # 0.4 die Wahrscheinlichkeit für kriegtNA == 1, also für NA,
                    # da es ja später umgedreht wird.
                }

                kriegtNA <- rbinom(n, 1, p)
                # n Bernoulli-Experimente, d.h. entweder 0 oder 1 als Ergebnis

                cat("Anteil NA in Variable ", j, "Datensatz ", i, sum(kriegtNA)/n, "\n")

                kriegtNA[kriegtNA == 1] <- NA
                kriegtNA <- 1 - kriegtNA
                # Damit bei * kriegtNA die ursprünglichen Werte erhalten bleiben,
                # wird der Vektor "umgedreht".

                x[, j] <- kriegtNA * x[, j]
            }
        }
        datsimul[[i]][, -1] <- x
        # speichert Variablen mit fehlenden Werten zurück
    }
    datsimul
}

# Test:
#setwd("Z:/Eigene Dateien/HiWi/HiWi-Code")
#setwd("C:/HiWi-Code")
#source("funktionenZumDGP.R")
#setwd("C:/HiWi-Code/zusatzNA")
#source("funktionenZuMV2class.R")
#d <- dgp1(niter=10, n=100)
#dNA <- deleteDEPy(d)
#setwd("C:/HiWi-Code/zusatzNA")
#source("funktionenZuMV2regres.R")
#d <- dgp2(niter=10, n=100)
#dNA <- deleteDEPy(d)

##############################################################################