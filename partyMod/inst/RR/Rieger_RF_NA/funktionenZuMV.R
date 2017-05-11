##############################################################################
#                                                                            #
#                           R-Code zur HiWi-Stelle                           #
#                                                                            #
#                                  Teil 2a:                                  #
#                        Funktionen zu Missing Values                        #
#                                                                            #
##############################################################################

#****************************************************************************#
# folgende Funktion übernommen aus Svejdar (2007), allerdings Notation geän- #
# dert                                                                       #
#****************************************************************************#

# Funktion zum Erzeugen von komplett zufälligen Werten *deleteMCAR*:
# erzeugt in einer Liste von Datensätzen pro Zufallsvariable beliebig viele
# fehlende Werte

deleteMCAR <- function(datsimul, mv = c(0.2, 0, 0.1, 0.2, 0), ...) {
    # Eingabe
    # datsimul: Liste von gleich großen Datensätzen (z.B. aus einer *dgp*-
    #   Funktion)
    # mv: Vektor mit Anteilen fehlender Werte (missing values), die pro Zufalls-
    #   variable eingestreut werden sollen (pro Zufallsvariable eine Anzahl)
    # Ausgabe
    # Liste von Datensätzen mit komplett zufällig fehlenden Werten

    niter <- length(datsimul)          # Anzahl erzeugte Datensätze
    nvar <- length(mv)                 # Anzahl der Variablen x1 - xp
    n <- nrow(datsimul[[1]])           # Anzahl Beob eines erzeugten Datensatzes
    mv <- mv * n
    
    for (i in 1:niter) {
        x <- datsimul[[i]][, -1]       # speichert nur die Variablen ab
        for (j in 1:nvar) {
            if(mv[j] > 0) {
                # streichen nur in den Variablen, die im Vektor mv angegeben sind
                x[sample(n, mv[j]), j] <- NA
                # ersetzt in den Beobachtungen (n) die jeweils angegebene Anzahl
                # (mv[j]) an zufällig gezogenen Stellen durch NA
            }
        }
        datsimul[[i]][, -1] <- x
        # speichert Variablen mit fehlenden Werten zurück
    }
    datsimul
}

#****************************************************************************#
# alle folgenden Funktionen übernommen aus Svejdar (2007), allerdings Nota-  #
# tion geändert und zusätzliche Funktionen eingefügt, z. B. Abhängigkeit vom #
# Response                                                                   #
#****************************************************************************#

# erste Funktion zum Erzeugen von zufällig fehlenden Werten *deleteMAR1*:
# erzeugt in einer Liste von Datensätzen pro Zufallsvariable beliebig viele
# fehlende Werte durch Bildung von Rängen

deleteMAR1 <- function(datsimul, mv = matrix(data=c(0.2, 1, 2,
                                                    0.2, 4, 5,
                                                    0.1, 3, 4), byrow=TRUE, ncol=3),
    depY = FALSE, ...) {
    # Eingabe
    # datsimul: Liste von gleich großen Datensätzen (z.B. aus einer *dgp*-
    #   Funktion)
    # mv: Matrix (pro Zufallsvariable eine Zeile) mit Spalten 1-3
    #   * Anteil fehlender Werte (missing values), die in dieser Zufallsvariable
    #     eingestreut werden sollen
    #   * Variable, in der die Werte gestrichen werden sollen
    #   * Variable, die als Beurteilungskriterium dient
    #   Es muss eine Korrelation zwischen den beiden Variablen vorliegen.
    # Ausgabe
    # Liste von Datensätzen mit zufällig fehlenden Werten nach MAR1

    niter <- length(datsimul)          # Anzahl erzeugte Datensätze
    n <- nrow(datsimul[[1]])           # Anzahl Beob eines erzeugten Datensatzes
    mv[, 1] <- mv[, 1] * n
    
    if (depY == TRUE) {
        mv <- mv[, -3]
        # würde auch funktionieren, wenn mv nur zwei Spalten hätte
    }

    for (i in 1:niter) {
        x <- datsimul[[i]][, -1]       # speichert nur die Variablen ab
        y <- datsimul[[i]][, 1]        # speichert nur den Response ab
        for (j in 1:nrow(mv)) {        # mv ist jetzt eine Matrix
            if (depY == FALSE) {
                z <- rank(x[, mv[j, 3] ])
                # Hilfsvektor, vergibt Ränge nach der in der dritten Spalte von
                # mv angegebenen Variable
            }
            if (depY == TRUE) {
                z <- rank(y)
                # Hilfsvektor, vergibt Ränge nach dem Response
            }
            p <- z/sum(1:n)
            # Wahrscheinlichkeitsvektor; je kleiner der Rang, desto kleiner die
            # Wahrscheinlichkeit, dass der Wert in der anderen Variable (mv,
            # Spalte 2) gestrichen zu werden
            x[sample(n, mv[j, 1], prob = p), mv[j, 2] ] <- NA
            # ersetzt in den Beobachtungen (n) der Variable, in der die Werte
            # fehlen sollen (mv[j], Spalte 2), die jeweils angegebene Anzahl
            # (mv[j], Spalte 1) an zufällig gezogenen Stellen mit Wahrschein-
            # lichkeit p durch NA
        }
        datsimul[[i]][, -1] <- x
        # speichert Variablen mit fehlenden Werten zurück
    }
    datsimul
}

# Test:
#setwd("C:/HiWi-Code")
#source("funktionenZumDGP.R")
#source("funktionenZuMV.R")
#d <- dgp1(niter=10, n=50)
#dNA <- deleteMAR1(d, depY = TRUE)      # funktioniert!
#d <- dgp2(niter=10, n=50)
#dNA <- deleteMAR1(d, depY = TRUE)

#****************************************************************************#

# zweite Funktion zum Erzeugen von zufällig fehlenden Werten *deleteMAR2*:
# erzeugt in einer Liste von Datensätzen pro Zufallsvariable beliebig viele
# fehlende Werte durch Bildung von zwei Risikogruppen

deleteMAR2 <- function(datsimul, mv = matrix(data=c(0.2, 1, 2,
                                                    0.2, 4, 5,
                                                    0.1, 3, 4), byrow=TRUE, ncol=3),
    depY = FALSE, ...) {
    # Eingabe
    # datsimul: Liste von gleich großen Datensätzen (z.B. aus einer *dgp*-
    #   Funktion)
    # mv: Matrix (pro Zufallsvariable eine Zeile) mit Spalten 1-3
    #   * Anteil fehlender Werte (missing values), die in dieser Zufallsvariable
    #     eingestreut werden sollen
    #   * Variable, in der die Werte gestrichen werden sollen
    #   * Variable, die als Beurteilungskriterium dient
    #   Es muss eine Korrelation zwischen den beiden Variablen vorliegen.
    # Ausgabe
    # Liste von Datensätzen mit zufällig fehlenden Werten nach MAR2

    niter <- length(datsimul)          # Anzahl erzeugte Datensätze
    n <- nrow(datsimul[[1]])           # Anzahl Beob eines erzeugten Datensatzes
    mv[, 1] <- mv[, 1] * n

    if (depY == TRUE) {
        mv <- mv[, -3]
        # würde auch funktionieren, wenn mv nur zwei Spalten hätte
    }

    for (i in 1:niter) {
        x <- datsimul[[i]][, -1]       # speichert nur die Variablen ab
        y <- datsimul[[i]][, 1]        # speichert nur den Response ab
        for (j in 1:nrow(mv)) {        # mv ist jetzt eine Matrix
            z <- rep(0, n)
            if (depY == FALSE) {
                z[ x[, mv[j, 3]] >= median(x[, mv[j, 3]]) ] <- 1
                # Hilfsvektor mit 0-1-Kodierung: falls der jeweilige Eintrag in
                # der Beurteilungsspalte von x (mv, 3. Spalte) größer/gleich als
                # der Median der Beurteilungsspalte ist, wird der Eintrag in z
                # "1", sonst bleibt er Null
            }
            if (depY == TRUE) {
                z[y >= median(y)] <- 1
                # Hilfsvektor mit 0-1-Kodierung: falls der jeweilige Eintrag im
                # Response größer/gleich dem Median des Response ist, wird der
                # Eintrag in z "1", sonst bleibt er Null
            }
            S <- sum(z)                # zählt die Einträge in z, die "1" sind
            p <- rep(0.1/(n - S), n)
            p[z == 1] <- 0.9/(n - S)
            # Wahrscheinlichkeitsvektor; falls der Eintrag in z "1" ist, wird
            # die Wahrscheinlichkeit, dass der Wert in der anderen Variable (mv,
            # Spalte 2) gestrichen wird, um einen festen Faktor (hier: 9)
            # größer. Sonst bleibt er klein. Damit haben die hohen Werte der Be-
            # urteilungsvariable eine größere Wahrscheinlichkeit auf fehlende
            # Werte in der Streich-Variable.
            # Wahrscheinlichkeit wird berechnet durch 0.1 bzw. 0.9 durch Anzahl
            # der zu streichenden Werte in dieser Gruppe, so dass die Gesamt-
            # Wahrscheinlichkeit in der Gruppe wieder 0.1 bzw. 0.9 ist.
            x[sample(n, mv[j, 1], prob=p), mv[j, 2] ] <- NA
            # ersetzt in den Beobachtungen (n) die jeweils angegebene Anzahl
            # (mv[j], Spalte 1) an zufällig gezogenen Stellen mit Wahrschein-
            # lichkeit p durch NA
        }
        datsimul[[i]][, -1] <- x
        # speichert Variablen mit fehlenden Werten zurück
    }
    datsimul
}

# Test:
#setwd("C:/HiWi-Code")
#source("funktionenZumDGP.R")
#source("funktionenZuMV.R")
#d <- dgp1(niter=10, n=50)
#dNA <- deleteMAR2(d, depY = TRUE)      # funktioniert nicht, braucht numerische Daten
#d <- dgp2(niter=10, n=50)
#dNA <- deleteMAR2(d, depY = TRUE)

#****************************************************************************#

# dritte Funktion zum Erzeugen von zufällig fehlenden Werten *deleteMAR3*:
# erzeugt in einer Liste von Datensätzen pro Zufallsvariable beliebig viele
# fehlende Werte durch rechtsseitige Trunkierung

deleteMAR3 <- function(datsimul, mv = matrix(data=c(0.2, 1, 2,
                                                    0.2, 4, 5,
                                                    0.1, 3, 5), byrow=TRUE, ncol=3),
    depY = FALSE, ...) {
    # Eingabe
    # datsimul: Liste von gleich großen Datensätzen (z.B. aus einer *dgp*-
    #   Funktion)
    # mv: Matrix (pro Zufallsvariable eine Zeile) mit Spalten 1-3
    #   * Anteil fehlender Werte (missing values), die in dieser Zufallsvariable
    #     eingestreut werden sollen
    #   * Variable, in der die Werte gestrichen werden sollen
    #   * Variable, die als Beurteilungskriterium dient
    #   Es muss eine Korrelation zwischen den beiden Variablen vorliegen.
    # Ausgabe
    # Liste von Datensätzen mit zufällig fehlenden Werten nach MAR3

    niter <- length(datsimul)          # Anzahl erzeugte Datensätze
    n <- nrow(datsimul[[1]])           # Anzahl Beob eines erzeugten Datensatzes
    mv[, 1] <- mv[, 1] * n

    if (depY == TRUE) {
        mv <- mv[, -3]
        # würde auch funktionieren, wenn mv nur zwei Spalten hätte
    }

    for (i in 1:niter) {
        x <- datsimul[[i]][, -1]       # speichert nur die Variablen ab
        y <- datsimul[[i]][, 1]        # speichert nur den Response ab
        for (j in 1:nrow(mv)) {        # mv ist jetzt eine Matrix
            if (depY == FALSE) {
                a <- quantile(x[, mv[j, 3] ], probs = (1 - (mv[j, 1]/n)))
                # berechnet das Quantil in der Beurteilungsvariable (mv, Spalte
                # 3), sodass die Anzahl an fehlenden Werten (mv, Spalte 1; nach
                # Berechnung innerhalb dieser Funktion) durch die Gesamt-Beo-
                # bachtungsanzahl eine Prozentzahl ergibt. Diese wird von 100%
                # = 1 abgezogen.
                # Achtung: Fehlermeldung, wenn in der Beurteilungsvariable MV's
                # sind und na.rm=FALSE (default)!
                # Es werden dann diejenigen Werte in der Streichvariable (mv,
                # Spalte 2) gestrichen, die in der Beurteilungsvariable über
                # diesem Quantil liegen:
                x[, mv[j, 2]][x[, mv[j, 3]] >= a] <- NA
            }
            if (depY == TRUE) {
                a <- quantile(y, probs = (1 - (mv[j, 1]/n)))
                # berechnet das Quantil im Response, sodass die Anzahl an feh-
                # lenden Werten (mv, Spalte 1; nach Berechnung innerhalb dieser
                # Funktion) durch die Gesamt-Beobachtungsanzahl eine Prozent-
                # zahl ergibt. Diese wird von 100% = 1 abgezogen.
                # Es werden dann diejenigen Werte in der Streichvariable (mv,
                # Spalte 2) gestrichen, die im Response über diesem Quantil
                # liegen:
                x[, mv[j, 2]][y >= a] <- NA
            }
        }
        datsimul[[i]][, -1] <- x
        # speichert Variablen mit fehlenden Werten zurück
    }
    datsimul
}

# Test:
#setwd("C:/HiWi-Code")
#source("funktionenZumDGP.R")
#source("funktionenZuMV.R")
#d <- dgp1(niter=10, n=50)
#dNA <- deleteMAR3(d, depY = TRUE)      # funktioniert nicht, gibt warnings aus
#d <- dgp2(niter=10, n=50)
#dNA <- deleteMAR3(d, depY = TRUE)

#****************************************************************************#

# vierte Funktion zum Erzeugen von zufällig fehlenden Werten *deleteMAR4*:
# (im Original (Svedjar, 2007): *deleteMAR5*)
# erzeugt in einer Liste von Datensätzen pro Zufallsvariable beliebig viele
# fehlende Werte durch symmetrische Turnkierung

deleteMAR4 <- function(datsimul, mv = matrix(data=c(0.2, 1, 2,
                                                    0.2, 4, 5,
                                                    0.1, 3, 5), byrow=TRUE, ncol=3),
    depY = FALSE, ...) {
    # Eingabe
    # datsimul: Liste von gleich großen Datensätzen (z.B. aus einer *dgp*-
    #   Funktion)
    # mv: Matrix (pro Zufallsvariable eine Zeile) mit Spalten 1-3
    #   * Anteil fehlender Werte (missing values), die in dieser Zufallsvariable
    #     eingestreut werden sollen
    #   * Variable, in der die Werte gestrichen werden sollen
    #   * Variable, die als Beurteilungskriterium dient
    #   Es muss eine Korrelation zwischen den beiden Variablen vorliegen.
    # Ausgabe
    # Liste von Datensätzen mit zufällig fehlenden Werten nach MAR4

    niter <- length(datsimul)          # Anzahl erzeugte Datensätze
    n <- nrow(datsimul[[1]])           # Anzahl Beob eines erzeugten Datensatzes
    mv[, 1] <- mv[, 1] * n

    if (depY == TRUE) {
        mv <- mv[, -3]
        # würde auch funktionieren, wenn mv nur zwei Spalten hätte
    }

    for (i in 1:niter) {
        x <- datsimul[[i]][, -1]       # speichert nur die Variablen ab
        y <- datsimul[[i]][, 1]        # speichert nur den Response ab
        for (j in 1:nrow(mv)) {        # mv ist jetzt eine Matrix
            if (depY == FALSE) {
                a <- quantile(x[, mv[j, 3]], prob = (1 - (0.5 * mv[j, 1]/n)))
                # berechnet das Quantil in der Beurteilungsvariable (mv,
                # Spalte 3), sodass die halbe Anzahl an fehlenden Werten (mv,
                # Spalte 1) durch die Gesamt-Beobachtungsanzahl eine Prozent-
                # zahl ergibt. Diese wird von 100% = 1 abgezogen.
                # Es werden dann diejenigen Werte in der Streichvariable (mv,
                # Spalte 2) gestrichen, die in der Beurteilungsvariable über
                # diesem Quantil liegen:
                x[, mv[j, 2]][x[, mv[j, 3]] >= a] <- NA

                b <- quantile(x[, mv[j, 3]], prob = (0.5 * mv[j, 1]/n))
                # berechnet das Quantil in der Beurteilungsvariable (mv,
                # Spalte 3), sodass die halbe Anzahl an fehlenden Werten (mv,
                # Spalte 1) durch die Gesamt-Beobachtungsanzahl eine Prozent-
                # zahl ergibt.
                # Es werden dann diejenigen Werte in der Streichvariable (mv,
                # Spalte 2) gestrichen, die in der Beurteilungsvariable unter
                # diesem Quantil liegen:
                x[, mv[j, 2]][x[, mv[j, 3]] <= b] <- NA
            }
            if (depY == TRUE) {
                a <- quantile(y, prob = (1 - (0.5 * mv[j, 1]/n)))
                # berechnet das Quantil im Response, sodass die halbe Anzahl
                # an fehlenden Werten (mv, Spalte 1) durch die Gesamt-Beobach-
                # tungsanzahl eine Prozentzahl ergibt. Diese wird von 100% = 1
                # abgezogen.
                # Es werden dann diejenigen Werte in der Streichvariable (mv,
                # Spalte 2) gestrichen, die im Response über diesem Quantil
                # liegen:
                x[, mv[j, 2]][y >= a] <- NA

                b <- quantile(y, prob = (0.5 * mv[j, 1]/n))
                # berechnet das Quantil im Response, sodass die halbe Anzahl
                # an fehlenden Werten (mv, Spalte 1) durch die Gesamt-Beobach-
                # tungsanzahl eine Prozentzahl ergibt.
                # Es werden dann diejenigen Werte in der Streichvariable (mv,
                # Spalte 2) gestrichen, die im Response unter diesem Quantil
                # liegen:
                x[, mv[j, 2]][y <= b] <- NA
            }

        }
        datsimul[[i]][, -1] <- x
        # speichert Variablen mit fehlenden Werten zurück
    }
    datsimul
}

# Test:
#setwd("C:/HiWi-Code")
#source("funktionenZumDGP.R")
#source("funktionenZuMV.R")
#d <- dgp1(niter=10, n=50)
#dNA <- deleteMAR4(d, depY = TRUE)      # funktioniert nicht, gibt warnings aus
#d <- dgp2(niter=10, n=50)
#dNA <- deleteMAR4(d, depY = TRUE)

##############################################################################
# Literatur:                                                                 #
#                                                                            #
# * Viola Svejdar: "Variablenselektion in Klassifikationsbäumen unter spezi- #
#   eller Berücksichtigung von fehlenden Werten", Diplomarbeit, Ludwig-Maxi- #
#   milians-Universität, München, 2007                                       #
##############################################################################