##############################################################################
#                                                                            #
#                           R-Code zur HiWi-Stelle                           #
#                                                                            #
#                                  Teil 1:                                   #
#                Funktionen zum datengenerierenden Prozess                   #
#                                                                            #
##############################################################################

#library(mvtnorm, lib.loc = "~/lib/oldmvtnorm/")

# Funktion zum datengenerierenden Prozess 1 *dgp1*:
# erzeugt niter Datensätze mit je n Beobachtungen aus 5 multivariat normalver-
# teilten Zufallsvariablen mit Erwartungswert Null, beliebiger Varianz
# ("sigma") und dem binären Response y (ein Faktor). Der Einfluss der einzel-
# nen Kovariablen kann mittels "coef" beliebig verändert werden. Zusätzlich
# werden w multivariat normalverteile Zufallsvariablen (ebenfalls je n Beo-
# bachtungen) gezogen, welche keinen Einfluss auf den Response haben.

dgp1 <- function(niter = 5, n = 100, sigma = diag(rep(1, 5)), coef = 1:5, w = 0) {
    #************************************************************************#
    # aufbauend auf der Funktion "create" aus Svejdar (2007)                 #
    #************************************************************************#
    #
    # Eingabe
    # niter: Anzahl der Datensätze
    # n: Anzahl der Beobachtungen je Datensatz
    # sigma: Kovarianzmatrix, default: Einheitsmatrix
    # coef: Einfluss der einzelnen Kovariablen x_i
    # w: Anzahl der Noise-Parameter
    # Ausgabe
    # datsimul: Liste der Datensätze

    stopifnot(require("mvtnorm"))
    # lädt Paket "mvtnorm", falls noch nicht geschehen
    datsimul <- list()

    for (i in 1:niter) {
        # wahre Parameter:
        X <- rmvnorm(n, mean = rep(0, 5), sigma = sigma)
        # der Erwartungswert ist Null
        
        # Response:
        pi <- as.vector( exp(X %*% coef)/(1 + exp(X %*% coef)) )
        # Wahrscheinlichkeit, in Klasse 2 zu fallen; nach Formel (4.1) aus
        # Svejdar (2007)
        y <- rbinom(n, 1, pi) + 1
        # +1, da nicht 0-1-kodiert, sondern 1-2
            
        # mit Noise-Parameter:
        if (w > 0) {
            S <- matrix(data = 0, ncol = max(w, 15), nrow = max(w, 15))
            S[1:5, 1:5] <- matrix(data = c(  1, 0.9, 0.9, 0.9, 0.9,
                                           0.9,   1, 0.9, 0.9, 0.9,
                                           0.9, 0.9,   1, 0.9, 0.9,
                                           0.9, 0.9, 0.9,   1, 0.9,
                                           0.9, 0.9, 0.9, 0.9,   1),
                byrow = TRUE, ncol = 5)
            S[6:10, 6:10] <- matrix(data = c(  1, 0.9, 0.9,   0,   0,
                                             0.9,   1, 0.9,   0,   0,
                                             0.9, 0.9,   1,   0,   0,
                                               0,   0,   0,   1, 0.9,
                                               0,   0,   0, 0.9,   1),
                byrow = TRUE, ncol = 5)
            S[11:15, 11:15] <- matrix(data = c(  1, 0.1, 0.1, 0.1, 0.1,
                                               0.1,   1, 0.1, 0.1, 0.1,
                                               0.1, 0.1,   1, 0.1, 0.1,
                                               0.1, 0.1, 0.1,   1, 0.1,
                                               0.1, 0.1, 0.1, 0.1,   1),
                byrow = TRUE, ncol = 5)
            S <- S[1:w, 1:w]
            W <- rmvnorm(n, mean = rep(0, w), sigma = S)

            # Datensatz:
            dat <- as.data.frame(cbind(y, X, W))
            dat$y <- as.factor(y)
            names(dat) <- c("y", paste("x", 1:5, sep = ""),
                paste("w", 1:w, sep = ""))
        }

        # ohne Noise-Parameter:
        if (w == 0) {
            # Datensatz:
            dat <- as.data.frame(cbind(y, X))
            dat$y <- as.factor(y)
            names(dat) <- c("y", paste("x", 1:5, sep = ""))
        }
        if (w < 0) {   stop("Falsche Anzahl an Noise-Parametern w: ", w)   }

        datsimul[[i]] <- dat
    }
    return(datsimul)
}

# Test:
#d <- dgp1(n=30)
#class(d[[1]]$y)
#rm(d)

#****************************************************************************#

# Funktion zum datengenerierenden Prozess 2 *dgp2*:
# erzeugt niter Datensätze mit je n Beobachtungen aus 5 multivariat gleichver-
# teilten Zufallsvariablen auf [0, 1], beliebiger Varianz ("sigma") und dem
# Response y, der nach dem Regressionsproblem "Friedman 1" aus Friedman (1991)
# berechnet wird. Zusätzlich werden u multivariat gleichverteile Zufallsvariab-
# len (ebenfalls je n Beobachtungen) gezogen, welche keinen Einfluss auf den
# Response haben.
# Für Details siehe:
# R> library(mlbench)
# R> ?mlbench.friedman1

dgp2 <- function(niter = 5, n = 100, sigma = diag(rep(1, 5)), u = 5) {
    # Eingabe
    # niter: Anzahl der Datensätze
    # n: Anzahl der Beobachtungen je Datensatz
    # sigma: Kovarianzmatrix, default: Einheitsmatrix
    # u: Anzahl der Noise-Parameter
    # Ausgabe
    # datsimul: die niter Datensätze

    stopifnot(require("mvtnorm"))
    # lädt Paket "mvtnorm", falls noch nicht geschehen
    datsimul <- list()
    if (u < 0) {   stop("Falsche Anzahl an Noise-Parametern u: ", u)   }

    for (i in 1:niter) {
        # Die gezogenen Zufallsvariablen sind gleichverteilt:
        # Vektor p speichert die p-Werte des KS-Tests. Falls ein Element von p
        # kleiner als 0.05 ist, wird die Nullhypothese U[,j] ~ U([0, 1]) abge-
        # lehnt. Damit mindestens einmal Daten gezogen werden, wird p auf 0.01
        # gesetzt.
        
        # wahre Parameter:
        p <- rep(0.01, 5)
        while(all(p <= 0.05)) {
            # hier Überprüfung, ob alle U[,j]'s gleichverteilt sind; falls ja,
            # wird die while-Schleife verlassen
            
            X <- rmvnorm(n, mean = rep(0, 5), sigma = sigma)
                # der Erwartungswert ist Null
            U <- apply(X, 2, pnorm)
                # 2 für spaltenweise Anwendung von apply
            for (j in 1:5) {
                p[j] <- ks.test(U[, j], y = punif, exact = FALSE)$p.value
                    # bei Bindungen ist exact = FALSE nötig
            }
        }
        # Response:
        y <- 10*sin(pi*U[, 1]*U[, 2]) + 20*(U[, 3] - 0.5)^2 + 10*U[, 4] + 5*U[, 5]
            # Berechnung nach Friedman1 (s. o.)
        # Noise-Parameter:
        p <- rep(0.01, 5)
        while(all(p <= 0.05)) {
            # hier Überprüfung, ob alle V[,j]'s gleichverteilt sind; falls ja,
            # wird die while-Schleife verlassen
            
            W <- rmvnorm( n, mean = rep(0, u), sigma = diag(rep(1, u)) )
            V <- apply(W, 2, pnorm)
                # 2 für spaltenweise Anwendung von apply
            for (j in 1:u) {
                p[j] <- ks.test(V[,j], y=punif)$p.value
            }
        }

        dat <- as.data.frame(cbind(y, U, V))
        names(dat) <- c("y", paste("x", 1:5, sep = ""),
            paste("u", 1:u, sep = ""))
        datsimul[[i]] <- dat
    }
    return(datsimul)
}

# Test:
#dgp2(n = 30, sigma = sigma[[2]], u = 25)

##############################################################################
# Literatur:                                                                 #
#                                                                            #
# * Viola Svejdar: "Variablenselektion in Klassifikationsbäumen unter spezi- #
#   eller Berücksichtigung von fehlenden Werten", Diplomarbeit, Ludwig-Maxi- #
#   milians-Universität, München, 2007                                       #
# * Friedman, Jerome H.: "Multivariate adaptive regression splines", 1991    #
#   The Annals of Statistics 19 (1), pages 1-67.                             #
##############################################################################