rm(list = ls())
R.Version()$version.string

library("impute", lib.loc = "~/lib/oldR")
packageDescription("impute")$Version
# altes impute lassen, da der eine Bug nicht gefixt ist.
library("mvtnorm")
packageDescription("mvtnorm")$Version
library("party")
packageDescription("party")$Version

##############################################################################
#                                                                            #
#                           R-Code zur HiWi-Stelle                           #
#                                                                            #
#                                  Teil 9c:                                  #
#                Berechnungen der zusätzlichen NA-Mechanismen                #
#                       *** auf dem biostat-Server ***                       #
#                                                                            #
##############################################################################

# Klassifikation:
#----------------

source("funktionenZuMV2class.R")
source("Simulationen2-alt-neu.R")
source("funktionenZuMV2class.R")

packageDescription("mvtnorm")$Version
packageDescription("impute")$Version
packageDescription("party")$Version

#NA2res1 <- foo(RF = RF1, dgp = dgp1, ntest = 5000, dgpfun.niter = 5)
NA2res1 <- foo(RF = RF1, dgp = dgp1, ntest = 5000, dgpfun.niter = 100)
# ntest im Test-Datensatz, 200 je Lern-Datensatz, dgpfun.niter Lern-Datensätze
NA2res1[, "dgp"] <- 1
save(NA2res1, file = "07NA2res1.Rda")

# Regression:
#------------

source("funktionenZuMV2regres.R")
source("Simulationen2-alt-neu.R")
source("funktionenZuMV2regres.R")

packageDescription("mvtnorm")$Version
packageDescription("impute")$Version
packageDescription("party")$Version

#NA2res2 <- foo(RF = RF2, dgp = dgp2, ntest = 5000, dgpfun.niter = 5)
NA2res2 <- foo(RF = RF2, dgp = dgp2, ntest = 5000, dgpfun.niter = 100)
# ntest im Test-Datensatz, 200 je Lern-Datensatz, dgpfun.niter Lern-Datensätze
NA2res2[, "dgp"] <- 2
save(NA2res2, file = "07NA2res2.Rda")

# Datenaufbereitung:
#-------------------

NA2res <- rbind(NA2res1, NA2res2)
NA2res <- as.data.frame(NA2res)
NA2res$dgp <- factor(NA2res$dgp, labels = c("class", "regres"))
NA2res$sigma <- factor(NA2res$sigma, labels = paste("s", 1:3, sep = ""))
NA2res$MV <- factor(NA2res$MV, labels = c("none", "LOG", "DEPy"))
NA2res$lernMV <- as.logical(NA2res$lernMV)
NA2res$testMV <- as.logical(NA2res$testMV)
NA2res$imput <- as.logical(NA2res$imput)
NA2res$bench <- factor(NA2res$bench, labels = c("test", "benchmark"))
save(NA2res, file = "07NA2res.Rda")

################################################################################