rm(list = ls())
R.Version()$version.string

library("impute", lib.loc = "~/lib/oldR")
packageDescription("impute")$Version
library("mvtnorm")
packageDescription("mvtnorm")$Version
library("party", lib.loc = "~/lib/oldR")
packageDescription("party")$Version

##############################################################################
#                                                                            #
#                           R-Code zur HiWi-Stelle                           #
#                                                                            #
#                                  Teil 9d:                                  #
#      Berechnungen der zusätzlichen NA-Mechanismen im Hochdimensionalen     #
#                       *** auf dem biostat-Server ***                       #
#                                                                            #
##############################################################################

# Klassifikation:
#----------------

source("funktionenZuMV2class.R")
source("Simulationen2.R")
source("funktionenZuMV2class.R")

packageDescription("mvtnorm")$Version
packageDescription("impute")$Version
packageDescription("party")$Version

#NA2hd1 <- foo(RF = RF1, dgp = dgp1, w = 5, ntest = 5000, dgpfun.niter = 5)
NA2hd1 <- foo(RF = RF1, dgp = dgp1, w = 45, ntest = 5000, dgpfun.niter = 100)
# ntest im Test-Datensatz, 200 je Lern-Datensatz, dgpfun.niter Lern-Datensätze
NA2hd1[, "dgp"] <- 1
save(NA2hd1, file = "04NA2hochdim1.Rda")

# Regression:
#------------

source("funktionenZuMV2regres.R")
source("Simulationen2.R")
source("funktionenZuMV2regres.R")

packageDescription("mvtnorm")$Version
packageDescription("impute")$Version
packageDescription("party")$Version

#NA2hd2 <- foo(RF = RF2, dgp = dgp2, u = 10, ntest = 5000, dgpfun.niter = 5)
NA2hd2 <- foo(RF = RF2, dgp = dgp2, u = 45, ntest = 5000, dgpfun.niter = 100)
# ntest im Test-Datensatz, 200 je Lern-Datensatz, dgpfun.niter Lern-Datensätze
NA2hd2[, "dgp"] <- 2
save(NA2hd2, file = "04NA2hochdim2.Rda")

# Datenaufbereitung:
#-------------------

NA2hochdim <- rbind(NA2hd1, NA2hd2)
NA2hochdim <- as.data.frame(NA2hochdim)
NA2hochdim$dgp <- factor(NA2hochdim$dgp, labels = c("class", "regres"))
NA2hochdim$sigma <- factor(NA2hochdim$sigma, labels = paste("s", 1:3, sep = ""))
NA2hochdim$MV <- factor(NA2hochdim$MV, labels = c("none", "LOG", "DEPy"))
NA2hochdim$lernMV <- as.logical(NA2hochdim$lernMV)
NA2hochdim$testMV <- as.logical(NA2hochdim$testMV)
NA2hochdim$imput <- as.logical(NA2hochdim$imput)
NA2hochdim$bench <- factor(NA2hochdim$bench, labels = c("test", "benchmark"))
save(NA2hochdim, file = "04NA2hochdim.Rda")

################################################################################