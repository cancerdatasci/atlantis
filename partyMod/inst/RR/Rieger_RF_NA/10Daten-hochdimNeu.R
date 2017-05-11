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
#                                  Teil 9b:                                  #
#                        hochdimensionale Berechnungen                       #
#                       *** auf dem biostat-Server ***                       #
#                                                                            #
##############################################################################

source("SimulationenNeu.R")

packageDescription("mvtnorm")$Version
packageDescription("impute")$Version
packageDescription("party")$Version

hochdim1 <- foo(RF = RF1, dgp = dgp1, w = 45, ntest = 5000, dgpfun.niter = 100)
# ntest im Test-Datensatz, 200 je Lern-Datensatz, dgpfun.niter Lern-Datensätze
hochdim1[, "dgp"] <- 1
save(hochdim1, file = "10hochdim1neu.Rda")

source("SimulationenNeu.R")

packageDescription("mvtnorm")$Version
packageDescription("impute")$Version
packageDescription("party")$Version

hochdim2 <- foo(RF = RF2, dgp = dgp2, u = 45, ntest = 5000, dgpfun.niter = 100)
# ntest im Test-Datensatz, 200 je Lern-Datensatz, dgpfun.niter Lern-Datensätze
hochdim2[, "dgp"] <- 2
save(hochdim2, file = "10hochdim2neu.Rda")

hochdim <- rbind(hochdim1, hochdim2)
hochdim <- as.data.frame(hochdim)
hochdim$dgp <- factor(hochdim$dgp, labels = c("class", "regres"))
hochdim$sigma <- factor(hochdim$sigma, labels = paste("s", 1:3, sep = ""))
hochdim$MV <- factor(hochdim$MV, labels = c("none", paste("MAR", 1:4, sep = ""),"MCAR"))
hochdim$lernMV <- as.logical(hochdim$lernMV)
hochdim$testMV <- as.logical(hochdim$testMV)
hochdim$imput <- as.logical(hochdim$imput)
hochdim$bench <- factor(hochdim$bench, labels = c("test", "benchmark"))
save(hochdim, file = "10hochdimneu.Rda")

##############################################################################