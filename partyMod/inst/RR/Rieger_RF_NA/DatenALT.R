################################################################################
#                                                                              #
#                       R-Code zu meiner Bachelor-Arbeit                       #
#                                                                              #
#    "Random-Forest-Regression bei fehlenden Werten in den Einflussgrößen"     #
#                                                                              #
#                                   Teil 9:                                    #
#                                 Berechnungen                                 #
#                                                                              #
################################################################################

source("Simulationen.R")
res1 <- foo(RF = RF1, dgp = dgp1, ntest = 5000, dgpfun.niter = 100)
# 5000 im Test-Datensatz, 200 je Lern-Datensatz, 100 Lern-Datensätze
res1[, "dgp"] <- 1
save(res1, file = "res1.Rda")

source("Simulationen.R")
res2 <- foo(RF = RF2, dgp = dgp2, ntest = 5000, dgpfun.niter = 100)
# 5000 im Test-Datensatz, 200 je Lern-Datensatz, 100 Lern-Datensätze
res2[, "dgp"] <- 2
save(res2, file = "res2.Rda")

res <- rbind(res1, res2)
res <- as.data.frame(res)
res$dgp <- factor(res$dgp, labels = c("class", "regres"))
res$sigma <- factor(res$sigma, labels = paste("s", 1:3, sep = ""))
res$MV <- factor(res$MV, labels = c("none", paste("MAR", 1:4, sep = ""),"MCAR"))
res$lernMV <- as.logical(res$lernMV)
res$testMV <- as.logical(res$testMV)
res$imput <- as.logical(res$imput)
res$bench <- factor(res$bench, labels = c("test", "benchmark"))
save(res, file = "res.Rda")

################################################################################