##############################################################################
#                                                                            #
#                           R-Code zur HiWi-Stelle                           #
#                                                                            #
#                  Berechnungen am realen Datensatz "Ozone"                  #
#                                                                            #
##############################################################################



library(mlbench)
data(Ozone)

# fünf fehlende Werte im Response V4 -> rausschmeißen
Ozone <- Ozone[-which(is.na(Ozone$V4)), ]

library("impute", lib.loc = "~/lib/oldR")
packageDescription("impute")$Version
#library("mvtnorm")
#packageDescription("mvtnorm")$Version
library("party", lib.loc = "~/lib/oldR")
packageDescription("party")$Version

nIter <- 100
MSE <- impMSE <- vector("numeric")

set.seed(856329)

for (i in 1:nIter) {
    rf <- try(cforest(V4 ~ ., data = dat[1,],
        control = cforest_control(maxsurrogate=3, ntree=50, minsplit=30)), silent=TRUE)
    while(inherits(rf, "try-error")) {

        bootstrap <- sample(1:dim(Ozone)[1], replace=TRUE)
        dat <- Ozone[bootstrap, ]   # genauso groß wie Original-Stpr.
        MSE[i] <- impMSE[i] <- NA

        # RF + Surrogat
        #--------------
        cat("Surrogates ")
    
        rf <- try(cforest(V4 ~ ., data = dat,
            control = cforest_control(maxsurrogate = 3, ntree = 50, minsplit = 30)))
    }
    print(i)


    oob <- Ozone[-bootstrap, ]

    f <- predict(rf, newdata = oob)
    MSE[i] <- mean((oob$V4 - f)^2)

    # Imputation + RF
    #----------------
    cat("Imputation ")

    dat[-1] <- as.data.frame(impute.knn(as.matrix(dat[-1]))$data)
    for(j in 1:3) {
        dat[, paste("V", j, sep="")] <- as.factor(dat[, paste("V", j, sep="")])
    }
    for(j in 4:13) {
        dat[, paste("V", j, sep="")] <- as.numeric(dat[, paste("V", j, sep="")])
    }

    rf <- cforest(V4 ~ ., data = dat,
        control = cforest_control(maxsurrogate = 3, ntree = 50, minsplit = 30))
    print(i)

    f <- predict(rf, newdata = oob)
    impMSE[i] <- mean((oob$V4 - f)^2)
}

boxplot(MSE - impMSE, main = "Differences between RF + Surrogate and Imputation + RF")
abline(h=0)

boxplot(MSE, impMSE, MSE-impMSE, names=c("RF + surrogate", "Imputation + RF",
    "differences"))
abline(h=0)

#------------------------------------------------------------------------------#
#save(MSE, impMSE, file = "Ozone.Rdata")
load("Ozone.Rdata")
#------------------------------------------------------------------------------#

MSE <- as.data.frame(MSE)
colnames(MSE) <- "risk"
ozon <- cbind(MSE, imp = FALSE)
head(ozon)
summary(ozon)

impMSE <- as.data.frame(impMSE)
colnames(impMSE) <- "risk"
ozon2 <- cbind(impMSE, imp = TRUE)
head(ozon2)
summary(ozon2)

ozon3 <- rbind(ozon, ozon2)
head(ozon3)
summary(ozon3)

ozon <- ozon3

#------------------------------------------------------------------------------#
#save(ozon, file = "ozon.Rdata")
load("ozon.Rdata")
#------------------------------------------------------------------------------#

library(lattice)

ozon$imput <- as.factor(ozon$imp)
levels(ozon$imput) <- c("sur", "knn")

pdf(file = "Ozone.pdf", width = 5, height = 5)
bwplot(risk ~ imput, data = ozon, cex.lab = 1.5, cex.axis = 1.5, font = 2)
dev.off()

t.test(risk ~ imput, data = ozon, paired = FALSE)

################################################################################