# Analyze the test again stuff
# tstAgain0.2520365.rds
# tstAgain0.3874482.rds
require(here)
here::i_am("code/analyzeTestAgain.R")
tstAgn06 <- readRDS(here::here("output", "tstAgain0.3874482.rds"))
nRep <- length(tstAgn06)
allParms06 <- array(dim=c(dim(tstAgn06[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  allParms06[,,i] <- as.matrix(tstAgn06[[i]]$popParmsByCyc)
}
popParmsByCyc06sd <- as.data.frame(apply(allParms06, 1:2, sd))
colnames(popParmsByCyc06sd) <- paste0(colnames(tstAgn06[[1]]$popParmsByCyc), "SD")
popParmsByCyc06 <- as.data.frame(apply(allParms06, 1:2, mean))
colnames(popParmsByCyc06) <- colnames(tstAgn06[[1]]$popParmsByCyc)
popParmsByCyc06$varCandMean <- popParmsByCyc06$varCandMean -
  popParmsByCyc06$breedPopMean[1]
popParmsByCyc06$breedPopMean <- popParmsByCyc06$breedPopMean -
  popParmsByCyc06$breedPopMean[1]

tstAgn12 <- readRDS(here::here("output", "tstAgain0.2520365.rds"))
nRep <- length(tstAgn12)
allParms12 <- array(dim=c(dim(tstAgn12[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  allParms12[,,i] <- as.matrix(tstAgn12[[i]]$popParmsByCyc)
}
popParmsByCyc12sd <- as.data.frame(apply(allParms12, 1:2, sd))
colnames(popParmsByCyc12sd) <- paste0(colnames(tstAgn12[[1]]$popParmsByCyc), "SD")
popParmsByCyc12 <- as.data.frame(apply(allParms12, 1:2, mean))
colnames(popParmsByCyc12) <- colnames(tstAgn12[[1]]$popParmsByCyc)
popParmsByCyc12$varCandMean <- popParmsByCyc12$varCandMean - popParmsByCyc12$breedPopMean[1]
popParmsByCyc12$breedPopMean <- popParmsByCyc12$breedPopMean - popParmsByCyc12$breedPopMean[1]

# Plot the genetic gains
yRange <- range(c(popParmsByCyc06$varCandMean,
                  popParmsByCyc12$varCandMean))
pdf(here::here("output", "AllGain.pdf"))
plot(0:12, popParmsByCyc06$varCandMean,
     pch=16, ylim=yRange,
     xlab="Year", ylab="Gain of variety candidates",
     main="Impact of allocation on cumulative gain",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, popParmsByCyc12$varCandMean,
       pch=16, col="dark red",
       cex=1.3)
lines(c(6,6), c(0,60), lwd=0.5)
lines(c(12,12), c(0,60), lwd=0.5)
legend(7.15, 25, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

# Plot the additive genetic Std Dev.  I think on log scale
yRange <- range(c(log2(popParmsByCyc06$breedPopAddSD),
                  log2(popParmsByCyc12$breedPopAddSD)))
pdf(here::here("output", "AdditiveStdDev.pdf"))
plot(0:12, log2(popParmsByCyc06$breedPopAddSD),
     pch=16, ylim=yRange,
     xlab="Year", ylab="log2 breeding pop std. dev.",
     main="Impact of allocation on diversity in 12 years",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, log2(popParmsByCyc12$breedPopAddSD),
       pch=16, col="dark red",
       cex=1.3)
legend(8, 2.0, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

plot(popParmsByCyc12$varCandMean - popParmsByCyc06$varCandMean,
     pch=16, ylab="Diff varCandMean")
plot(popParmsByCyc12$breedPopMean - popParmsByCyc06$breedPopMean,
     pch=16, ylab="Diff breedPopMean")
