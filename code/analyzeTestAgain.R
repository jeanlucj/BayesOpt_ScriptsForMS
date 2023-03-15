require(here)
library(tidyverse)
here::i_am("code/analyzeTestAgain.R")

# Get the testAgain files for a particular PIC budget
getTestAgain <- function(picBudgString="0.3874482"){
  rdsFiles <- list.files(path=here::here("output", "testAgain"),
                         pattern=picBudgString)
  allList <- list()
  for (i in 1:length(rdsFiles)){
    allList <- c(allList,
                 readRDS(here::here("output", "testAgain", rdsFiles[i])))
  }

  # Put the list into an array: easier to manipulate
  nRep <- length(allList)
  allParms <- array(dim=c(dim(allList[[1]]$popParmsByCyc), nRep))
  for (i in 1:nRep){
    ppc <- allList[[i]]$popParmsByCyc
    ppc$varCandMean <- ppc$varCandMean - ppc$breedPopMean[1]
    ppc$breedPopMean <- ppc$breedPopMean - ppc$breedPopMean[1]
    allParms[,,i] <- as.matrix(ppc)
  }
  parmsMean <- as.data.frame(apply(allParms, 1:2, mean))
  colnames(parmsMean) <- colnames(allList[[1]]$popParmsByCyc)
  parmsSE <- as.data.frame(apply(allParms, 1:2, sd) / sqrt(nRep))
  colnames(parmsSE) <- paste0(colnames(allList[[1]]$popParmsByCyc), "se")
  return(list(parmsMean=parmsMean, parmsSE=parmsSE))
}
# Analyzing small test again
smlTA <- getTestAgain(picBudgString="0.36805244")
View(smlTA$parmsMean)
View(smlTA$parmsSE)
varCandMeanSD <- dplyr::last(smlTA$parmsSE$varCandMeanse)*20 # 5.107004

# Lots of assembling of data happens throughout
# The final plots that I want to do are at the bottom, starting around line 210
# Analyze the test again stuff
# PIC budget, mean of separate optimiations
# 06 cycles: 0.2520365
# 12 cycles: 0.3874482
# Make plots for the budgets that are the mean of separate optimizations

sep_popParmsByCyc06 <- getTestAgain(picBudgString="0.3874482")
sep_popParmsByCyc06se <- sep_popParmsByCyc06[[2]]
sep_popParmsByCyc06 <- sep_popParmsByCyc06[[1]]
sep_popParmsByCyc12 <- getTestAgain(picBudgString="0.2520365")
sep_popParmsByCyc12se <- sep_popParmsByCyc12[[2]]
sep_popParmsByCyc12 <- sep_popParmsByCyc12[[1]]

sp8_popParmsByCyc06 <- getTestAgain(picBudgString="0.3858742")
sp8_popParmsByCyc06se <- sp8_popParmsByCyc06[[2]]
sp8_popParmsByCyc06 <- sp8_popParmsByCyc06[[1]]
sp8_popParmsByCyc12 <- getTestAgain(picBudgString="0.2620301")
sp8_popParmsByCyc12se <- sp8_popParmsByCyc12[[2]]
sp8_popParmsByCyc12 <- sp8_popParmsByCyc12[[1]]

jnt_popParmsByCyc06 <- getTestAgain(picBudgString="0.41512103")
jnt_popParmsByCyc06se <- jnt_popParmsByCyc06[[2]]
jnt_popParmsByCyc06 <- jnt_popParmsByCyc06[[1]]
jnt_popParmsByCyc12 <- getTestAgain(picBudgString="0.2523531")
jnt_popParmsByCyc12se <- jnt_popParmsByCyc12[[2]]
jnt_popParmsByCyc12 <- jnt_popParmsByCyc12[[1]]

pc1_popParmsByCyc06 <- getTestAgain(picBudgString="0.3235522")
pc1_popParmsByCyc06se <- pc1_popParmsByCyc06[[2]]
pc1_popParmsByCyc06 <- pc1_popParmsByCyc06[[1]]
pc1_popParmsByCyc12 <- getTestAgain(picBudgString="0.3836027")
pc1_popParmsByCyc12se <- pc1_popParmsByCyc12[[2]]
pc1_popParmsByCyc12 <- pc1_popParmsByCyc12[[1]]

allVarCand <- tibble(year=0:12,
                     sep6=sep_popParmsByCyc06$varCandMean,
                     sep12=sep_popParmsByCyc12$varCandMean,
                     sp86=sp8_popParmsByCyc06$varCandMean,
                     sp812=sp8_popParmsByCyc12$varCandMean,
                     jnt6=jnt_popParmsByCyc06$varCandMean,
                     jnt12=jnt_popParmsByCyc12$varCandMean,
                     pc16=pc1_popParmsByCyc06$varCandMean,
                     pc112=pc1_popParmsByCyc12$varCandMean)
avc <- allVarCand %>% filter(year %in% c(6, 12))

# Plot the genetic gains
yRange <- range(c(sep_popParmsByCyc06$varCandMean,
                  sep_popParmsByCyc12$varCandMean))
pdf(here::here("output", "AllGain.pdf"))
plot(0:12, sep_popParmsByCyc06$varCandMean,
     pch=16, ylim=yRange,
     xlab="Year", ylab="Gain of variety candidates",
     main="Impact of allocation on cumulative gain",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, sep_popParmsByCyc12$varCandMean,
       pch=16, col="dark red",
       cex=1.3)
lines(c(6,6), c(0,60), lwd=0.5)
lines(c(12,12), c(0,60), lwd=0.5)
legend(7.15, 25, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

# Plot the additive genetic Std Dev.  I think on log scale
yRange <- range(c(log2(sep_popParmsByCyc06$breedPopAddSD),
                  log2(sep_popParmsByCyc12$breedPopAddSD)))
pdf(here::here("output", "AdditiveStdDev.pdf"))
plot(0:12, log2(sep_popParmsByCyc06$breedPopAddSD),
     pch=16, ylim=yRange,
     xlab="Year", ylab="log2 breeding pop std. dev.",
     main="Impact of allocation on diversity in 12 years",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, log2(sep_popParmsByCyc12$breedPopAddSD),
       pch=16, col="dark red",
       cex=1.3)
legend(8, 2.0, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

plot(sep_popParmsByCyc12$varCandMean - sep_popParmsByCyc06$varCandMean,
     pch=16, ylab="Diff varCandMean", main="Separate Budget")
plot(sep_popParmsByCyc12$breedPopMean - sep_popParmsByCyc06$breedPopMean,
     pch=16, ylab="Diff breedPopMean", main="Separate Budget")

# tstAgain_1_0.41512103.rds
# tstAgain_2_0.2523531.rds
# Make plots for the budgets that are from one joint optimum estimate
jntOpt06 <- list.files(path=here::here("output", "testAgain"),
                       pattern="0.41512103")
tstAgn06 <- list()
for (i in 1:length(jntOpt06)){
  tstAgn06 <- c(tstAgn06,
                readRDS(here::here("output", "testAgain", jntOpt06[i])))
}

jntOpt12 <- list.files(path=here::here("output", "testAgain"),
                       pattern="0.2523531")
tstAgn12 <- list()
for (i in 1:length(jntOpt12)){
  tstAgn12 <- c(tstAgn12,
                readRDS(here::here("output", "testAgain", jntOpt12[i])))
}

# Put the list into an array: easier to manipulate
nRep <- length(tstAgn06)
jnt_allParms06 <- array(dim=c(dim(tstAgn06[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  ppc <- tstAgn06[[i]]$popParmsByCyc
  ppc$varCandMean <- ppc$varCandMean - ppc$breedPopMean[1]
  ppc$breedPopMean <- ppc$breedPopMean - ppc$breedPopMean[1]
  jnt_allParms06[,,i] <- as.matrix(ppc)
}

jnt_popParmsByCyc06se <- as.data.frame(apply(jnt_allParms06, 1:2, sd) / sqrt(nRep))
colnames(jnt_popParmsByCyc06se) <- paste0(colnames(tstAgn06[[1]]$popParmsByCyc), "se")
jnt_popParmsByCyc06 <- as.data.frame(apply(jnt_allParms06, 1:2, mean))
colnames(jnt_popParmsByCyc06) <- colnames(tstAgn06[[1]]$popParmsByCyc)

nRep <- length(tstAgn12)
jnt_allParms12 <- array(dim=c(dim(tstAgn12[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  ppc <- tstAgn12[[i]]$popParmsByCyc
  ppc$varCandMean <- ppc$varCandMean - ppc$breedPopMean[1]
  ppc$breedPopMean <- ppc$breedPopMean - ppc$breedPopMean[1]
  jnt_allParms12[,,i] <- as.matrix(ppc)
}
jnt_popParmsByCyc12se <- as.data.frame(apply(jnt_allParms12, 1:2, sd) / sqrt(nRep))
colnames(jnt_popParmsByCyc12se) <- paste0(colnames(tstAgn12[[1]]$popParmsByCyc), "se")
jnt_popParmsByCyc12 <- as.data.frame(apply(jnt_allParms12, 1:2, mean))
colnames(jnt_popParmsByCyc12) <- colnames(tstAgn12[[1]]$popParmsByCyc)

# Plot the genetic gains
yRange <- range(c(jnt_popParmsByCyc06$varCandMean,
                  jnt_popParmsByCyc12$varCandMean))
pdf(here::here("output", "Jnt_AllGain.pdf"))
plot(0:12, jnt_popParmsByCyc06$varCandMean,
     pch=16, ylim=yRange,
     xlab="Year", ylab="Gain of variety candidates",
     main="Impact of allocation on cumulative gain",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, jnt_popParmsByCyc12$varCandMean,
       pch=16, col="dark red",
       cex=1.3)
lines(c(6,6), c(0,60), lwd=0.5)
lines(c(12,12), c(0,60), lwd=0.5)
legend(7.15, 25, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

# Plot the additive genetic Std Dev.  I think on log scale
yRange <- range(c(log2(jnt_popParmsByCyc06$breedPopAddSD),
                  log2(jnt_popParmsByCyc12$breedPopAddSD)))
pdf(here::here("output", "Jnt_AdditiveStdDev.pdf"))
plot(0:12, log2(jnt_popParmsByCyc06$breedPopAddSD),
     pch=16, ylim=yRange,
     xlab="Year", ylab="log2 breeding pop std. dev.",
     main="Impact of allocation on diversity in 12 years",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, log2(jnt_popParmsByCyc12$breedPopAddSD),
       pch=16, col="dark red",
       cex=1.3)
legend(8, 2.0, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

plot(jnt_popParmsByCyc12$varCandMean - jnt_popParmsByCyc06$varCandMean,
     pch=16, ylab="Diff varCandMean", main="Joint Budget")
plot(jnt_popParmsByCyc12$breedPopMean - jnt_popParmsByCyc06$breedPopMean,
     pch=16, ylab="Diff breedPopMean", main="Joint Budget")

# The test of whether joint or separate is better:
# Six cycle optimum at six cycles
cat("Separate", sep_popParmsByCyc06$varCandMean[7], "\n")
cat("Joint   ", jnt_popParmsByCyc06$varCandMean[7], "\n")
# Twelve cycle optimum at twelve cycles
cat("Separate", sep_popParmsByCyc12$varCandMean[13], "\n")
cat("Joint   ", jnt_popParmsByCyc12$varCandMean[13], "\n")

# At the wrong time point
# Six cycle optimum at twelve cycles
cat("Separate", sep_popParmsByCyc06$varCandMean[13], "\n")
cat("Joint   ", jnt_popParmsByCyc06$varCandMean[13], "\n")
# Twelve cycle optimum at six cycles
cat("Separate", sep_popParmsByCyc12$varCandMean[7], "\n")
cat("Joint   ", jnt_popParmsByCyc12$varCandMean[7], "\n")

# Final plotting
# Plot the genetic gains
yRange <- range(c(sep_popParmsByCyc06$varCandMean,
                  jnt_popParmsByCyc12$varCandMean))
pdf(here::here("output", "Hyb_AllGain.pdf"))
plot(0:12, sep_popParmsByCyc06$varCandMean,
     pch=16, ylim=yRange,
     xlab="Year", ylab="Gain of variety candidates",
     main="Impact of allocation on cumulative gain",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, jnt_popParmsByCyc12$varCandMean,
       pch=16, col="dark red",
       cex=1.3)
lines(c(6,6), yRange, lwd=0.5)
lines(c(12,12), yRange, lwd=0.5)
# This shows that the std err is smaller than half the width of the symbol
#lines(c(11.5,12.5),
#      rep(jnt_popParmsByCyc12$varCandMean[13] -
#            jnt_popParmsByCyc12se$varCandMeanse[13], 2), lwd=0.5)
#lines(c(11.5, 12.5),
#      rep(jnt_popParmsByCyc12$varCandMean[13] +
#            jnt_popParmsByCyc12se$varCandMeanse[13], 2), lwd=0.5)
legend(7.15, 25, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

# Plot the additive genetic Std Dev.  I think on log scale
yRange <- range(c(log2(sep_popParmsByCyc06$breedPopAddSD),
                  log2(jnt_popParmsByCyc12$breedPopAddSD)))
pdf(here::here("output", "Hyb_AdditiveStdDev.pdf"))
plot(0:12, log2(sep_popParmsByCyc06$breedPopAddSD),
     pch=16, ylim=yRange,
     xlab="Year", ylab="log2 breeding pop std. dev.",
     main="Impact of allocation on diversity in 12 years",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, log2(jnt_popParmsByCyc12$breedPopAddSD),
       pch=16, col="dark red",
       cex=1.3)
# These lines are not really needed here...
#lines(c(6,6), yRange, lwd=0.5)
#lines(c(12,12), yRange, lwd=0.5)
# This shows that the std err is smaller than half the width of the symbol
#lines(c(5.5,6.5),
#      rep(log2(jnt_popParmsByCyc12$breedPopAddSD[7] -
#            jnt_popParmsByCyc12se$breedPopAddSDse[7]), 2), lwd=0.5)
#lines(c(5.5, 6.5),
#      rep(log2(jnt_popParmsByCyc12$breedPopAddSD[7] +
#            jnt_popParmsByCyc12se$breedPopAddSDse[7]), 2), lwd=0.5)
legend(8, 2.0, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()


pic1Opt06 <- list.files(path=here::here("output", "testAgain"),
                        pattern="0.32355")
tstAgn06 <- list()
for (i in 1:length(pic1Opt06)){
  tstAgn06 <- c(tstAgn06,
                readRDS(here::here("output", "testAgain", pic1Opt06[i])))
}

pic1Opt12 <- list.files(path=here::here("output", "testAgain"),
                        pattern="0.38360")
tstAgn12 <- list()
for (i in 1:length(pic1Opt12)){
  tstAgn12 <- c(tstAgn12,
                readRDS(here::here("output", "testAgain", pic1Opt12[i])))
}

# Put the list into an array: easier to manipulate
nRep <- length(tstAgn06)
pc1_allParms06 <- array(dim=c(dim(tstAgn06[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  ppc <- tstAgn06[[i]]$popParmsByCyc
  ppc$varCandMean <- ppc$varCandMean - ppc$breedPopMean[1]
  ppc$breedPopMean <- ppc$breedPopMean - ppc$breedPopMean[1]
  pc1_allParms06[,,i] <- as.matrix(ppc)
}

pc1_popParmsByCyc06se <- as.data.frame(apply(pc1_allParms06, 1:2, sd) / sqrt(nRep))
colnames(pc1_popParmsByCyc06se) <- paste0(colnames(tstAgn06[[1]]$popParmsByCyc), "se")
pc1_popParmsByCyc06 <- as.data.frame(apply(pc1_allParms06, 1:2, mean))
colnames(pc1_popParmsByCyc06) <- colnames(tstAgn06[[1]]$popParmsByCyc)

nRep <- length(tstAgn12)
pc1_allParms12 <- array(dim=c(dim(tstAgn12[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  ppc <- tstAgn12[[i]]$popParmsByCyc
  ppc$varCandMean <- ppc$varCandMean - ppc$breedPopMean[1]
  ppc$breedPopMean <- ppc$breedPopMean - ppc$breedPopMean[1]
  pc1_allParms12[,,i] <- as.matrix(ppc)
}

pc1_popParmsByCyc12se <- as.data.frame(apply(pc1_allParms12, 1:2, sd) / sqrt(nRep))
colnames(pc1_popParmsByCyc12se) <- paste0(colnames(tstAgn12[[1]]$popParmsByCyc), "se")
pc1_popParmsByCyc12 <- as.data.frame(apply(pc1_allParms12, 1:2, mean))
colnames(pc1_popParmsByCyc12) <- colnames(tstAgn12[[1]]$popParmsByCyc)

yRange <- range(c(sep_popParmsByCyc06$varCandMean,
                  jnt_popParmsByCyc12$varCandMean,
                  pc1_popParmsByCyc06$varCandMean,
                  pc1_popParmsByCyc12$varCandMean))
pdf(here::here("output", "Hyb+PC1_AllGain.pdf"))
plot(0:12, sep_popParmsByCyc06$varCandMean,
     pch=16, ylim=yRange,
     xlab="Year", ylab="Gain of variety candidates",
     main="Impact of allocation on cumulative gain",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, jnt_popParmsByCyc12$varCandMean,
       pch=16, col="dark red",
       cex=1.3)
points(0:12, pc1_popParmsByCyc06$varCandMean,
       pch=18, cex=1.3)
points(0:12, pc1_popParmsByCyc12$varCandMean,
       pch=18, col="dark red",
       cex=1.3)
lines(c(6,6), yRange, lwd=0.5)
lines(c(12,12), yRange, lwd=0.5)
dev.off()


# Plot the additive genetic Std Dev.  I think on log scale
yRange <- range(c(log2(sep_popParmsByCyc06$breedPopAddSD),
                  log2(jnt_popParmsByCyc12$breedPopAddSD),
                  log2(pc1_popParmsByCyc06$breedPopAddSD),
                  log2(pc1_popParmsByCyc12$breedPopAddSD)))
pdf(here::here("output", "Hyb+PC1_AdditiveStdDev.pdf"))
plot(0:12, log2(sep_popParmsByCyc06$breedPopAddSD),
     pch=16, ylim=yRange,
     xlab="Year", ylab="log2 breeding pop std. dev.",
     main="Impact of allocation on diversity in 12 years",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(0:12, log2(jnt_popParmsByCyc12$breedPopAddSD),
       pch=16, col="dark red",
       cex=1.3)
points(0:12, log2(pc1_popParmsByCyc06$breedPopAddSD),
       pch=18, cex=1.3)
points(0:12, log2(pc1_popParmsByCyc12$breedPopAddSD),
       pch=18, col="dark red",
       cex=1.3)
# These lines are not really needed here...
#lines(c(6,6), yRange, lwd=0.5)
#lines(c(12,12), yRange, lwd=0.5)
# This shows that the std err is smaller than half the width of the symbol
#lines(c(5.5,6.5),
#      rep(log2(jnt_popParmsByCyc12$breedPopAddSD[7] -
#            jnt_popParmsByCyc12se$breedPopAddSDse[7]), 2), lwd=0.5)
#lines(c(5.5, 6.5),
#      rep(log2(jnt_popParmsByCyc12$breedPopAddSD[7] +
#            jnt_popParmsByCyc12se$breedPopAddSDse[7]), 2), lwd=0.5)
legend(8, 2.0, c("12-Yr Opt", " 6-Yr Opt"), pch=16,
       col=c("dark red", "black"),
       cex=1.5)
dev.off()

################################################################
# Do some selectiongain VDP analyses to see how the VDP performs for each
# of these budgets
# Extract the genetic variance at the start of genomic selection
pc2Cy06 <- sep_popParmsByCyc06$breedPopSD[1]
pc2Cy12 <- jnt_popParmsByCyc12$breedPopSD[1]
pc1Cy06 <- pc1_popParmsByCyc06$breedPopSD[1]
pc1Cy12 <- pc1_popParmsByCyc12$breedPopSD[1]
startVar <- mean(c(pc2Cy06, pc2Cy12, pc1Cy06, pc1Cy12))^2 # 26.816

if (exists("init_num")) rm(init_num)
source(here::here("code", "breedSchemeInitialize.R"))
print(bsd$initBudget)
percentages <- bsd$initBudget[grep("perc_", names(bsd$initBudget))]
# tst <- BreedSimCost::budgetToScheme(percentages=percentages, bsd)
bsd$genVar <- startVar
vdpCovMat <- BreedSimCost::calcVDPcovMat(bsd)
corMat <- cov2cor(vdpCovMat)
selFrac <- c(bsd$nEntries[-1], bsd$nToMarketingDept) / bsd$nEntries
truncPts <- BreedSimCost::multistageTruncPt(alpha=selFrac, corr=corMat)
gain <- BreedSimCost::multistageGain(corMat, truncPts)
cat("Burn-in Gain", gain, "\n")
# Burn-in Gain 1.8761

# Best budgets
# At the moment in the MS
#  6 1 0.324 (0.033) 0.416 (0.025) 0.236 (0.031) 0.024 (0.005) # OK
#  6 2 0.386 (0.017) 0.423 (0.028) 0.139 (0.017) 0.052 (0.013) # different
# 12 1 0.384 (0.038) 0.418 (0.029) 0.179 (0.036) 0.019 (0.002) # OK
# 12 2 0.252 (0.021) 0.342 (0.020) 0.384 (0.036) 0.022 (0.010) # OK
# Put this one in the manuscript
# budg6Cyc <- c(0.3874482, 0.3870322, 0.1624069, 0.06311256) # Used

# No change on which used as a result of doing 8 reps
budgAT6Cyc <- c(0.41512103, 0.41056531, 0.15239730, 0.02191636)
budgAT12Cyc <- c(0.2523531, 0.3416312, 0.38409934, 0.02191636) # Used
budg6Cyc <- c(0.3874482, 0.3870322, 0.1624069, 0.06311256) # Used
budg12Cyc <- c(0.2520365, 0.3929344, 0.2805345, 0.07449453)
budg6Cyc8 <- c(0.3858742, 0.4231775, 0.1392178, 0.05173058) # 8 reps
budg12Cyc8 <- c(0.2620301, 0.3966589, 0.2810439, 0.06026706) # 8 reps
budgP1Cyc6 <- c(0.3235522, 0.4163751, 0.2357820, 0.02429063) # Used
budgP1Cyc12 <- c(0.3836027, 0.4179224, 0.1791765, 0.01929836) # Used
allBudg <- list(budgAT6Cyc, budgAT12Cyc,
                budg6Cyc, budg12Cyc,
                budgP1Cyc6, budgP1Cyc12,
                budg6Cyc8, budg12Cyc8)

# Figure out the selectiongain predicted gain from the VDP
for (whichBudg in c(3, 2, 5, 6)){
  perc <- allBudg[[whichBudg]]
  tst <- BreedSimCost::budgetToScheme(percentages=perc, bsd)
  print(tst$nEntries)
  selFrac <- c(tst$nEntries[-1], tst$nToMarketingDept) / tst$nEntries
  truncPts <- BreedSimCost::multistageTruncPt(alpha=selFrac, corr=corMat)
  gain <- BreedSimCost::multistageGain(corMat, truncPts)
  cat("Which budget:", whichBudg, " Gain:", gain, "\n")
}
# SDN CET PYT Cyc06 PIC2
# 864 162  14
# Which budget: 3  Gain: 1.850112
# Which budget: 3  Maximum gain: 1.894315 [97.7%]
# Budget at max: 0.3874482 0.3166099 0.1785887 0.1173532
# SDN CET PYT
# 707 178  27

# SDN CET PYT Cyc12 PIC2
# 762 382   5
# Which budget: 2  Gain: 1.607731
# Which budget: 2  Maximum gain: 1.950499 [82.4%]
# Budget at max: 0.2523531 0.3794874 0.2333103 0.1348491
# SDN CET PYT
# 847 232  31

# SDN CET PYT Cyc06 PIC1
# 929 235   6
# Which budget: 5  Gain: 1.699636
# Which budget: 5  Maximum gain: 1.921084 [88.5%]
# Budget at max: 0.3235522 0.3523052 0.1940094 0.1301332
# SDN CET PYT
# 786 193  30

# SDN CET PYT Cyc12 PIC1
# 933 178   4
# Which budget: 6  Gain: 1.469301
# Which budget: 6  Maximum gain: 1.896636 [77.5%]
# Budget at max: 0.3836027 0.3124614 0.1909625 0.1129735
# SDN CET PYT
# 697 190  26

# For each budget, figure out what would be the optimum for the VDP using a grid
# search with selectiongain
# Make search grids
source(here::here("code", "postGridFunctions.R"))
gridBudgets <- makeGridForTernary()
readr::write_csv(as.data.frame(gridBudgets),
                 file=here::here("output", "GridFullStandard.csv"),
                 col_names=F, quote="none")
# NOTE: before doing the following, in the OptimizationCtrlFile.txt make
# sure you have (without the #)
# minPercentage
# 0.22 0.33 0.10 0.02
# maxPercentage
# 0.42 0.43 0.40 0.08
gridOptBudgets <- makeGridForTernary(doMax=F)
readr::write_csv(as.data.frame(gridOptBudgets),
                 file=here::here("output", "GridCloseOptStandard.csv"),
                 col_names=F, quote="none")
# Make grids for each optimum (6 vs 12) x (PIC1 vs PIC2)
nCycTxt <- c("06", "12")
nPICtxt <- c("PC2", "PC1")
nameNum <- 0
for (whichBudg in c(3, 2, 5, 6)){
  budgName <- paste0(nCycTxt[nameNum %% 2 + 1], nPICtxt[nameNum %/% 2 + 1])
  perc <- allBudg[[whichBudg]]
  picBudget <- perc[1]
  gridBudgets <- makeGridForTernary(picBudget=picBudget, justVDP=T, doMax=T)
  readr::write_csv(as.data.frame(gridBudgets),
                   file=here::here("output", paste0("GridFull", budgName, ".csv")),
                   col_names=F, quote="none")
  # minPercentage
  # 0.22 0.23 0.10 0.02
  # maxPercentage
  # 0.42 0.53 0.50 0.22
  gridOptBudgets <- makeGridForTernary(picBudget=picBudget, justVDP=T, doMax=F)
  readr::write_csv(as.data.frame(gridOptBudgets),
                   file=here::here("output", paste0("GridCloseOpt", budgName, ".csv")),
                   col_names=F, quote="none")
  nameNum <- nameNum+1
}

# Now run the selectiongain search grid on these budgets...
nameNum <- 0
for (whichBudg in c(3, 2, 5, 6)){
  budgName <- paste0(nCycTxt[nameNum %% 2 + 1], nPICtxt[nameNum %/% 2 + 1])
  gridBudgets <- readr::read_csv(
    file=here::here("output", paste0("GridFull", budgName, ".csv")),
    col_names=F)

  allGain <- NULL
  for (i in 1:nrow(gridBudgets)){
    perc <- gridBudgets[i,] %>% unlist
    perc <- c(perc, 1-sum(perc))
    tst <- BreedSimCost::budgetToScheme(percentages=perc, bsd)
    selFrac <- c(tst$nEntries[-1], tst$nToMarketingDept) / tst$nEntries
    truncPts <- BreedSimCost::multistageTruncPt(alpha=selFrac, corr=corMat)
    gain <- BreedSimCost::multistageGain(corMat, truncPts)
    allGain <- c(allGain, gain)
  }
  cat("Which budget:", whichBudg, " Maximum gain:", max(allGain), "\n")
  perc <- gridBudgets[which.max(allGain),] %>% unlist
  perc <- c(perc, 1-sum(perc))
  tst <- BreedSimCost::budgetToScheme(percentages=perc, bsd)
  cat("Budget at max:", perc, "\n")
  print(tst$nEntries)

  gridBudgets <- readr::read_csv(
    file=here::here("output", paste0("GridCloseOpt", budgName, ".csv")),
    col_names=F)

  allGain <- NULL
  for (i in 1:nrow(gridBudgets)){
    perc <- gridBudgets[i,] %>% unlist
    perc <- c(perc, 1-sum(perc))
    tst <- BreedSimCost::budgetToScheme(percentages=perc, bsd)
    selFrac <- c(tst$nEntries[-1], tst$nToMarketingDept) / tst$nEntries
    truncPts <- BreedSimCost::multistageTruncPt(alpha=selFrac, corr=corMat)
    gain <- BreedSimCost::multistageGain(corMat, truncPts)
    allGain <- c(allGain, gain)
  }
  cat("Which budget:", whichBudg, " Maximum gain:", max(allGain), "\n")
  perc <- gridBudgets[which.max(allGain),] %>% unlist
  perc <- c(perc, 1-sum(perc))
  tst <- BreedSimCost::budgetToScheme(percentages=perc, bsd)
  cat("Budget at max:", perc, "\n")
  print(tst$nEntries)

  nameNum <- nameNum+1
}

# Which budget: 3  Maximum gain: 1.894315
# Which budget: 2  Maximum gain: 1.950499
# Which budget: 5  Maximum gain: 1.921084
# Which budget: 6  Maximum gain: 1.896636
