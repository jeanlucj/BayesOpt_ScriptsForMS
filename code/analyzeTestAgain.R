require(here)
here::i_am("code/analyzeTestAgain.R")

# Analyze the test again stuff
# PIC budget, mean of separate optimiations
# 06 cycles: 0.2520365
# 12 cycles: 0.3874482
# Make plots for the budgets that are the mean of separate optimizations
sepOpt06 <- list.files(path=here::here("output", "PostGrid"),
                       pattern="0.3874482")
tstAgn06 <- list()
for (i in 1:length(sepOpt06)){
  tstAgn06 <- c(tstAgn06,
                readRDS(here::here("output", "PostGrid", sepOpt06[i])))
}

sepOpt12 <- list.files(path=here::here("output", "PostGrid"),
                       pattern="0.2520365")
tstAgn12 <- list()
for (i in 1:length(sepOpt12)){
  tstAgn12 <- c(tstAgn12,
                readRDS(here::here("output", "PostGrid", sepOpt12[i])))
}

# Put the list into an array: easier to manipulate
nRep <- length(tstAgn06)
sep_allParms06 <- array(dim=c(dim(tstAgn06[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  sep_allParms06[,,i] <- as.matrix(tstAgn06[[i]]$popParmsByCyc)
}
sep_popParmsByCyc06sd <- as.data.frame(apply(sep_allParms06, 1:2, sd))
colnames(sep_popParmsByCyc06sd) <- paste0(colnames(tstAgn06[[1]]$popParmsByCyc), "SD")
sep_popParmsByCyc06 <- as.data.frame(apply(sep_allParms06, 1:2, mean))
colnames(sep_popParmsByCyc06) <- colnames(tstAgn06[[1]]$popParmsByCyc)
sep_popParmsByCyc06$varCandMean <- sep_popParmsByCyc06$varCandMean -
  sep_popParmsByCyc06$breedPopMean[1]
sep_popParmsByCyc06$breedPopMean <- sep_popParmsByCyc06$breedPopMean -
  sep_popParmsByCyc06$breedPopMean[1]

nRep <- length(tstAgn12)
sep_allParms12 <- array(dim=c(dim(tstAgn12[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  sep_allParms12[,,i] <- as.matrix(tstAgn12[[i]]$popParmsByCyc)
}
sep_popParmsByCyc12sd <- as.data.frame(apply(sep_allParms12, 1:2, sd))
colnames(sep_popParmsByCyc12sd) <- paste0(colnames(tstAgn12[[1]]$popParmsByCyc), "SD")
sep_popParmsByCyc12 <- as.data.frame(apply(sep_allParms12, 1:2, mean))
colnames(sep_popParmsByCyc12) <- colnames(tstAgn12[[1]]$popParmsByCyc)
sep_popParmsByCyc12$varCandMean <- sep_popParmsByCyc12$varCandMean - sep_popParmsByCyc12$breedPopMean[1]
sep_popParmsByCyc12$breedPopMean <- sep_popParmsByCyc12$breedPopMean - sep_popParmsByCyc12$breedPopMean[1]

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
     pch=16, ylab="Diff varCandMean")
plot(sep_popParmsByCyc12$breedPopMean - sep_popParmsByCyc06$breedPopMean,
     pch=16, ylab="Diff breedPopMean")

# tstAgain_1_0.41512103.rds
# tstAgain_2_0.2523531.rds
# Make plots for the budgets that are from one joint optimum estimate
jntOpt06 <- list.files(path=here::here("output", "PostGrid"),
                       pattern="0.41512103")
tstAgn06 <- readRDS(here::here("output", "PostGrid", jntOpt06[1]))

jntOpt12 <- list.files(path=here::here("output", "PostGrid"),
                       pattern="0.2523531")
tstAgn12 <- readRDS(here::here("output", "PostGrid", jntOpt12[1]))

# Put the list into an array: easier to manipulate
nRep <- length(tstAgn06)
jnt_allParms06 <- array(dim=c(dim(tstAgn06[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  jnt_allParms06[,,i] <- as.matrix(tstAgn06[[i]]$popParmsByCyc)
}
jnt_popParmsByCyc06sd <- as.data.frame(apply(jnt_allParms06, 1:2, sd))
colnames(jnt_popParmsByCyc06sd) <- paste0(colnames(tstAgn06[[1]]$popParmsByCyc), "SD")
jnt_popParmsByCyc06 <- as.data.frame(apply(jnt_allParms06, 1:2, mean))
colnames(jnt_popParmsByCyc06) <- colnames(tstAgn06[[1]]$popParmsByCyc)
jnt_popParmsByCyc06$varCandMean <- jnt_popParmsByCyc06$varCandMean -
  jnt_popParmsByCyc06$breedPopMean[1]
jnt_popParmsByCyc06$breedPopMean <- jnt_popParmsByCyc06$breedPopMean -
  jnt_popParmsByCyc06$breedPopMean[1]

nRep <- length(tstAgn12)
jnt_allParms12 <- array(dim=c(dim(tstAgn12[[1]]$popParmsByCyc), nRep))
for (i in 1:nRep){
  jnt_allParms12[,,i] <- as.matrix(tstAgn12[[i]]$popParmsByCyc)
}
jnt_popParmsByCyc12sd <- as.data.frame(apply(jnt_allParms12, 1:2, sd))
colnames(jnt_popParmsByCyc12sd) <- paste0(colnames(tstAgn12[[1]]$popParmsByCyc), "SD")
jnt_popParmsByCyc12 <- as.data.frame(apply(jnt_allParms12, 1:2, mean))
colnames(jnt_popParmsByCyc12) <- colnames(tstAgn12[[1]]$popParmsByCyc)
jnt_popParmsByCyc12$varCandMean <- jnt_popParmsByCyc12$varCandMean - jnt_popParmsByCyc12$breedPopMean[1]
jnt_popParmsByCyc12$breedPopMean <- jnt_popParmsByCyc12$breedPopMean - jnt_popParmsByCyc12$breedPopMean[1]

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
     pch=16, ylab="Diff varCandMean")
plot(jnt_popParmsByCyc12$breedPopMean - jnt_popParmsByCyc06$breedPopMean,
     pch=16, ylab="Diff breedPopMean")
