# Uh oh.  The two replicates of smlShrtMany are giving exactly the same thing
# I will have to check if somehow I set the seed to the same thing twice
# or if there was a duplication downstream.
library(tidyverse)
here::i_am("code/analyzeSmlShrtMany.R")
setwd("/Users/jj332/Documents/simOut/bso/output/PostGrid")
allSml <- list.files(pattern="sml")
# Parse allSml of this form 188690362_0_1036sml06PIC.rds
where_ <- regexpr("_", allSml)
wheres <- regexpr("s", allSml)
subSeed <- substring(allSml, 1, where_ - 1)
optimRepl <- substring(allSml, where_ + 1, where_ + 1)
nIter <- as.numeric(substring(allSml, where_ + 3, wheres - 1))

allRes <- NULL
for (f in allSml){
  res <- readRDS(f)
  allRes <- rbind(allRes, c(res$bestBudget, res$maxPredGain, res$postVarAtMax))
}
colnames(allRes) <- c("budgPIC", "budgSDN", "budgCET",
                      "maxPredGain", "postVarAtMax")
allSml <- tibble(fileName=allSml, subSeed=subSeed,
                 optimRepl=optimRepl, nOptIter=nIter)
allSml <- dplyr::bind_cols(allSml, allRes)
allIter <- allSml %>% pull(nOptIter) %>% unique %>% sort

fitSubSeed <- lm(budgPIC ~ subSeed, data=allSml, subset=(nOptIter==1536))
anova(fitSubSeed)
fitSubSeed <- lm(budgSDN ~ subSeed, data=allSml, subset=(nOptIter==1536))
anova(fitSubSeed)
fitSubSeed <- lm(budgCET ~ subSeed, data=allSml, subset=(nOptIter==1536))
anova(fitSubSeed)
# Response: budgPIC
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# subSeed   11 0.022027 0.0020025  0.6218 0.7802
# Residuals 12 0.038646 0.0032205
# Response: budgSDN
#           Df   Sum Sq   Mean Sq F value    Pr(>F)
# subSeed   11 0.176479 0.0160435   29.25 5.926e-07 ***
# Residuals 12 0.006582 0.0005485
# Response: budgCET
#           Df   Sum Sq   Mean Sq F value  Pr(>F)
# subSeed   11 0.126047 0.0114588  3.6893 0.01692 *
# Residuals 12 0.037271 0.0031059

man <- manova(cbind(budgPIC, budgSDN, budgCET) ~ subSeed,
              data=allSml, subset=(nOptIter==1536))
summary(man)
#           Df Pillai approx F num Df den Df  Pr(>F)
# subSeed   13 1.8616   1.7611     39     42 0.03699 *
# Residuals 14

allSml <- allSml %>% arrange(subSeed, optimRepl, nOptIter)
allSml <- allSml %>% dplyr::mutate(diffPIC=budgPIC - lag(budgPIC))
allSml <- allSml %>%
  dplyr::mutate(diffPIC=if_else(nOptIter==36, as.numeric(NA), diffPIC))
plot(allSml$nOptIter, allSml$diffPIC)
asr <- allSml %>% group_by(nOptIter) %>% filter(nOptIter != 36) %>% summarize(diffPICrange=range(diffPIC))
badConv <- allSml %>%
  filter(nOptIter > 1300) %>%
  filter(abs(diffPIC) > 0.02) %>%
  pull(subSeed) %>% unique

goodConv <- allSml %>% filter(!(subSeed %in% badConv))
man <- manova(cbind(budgPIC, budgSDN, budgCET) ~ subSeed,
              data=goodConv, subset=(nOptIter==1536))
summary(man)
fitSubSeed <- lm(budgPIC ~ subSeed, data=goodConv, subset=(nOptIter==1536))
anova(fitSubSeed)
fitSubSeed <- lm(budgSDN ~ subSeed, data=goodConv, subset=(nOptIter==1536))
anova(fitSubSeed)
fitSubSeed <- lm(budgCET ~ subSeed, data=goodConv, subset=(nOptIter==1536))
anova(fitSubSeed)

library(lme4)
varComp <- NULL
for (nIter in allIter){
  tst <- lmer(budgPIC ~ (1 | subSeed), data=allSml, subset=(nOptIter==nIter))
  vc <- VarCorr(tst)
  varComp <- rbind(varComp, c(attr(vc, "sc"), sqrt(vc$subSeed)))
}
maxVC <- max(varComp)
plot(varComp[,1], ylim=c(0, maxVC), pch=16)
points(varComp[,2], pch=16, col="dark green")

############# All Big below here ##################

setwd("/Users/jj332/Documents/simOut/bso/output/PostGrid")
allBig <- list.files(pattern="big")
# Parse allBig of this form 188690362_0_1036sml06PIC.rds
where_ <- gregexpr("_", allBig)
whereb <- regexpr("b", allBig)
parseAllBig <- tibble()
for (i in 1:length(allBig)){
  w_ <- where_[[i]]
  ss <- substring(allBig[i], 1, w_[1] - 1)
  if (length(w_) > 1){
    opr <- substring(allBig[i], w_[1] + 1, w_[1] + 1)
    w_ <- w_[2]
  } else opr <- "0"
  ni <- as.numeric(substring(allBig[i], w_ + 1, whereb[i] - 1))
  nc <- substring(allBig[i], whereb[i] + 3, whereb[i] + 4)
  parseAllBig <- dplyr::bind_rows(parseAllBig,
                                  tibble(fileName=allBig[i],
                                         subSeed=ss,
                                         optimRepl=opr,
                                         nOptIter=ni,
                                         nCyc=nc)
                                  )
}

allRes <- NULL
for (f in allBig){
  res <- readRDS(f)
  allRes <- rbind(allRes, c(res$bestBudget, res$maxPredGain, res$postVarAtMax))
}
colnames(allRes) <- c("budgPIC", "budgSDN", "budgCET",
                      "maxPredGain", "postVarAtMax")
allBig <- dplyr::bind_cols(parseAllBig, allRes)
allIter <- allBig %>% pull(nOptIter) %>% unique %>% sort

bigMeans <- allBig %>% group_by(nOptIter, nCyc) %>%
  summarize(
    meanPIC=mean(budgPIC),
    meanSDN=mean(budgSDN),
    meanCET=mean(budgCET),
    meanPYT=mean(1-budgPIC-budgSDN-budgCET),
    meanGain=mean(maxPredGain),
    meanVar=mean(postVarAtMax),
    sdPIC=sd(budgPIC),
    sdSDN=sd(budgSDN),
    sdCET=sd(budgCET),
    sdPYT=sd(budgPIC+budgSDN+budgCET),
    sdGain=sd(maxPredGain),
    sdVar=sd(postVarAtMax),
  )

yRange <- bigMeans %>% ungroup %>%
  select(contains("mean")) %>% select(-meanGain, -meanVar) %>%
  range

bm06 <- bigMeans %>% filter(nCyc=="06")
pdf(here::here("output", "Allocation06Cyc.pdf"))
plot(bm06$nOptIter, bm06$meanPIC, pch=16, ylim=yRange,
     xlab="Optimization iteration", ylab="Budget allocation",
     main="Budget allocation for a 6 cycle horizon",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(bm06$nOptIter, bm06$meanSDN, pch=16, col="dark red", cex=1.3)
points(bm06$nOptIter, bm06$meanCET, pch=16, col="dark green", cex=1.3)
points(bm06$nOptIter, bm06$meanPYT, pch=16, col="blue", cex=1.3)
legend(280, 0.34, c("PIC", "SDN", "CET", "PYT"), pch=16,
       col=c("black", "dark red", "dark green", "blue"),
       cex=1.3)
dev.off()

bm12 <- bigMeans %>% filter(nCyc=="12")
pdf(here::here("output", "Allocation12Cyc.pdf"))
plot(bm12$nOptIter, bm12$meanPIC, pch=16, ylim=yRange,
     xlab="Optimization iteration", ylab="Budget allocation",
     main="Budget allocation for a 12 cycle horizon",
     cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(bm12$nOptIter, bm12$meanSDN, pch=16, col="dark red", cex=1.3)
points(bm12$nOptIter, bm12$meanCET, pch=16, col="dark green", cex=1.3)
points(bm12$nOptIter, bm12$meanPYT, pch=16, col="blue", cex=1.3)
legend(280, 0.21, c("PIC", "SDN", "CET", "PYT"), pch=16,
       col=c("black", "dark red", "dark green", "blue"),
       cex=1.3)
dev.off()

# Plot the difference in the means between 6-cycle versus 12-cycle simulation
# The difference seems to have more or less stabilized by iteration 200
picDiff <- bigMeans %>% filter(nCyc=="06") %>% pull(meanPIC) -
  bigMeans %>% filter(nCyc=="12") %>% pull(meanPIC)
sdnDiff <- bigMeans %>% filter(nCyc=="06") %>% pull(meanSDN) -
  bigMeans %>% filter(nCyc=="12") %>% pull(meanSDN)
cetDiff <- bigMeans %>% filter(nCyc=="06") %>% pull(meanCET) -
  bigMeans %>% filter(nCyc=="12") %>% pull(meanCET)
pytDiff <- bigMeans %>% filter(nCyc=="06") %>% pull(meanPYT) -
  bigMeans %>% filter(nCyc=="12") %>% pull(meanPYT)
budgDiffs <- tibble(nOptIter=36+seq(0, 300, 50), picDiff=picDiff,
                    sdnDiff=sdnDiff, cetDiff=cetDiff, pytDiff=pytDiff)
rangeDiff <- budgDiffs %>% dplyr::select(-nOptIter) %>% unlist %>% range
plot(budgDiffs$nOptIter, budgDiffs$picDiff, pch=16, ylim=rangeDiff)
points(budgDiffs$nOptIter, budgDiffs$sdnDiff, pch=16, col="dark red")
points(budgDiffs$nOptIter, budgDiffs$cetDiff, pch=16, col="dark green")
points(budgDiffs$nOptIter, budgDiffs$pytDiff, pch=16, col="blue")

# Plot the std dev of the difference between 6-cycle versus 12-cycle simulation
# The take to 250 iterations to stabilize.  They do not go to zero which again
# suggests a founder effect: optimizations on different founder populations do
# not converge to the same budgets
picSD <- ((bigMeans %>% filter(nCyc=="06") %>% pull(sdPIC))^2 +
  (bigMeans %>% filter(nCyc=="12") %>% pull(sdPIC))^2) %>% sqrt
sdnSD <- ((bigMeans %>% filter(nCyc=="06") %>% pull(sdSDN))^2 +
  (bigMeans %>% filter(nCyc=="12") %>% pull(sdSDN))^2) %>% sqrt
cetSD <- ((bigMeans %>% filter(nCyc=="06") %>% pull(sdCET))^2 +
            (bigMeans %>% filter(nCyc=="12") %>% pull(sdCET))^2) %>% sqrt
pytSD <- ((bigMeans %>% filter(nCyc=="06") %>% pull(sdPYT))^2 +
            (bigMeans %>% filter(nCyc=="12") %>% pull(sdPYT))^2) %>% sqrt
budgDiffs <- budgDiffs %>% mutate(picSD=picSD,
                    sdnSD=sdnSD, cetSD=cetSD, pytSD=pytSD)
rangeDiff <- budgDiffs %>% dplyr::select(contains("SD", ignore.case=F)) %>%
  unlist %>% range
rangeDiff[1] <- 0
plot(budgDiffs$nOptIter, budgDiffs$picSD, pch=16, ylim=rangeDiff)
points(budgDiffs$nOptIter, budgDiffs$sdnSD, pch=16, col="dark red")
points(budgDiffs$nOptIter, budgDiffs$cetSD, pch=16, col="dark green")
points(budgDiffs$nOptIter, budgDiffs$pytSD, pch=16, col="blue")

# Collect the p-values over iterations to see when it becomes significant
# I should just collect the effect to see when it stabilizes...
# And that I can just get from bigMeans
pvalOverIter <- NULL
for (nIter in allIter){
  pic <- lm(budgPIC ~ nCyc, data=allBig, subset=(nOptIter==nIter))
  sdn <- lm(budgSDN ~ nCyc, data=allBig, subset=(nOptIter==nIter))
  cet <- lm(budgCET ~ nCyc, data=allBig, subset=(nOptIter==nIter))
  man <- manova(cbind(budgPIC, budgSDN, budgCET) ~ nCyc,
                data=allBig, subset=(nOptIter==nIter))
  pvalOverIter <- rbind(pvalOverIter,
                        c(summary(pic)$coefficients[2, 1],
                          summary(sdn)$coefficients[2, 1],
                          summary(cet)$coefficients[2, 1],
                          summary(man)$stats[1, 6]))
}

plot()

#########################################################
# Everything below here is from Oct 2021
#########################################################
# Parse to show which are repeat optimizations versus founder sets
allPostSm <- allPostSm %>% as_tibble(allPostSm)
allPostSm <- allPostSm %>% filter(!(initNum == 7 | initNum == 14)) %>%
  dplyr::mutate(initNum=if_else(initNum > 6, initNum-1, initNum)) %>%
  dplyr::mutate(replication=as.factor((initNum-1) %% 2),
                realInit=as.factor((initNum-1) %/% 2)) %>%
  dplyr::mutate(PYT=1 - (PIC+SDN+CET)) %>% relocate(PYT, .after=CET)

# Calculate the distance between the the maximum for the optimization at nIter
# iterations versus the final maximum after 1500 iterations, for both same and
# independent optimizations
allPostSm <- dplyr::mutate(allPostSm,
                           distSameOpt=numeric(nrow(allPostSm)),
                           distOthrOpt=numeric(nrow(allPostSm)))
for (i in unique(allPostSm$initNum)){
  sameMax <- finalPostMax %>%
    dplyr::filter(realInit==(i-1) %/% 2 & replication==(i-1) %% 2) %>%
    dplyr::select(PIC, SDN, CET, PYT)
  othrMax <- finalPostMax %>%
    dplyr::filter(realInit==(i-1) %/% 2 & replication==i %% 2) %>%
    dplyr::select(PIC, SDN, CET, PYT)

  distFori <- rowwise(allPostSm) %>%
    dplyr::transmute(d=(c(PIC, SDN, CET, PYT) - sameMax) %>% abs %>% sum) %>%
    pull(d)
  distForNoti <- pull(allPostSm, distSameOpt)
  allPostSm <- dplyr::mutate(allPostSm, distSameOpt=
                               if_else(initNum==i,distFori, distForNoti))

  distFori <- rowwise(allPostSm) %>%
    dplyr::transmute(d=(c(PIC, SDN, CET, PYT) - othrMax) %>% abs %>% sum) %>%
    pull(d)
  distForNoti <- pull(allPostSm, distOthrOpt)
  allPostSm <- dplyr::mutate(allPostSm, distOthrOpt=
                               if_else(initNum==i,distFori, distForNoti))
}

finalPostMax <- allPostSm %>% filter(nIter==max(nIter))

meanPost <- allPostSm %>% dplyr::group_by(nIter) %>% summarize(same=mean(distSameOpt), othr=mean(distOthrOpt))
plot(meanPost[,1:2], pch=16, ylim=range(meanPost[,2:3]))
points(meanPost[,c(1,3)], pch=16, col="red", ylim=range(meanPost[,2:3]))

checkInit <- function(afterNIter){
  postAfterIter <- allPostSm %>% filter(nIter == afterNIter)
  tstInitCond <- lm(PIC ~ as.factor(realInit), data=postAfterIter)
  return(as.matrix(anova(tstInitCond))[1, 5])
}

whenInitMatters <- sapply(1:30*50, checkInit)
plot(whenInitMatters)

# This is just a dumb ANOVA to see if there is difference in PIC across
# different founder populations
tst <- lm(PIC ~ realInit, data=allPostSm, subset=(nIter==1500))
summary(tst)
aov(tst)
# Apparently: p-value: 0.01354
#                   realInit  Residuals
# Sum of Squares  0.07392670 0.01146312
# Deg. of Freedom          5          6
library(lme4)
varComp <- NULL
for (nOptIter in 50*1:30){
  tst <- lmer(PIC ~ (1 | realInit), data=allPostSm, subset=(nIter==nOptIter))
  vc <- VarCorr(tst)
  varComp <- rbind(varComp, c(attr(vc, "sc"), sqrt(vc$realInit)))
}
maxVC <- max(varComp)
plot(varComp[,1], ylim=c(0, maxVC), pch=16)
points(varComp[,2], pch=16, col="dark green")

# This is a MANOVA for the same purpose: also no effect
tst <- manova(cbind(PIC, SDN, CET) ~ realInit,
              data=allPostSm, subset=(nIter==1500))
summary(tst)
# Also for MANOVA
#           Df Pillai approx F num Df den Df   Pr(>F)
# realInit   5 2.2253    3.447     15     18 0.007118 **

