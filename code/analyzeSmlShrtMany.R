library(tidyverse)
setwd("/Users/jj332/Documents/GitRepo/BayesOpt_ScriptsForMS")
here::i_am("code/analyzeSmlShrtMany.R")

# Get information on outputs at Ceres
has499Iter <- read_csv(here::here("output", "has499Iter.csv"))
allSubSeed <- has499Iter %>% pull(subSeed) %>% unique

# Get information on the outputs that I actually have
allSml <- list.files(path=here::here("output", "PostGrid"), pattern="sml")
# Parse allSml of this form 188690362_0_1036sml06PIC.rds
where_ <- regexpr("_", allSml)
wheres <- regexpr("s", allSml)
subSeed <- substring(allSml, 1, where_ - 1)
optimRepl <- substring(allSml, where_ + 1, where_ + 1)
nIter <- as.numeric(substring(allSml, where_ + 3, wheres - 1))

allRes <- NULL
for (f in allSml){
  res <- readRDS(here::here("output", "PostGrid", f))
  allRes <- rbind(allRes, c(res$bestBudget, res$maxPredGain, res$postVarAtMax))
}
colnames(allRes) <- c("budgPIC", "budgSDN", "budgCET",
                      "maxPredGain", "postVarAtMax")
allSml <- tibble(fileName=allSml, subSeed=subSeed,
                 optimRepl=optimRepl, nOptIter=nIter)
allSml <- dplyr::bind_cols(allSml, allRes)
allSml <- allSml %>% mutate(budgPYT = 1 - (budgPIC + budgSDN + budgCET)) %>%
  relocate(budgPYT, .after=budgCET) %>%
  arrange(subSeed, optimRepl, nOptIter)

print("All the random seeds")
hereSubSeed <- allSml %>% pull(subSeed) %>% unique
print(hereSubSeed)
print("Missing subSeeds")
has499Iter %>% filter(!(subSeed %in% hereSubSeed)) %>%
  pull(subSeed) %>% unique %>% print

fitSubSeed <- lm(budgPIC ~ subSeed, data=allSml, subset=(nOptIter==1536))
anova(fitSubSeed)
fitSubSeed <- lm(budgSDN ~ subSeed, data=allSml, subset=(nOptIter==1536))
anova(fitSubSeed)
fitSubSeed <- lm(budgCET ~ subSeed, data=allSml, subset=(nOptIter==1536))
anova(fitSubSeed)
fitSubSeed <- lm(budgPYT ~ subSeed, data=allSml, subset=(nOptIter==1536))
anova(fitSubSeed)
# Response: budgPIC
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# subSeed   15 0.037751 0.0025167  0.8903 0.5868
# Residuals 16 0.045229 0.0028268
# Response: budgSDN
#           Df   Sum Sq   Mean Sq F value    Pr(>F)
# subSeed   15 0.181928 0.0121285    7.62 0.0001086 ***
# Residuals 16 0.025467 0.0015917
# Response: budgCET
#           Df   Sum Sq   Mean Sq F value  Pr(>F)
# subSeed   15 0.143559 0.0095706  2.8105 0.02413 *
# Residuals 16 0.054484 0.0034053
# Response: budgPYT
#           Df    Sum Sq    Mean Sq F value Pr(>F)
# subSeed   15 0.0137745 0.00091830  1.5915 0.1832
# Residuals 16 0.0092322 0.00057701

man <- manova(cbind(budgPIC, budgSDN, budgCET) ~ subSeed,
              data=allSml, subset=(nOptIter==1536))
summary(man)
#           Df Pillai approx F num Df den Df  Pr(>F)
# subSeed   15 1.8536   1.7246     45     48 0.03255 *
# Residuals 16
summary(man, test="Wilks")
#           Df    Wilks approx F num Df den Df   Pr(>F)
# subSeed   15 0.029978   2.1246     45 42.371 0.007449 **
# Residuals 16

# There is evidence of difference by Founder pop.  It is pretty weak
# If I remove Founder pops where there is some evidence of bad optimization
# convergence, does that change the picture? Essentially, no. So leave them in
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

# After how many iterations of optimization does the Founder effect appear?
library(lme4)
approxF <- NULL
allIter <- allSml %>% pull(nOptIter) %>% unique %>% sort
for (nIter in allIter){
  man <- manova(cbind(budgPIC, budgSDN, budgCET) ~ subSeed,
                data=allSml, subset=(nOptIter==nIter))
  sman <- summary(man)
  # sman <- summary(man, test="Wilks")
  approxF <- c(approxF, sman$stats["subSeed", "approx F"])
}
pdf(here::here("output", "figures", "FounderEffectApproxFoverIterations.pdf"))
plot(allIter, approxF, pch=16, xlab="Optimization iteration", ylab="Pillai test approx. F stat.", cex.lab=1.3, cex.axis=1.3, cex=1.3)
dev.off()

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

