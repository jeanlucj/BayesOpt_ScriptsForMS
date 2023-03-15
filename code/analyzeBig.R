
library(tidyverse)
setwd("/Users/jj332/Documents/GitRepo/BayesOpt_ScriptsForMS")
here::i_am("code/analyzeBig.R")

# Comes out of SortThroughOutputs but doctored in XL to show which were PIC1 etc
allOpt <- read_csv(here::here("output", "has99IterXL.csv"))
allOpt$subSeed <- as.character(allOpt$subSeed)
# Have to get rid of the b06 and b12
setwd(here::here("output", "PostGrid", "Feb17"))
allBig <- list.files(pattern="big")
allBig <- allBig[-grep("b06", allBig)]
allBig <- allBig[-grep("b12", allBig)]
# Parse allBig of this form 970854781_336big06PIC.rds
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
allBig <- dplyr::full_join(allBig,
         allOpt %>% filter(Budget_Full_Shrt=="Full") %>%
           select(subSeed, nPICcycles, Budget_Full_Shrt))

# Cool. Shows that the relative budgets depend on the number of PIC cycles
man <- manova(cbind(budgPIC, budgSDN, budgCET) ~ nCyc*nPICcycles,
              data=allBig, subset=(nOptIter==336))
summary(man)
#                 Df  Pillai approx F num Df den Df  Pr(>F)
# nCyc             1 0.10308   1.1492      3     30 0.34532
# nPICcycles       1 0.21163   2.6845      3     30 0.06438 .
# nCyc:nPICcycles  1 0.28711   4.0273      3     30 0.01607 *
# Residuals       32

allIter <- allBig %>% pull(nOptIter) %>% unique %>% sort
allBig %>% pull(nOptIter) %>% table
bigMeans <- allBig %>% group_by(nOptIter, nCyc, nPICcycles) %>%
  summarize(
    nOpt=n(),
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

# Save information for Results table
big336 <- bigMeans %>% select(-contains("Var")) %>%
  mutate(nCyc=as.numeric(nCyc)) %>%
  filter(nOptIter == 336) %>% round(digits=3)
write_csv(big336, file=here::here("output", "_PubOut", "BigOptBudg.csv"))

##################################################
stop("Other stuff below")
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
