# Visualize postGrid
require(here)
library(tidyverse)
here::i_am("code/visualizePostGrid.R")
has499 <- read_csv(here::here("data", "has499Iter.csv"))
allSubSd <- has499 %>% pull(subSeed) %>% unique

# Find out if I have postGrid of all those subseeds
for (ss in allSubSd){
  ssFiles <- list.files(here::here("output", "PostGrid"), pattern=as.character(ss))
  if (length(ssFiles)==0) cat("No files for", ss, "\n")
}
# YES: I have them all

# For each subSeed, find postGrid points that are close to the max - SD/10
# The standard deviation is SD=5.107004
stdDevAtMax.1 <- 0.5107
budgSml <- c(0.36805244, 0.33554188, 0.19840319, 0.09800249)
files1536 <- list.files(here::here("output", "PostGrid"), pattern="1536")
for (fName in files1536){
  pg <- readRDS(here::here("output", "PostGrid", fName))
  distToMax <- abs(c(pg$predGains) - c(pg$maxPredGain))
  allocBelowMax <- order(distToMax) %>% .[1:20] %>% sort
  predAlcBelMx <- pg$predBudgets[allocBelowMax,]
  predAlcBelMx <- cbind(predAlcBelMx, apply(predAlcBelMx, 1, function(v) 1 - sum(v)))
}
