# When I pull down simulation outputs from ceres, I need to sort through what
# has been done and put into different sets for analyses.
# setwd("/Users/jj332/Documents/GitRepo/BayesOpt_ScriptsForMS")
here::i_am("code/SortThroughOutputs.R")
# setwd("/Users/jj332/Documents/simOut/bso/output")
library(tidyverse)

ceresJobs <- readLines(here::here("ceresSLURMjobs.txt"))
whereSkip <- which(ceresJobs == "") %>% max # WARNING somewhat brittle

ceresJobs <- read_csv(here::here("ceresSLURMjobs.txt"),
                      col_names=F,
                      col_types="iciciccc",
                      skip=whereSkip)
colnames(ceresJobs) <- c("initSeed", "subSeed", "nCycSet", "timePrep",
                         "slurmJobID", "status", "time", "notes")
# Retain only the last instance of an initSeed
ceresJobs <- ceresJobs %>% filter(!duplicated(initSeed, fromLast=T))

# The NBO files ("new batch out") just looks at the first 36 budgets
# It records the complete set of breeding cycles, and so counts those cycles
# I'm using it to check on the number of cycles were run to ensure that I got
# that right. NOTE: only named relative to the subSeed. None will be duplicated
allNBO <- list.files(path=here::here("output"), pattern="NBO")
allNcyc <- NULL
for (nbo in allNBO){
  tst <- readRDS(here::here("output", nbo))
  allNcyc <- c(allNcyc, nrow(tst[[1]][[2]]) - 1)
}

allNBOsubSeed <- gsub("bgNBO", "", allNBO)
allNBOsubSeed <- gsub(".rds", "", allNBOsubSeed)
allNBO <- tibble(fileNBO=allNBO, subSeed=allNBOsubSeed, nCycNBO=allNcyc)

# outputInfo will collect all the information being sorted here
outputInfo <- full_join(allNBO, ceresJobs, by="subSeed")

# The BoTorch files are the final output of the optimization
# Checking here for the number of optimization iterations that were done
# NOTE: BoTorch files can be duplicated because a subSeed might be optimized
# independently twice.  If a file has NA for the optimization replication while
# another file with the same subSeed is not NA for the optimRepl, the NA file
# should be discarded
allNiter <- NULL
allBoTorch <- list.files(path=here::here("output"), pattern="BoTorch")
for (bt in allBoTorch){
  nIter <- try(readRDS(here::here("output", bt))[[4]], silent=T)
  if (class(nIter) != "try-error") allNiter <- c(allNiter, nIter)
}

allBTsubSeed <- gsub("BoTorchOut", "", allBoTorch)
allBTsubSeed <- gsub(".rds", "", allBTsubSeed)
where_ <- regexpr("_", allBTsubSeed)
allBTsubSeed1 <- mapply(
  FUN=function(w_, s) if_else(w_ < 0, s, substring(s, 1, w_-1)),
  where_, allBTsubSeed)
allBToptimRepl <- mapply(
  FUN=function(w_, s) if_else(w_ < 0, "", substring(s, w_+1, w_+1)),
  where_, allBTsubSeed)

allBT <- tibble(fileBoTorch=allBoTorch,
                subSeed=allBTsubSeed1,
                optimRepl=allBToptimRepl,
                nOptIter=allNiter)
# Get rid of files that don't have the optimRepl if they should have it
allBT <- allBT %>% filter(
  !((duplicated(subSeed) | duplicated(subSeed, fromLast=T)) &
    optimRepl == "")
  )

outputInfo <- full_join(outputInfo, allBT, by="subSeed")

# I have saved the simulation parameters in directories that are named according
# to the subSeed.  Useful to know if those exist
areDirectories <- list.files(path=here::here("output"))
areDirectories <- areDirectories[-grep(".", areDirectories, fixed=T)]
hasDirectory <- tibble(subSeed=areDirectories, hasDirectory=T)

outputInfo <- full_join(outputInfo, hasDirectory, by="subSeed")

# Now, the optimization prints out all the parameters in the optOutFile
# Check if that's available and pull the most critical (changing) parameters
# The optOutFile will only be named by the subSeed so none should be duplicated
# even if the optimization was replicated
allOOF <- list.files(path=here::here(), pattern="optOutFile")
allOOFparms <- NULL
for (oof in allOOF){
  tmp <- readLines(here::here(oof))
  getParms <- function(parmName){ # Not used at the moment
    pline <- tmp[grep(parmName, tmp)]
    parms <- strsplit(pline, " ") %>% unlist
    return(parms[-1])
  }
  if (length(grep("nBreedingProg", tmp)) > 0){
    parmsOfInterest <- NULL
    for (parmName in c("nCyclesToRun", "nEntries", "nBreedingProg")){
      parmsOfInterest <- paste0(parmsOfInterest, tmp[grep(parmName, tmp)], " ")
    }
  } else parmsOfInterest <- "NA"
  allOOFparms <- c(allOOFparms, parmsOfInterest)
}
allOOFsubSeed <- gsub("optOutFile", "", allOOF)
allOOFsubSeed <- gsub(".txt", "", allOOFsubSeed)
allOOF <- tibble(fileOOF=allOOF, subSeed=allOOFsubSeed, parms=allOOFparms)

outputInfo <- full_join(outputInfo, allOOF, by="subSeed")

has99Iter <- outputInfo %>%
  filter(nOptIter == 99) %>%
  filter(!duplicated(subSeed))
write_csv(has99Iter, here::here("output", "has99Iter.csv"), quote="none")

has499Iter <- outputInfo %>% filter(nOptIter == 499) %>%
  slice(grep("100 40 16", parms))
# Make sure all the subSeeds you keep have the two optimizations
checkHas01optim <- function(chkSub, tb499){
  or <- tb499 %>% dplyr::filter(subSeed==chkSub)
  return(all(0:1 %in% or$optimRepl))
}
keepSub <- has499Iter$subSeed %>% unique
keepSub <- keepSub[sapply(keepSub, FUN=checkHas01optim, has499Iter)]
has499Iter <- has499Iter %>% dplyr::filter(subSeed %in% keepSub)

write_csv(has499Iter, here::here("output", "has499Iter.csv"), quote="none")
