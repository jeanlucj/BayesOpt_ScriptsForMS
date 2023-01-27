require(here)
here::i_am("code/RunShortLongSimulationsAtOptima.R")
library(BreedSimCost)

runBurnIn <- function(dummy){
  ## ----Initialize program-----------------------------------------------
  bsd <- initializeProgram("data/FounderCtrlFileSml.txt",
                           "data/SchemeCtrlFileSml.txt",
                           "data/CostsCtrlFile.txt",
                           "data/OptimizationCtrlFileSml.txt")

  ## ----Fill variety development pipeline------------------
  # Year 1
  bsd$year <- bsd$year+1
  bsd <- makeVarietyCandidates(bsd)

  bsd$entries <- bsd$varietyCandidates@id
  bsd <- runVDPtrial(bsd, "SDN")

  parents <- selectParentsBurnIn(bsd)
  bsd <- makeCrossesBurnIn(bsd, parents)

  # Year 2
  bsd$year <- bsd$year+1
  bsd <- makeVarietyCandidates(bsd)

  bsd <- chooseTrialEntries(bsd, toTrial="SDN")
  bsd <- runVDPtrial(bsd, "SDN")
  bsd <- chooseTrialEntries(bsd, fromTrial="SDN", toTrial="CET")
  bsd <- runVDPtrial(bsd, "CET")

  parents <- selectParentsBurnIn(bsd)
  bsd <- makeCrossesBurnIn(bsd, parents)

  # Year 3 and onward
  for (burnIn in 1:bsd$nBurnInCycles){
    bsd$year <- bsd$year+1
    bsd <- makeVarietyCandidates(bsd)

    bsd <- chooseTrialEntries(bsd, toTrial="SDN")
    bsd <- runVDPtrial(bsd, "SDN")
    bsd <- chooseTrialEntries(bsd, fromTrial="SDN", toTrial="CET")
    bsd <- runVDPtrial(bsd, "CET")
    bsd <- chooseTrialEntries(bsd, fromTrial="CET", toTrial="PYT")
    bsd <- runVDPtrial(bsd, "PYT")

    parents <- selectParentsBurnIn(bsd)
    bsd <- makeCrossesBurnIn(bsd, parents)
  }
  return(bsd)
}

burnedInPops <- lapply(1:4, runBurnIn)

bsd <- burnedInPops[[1]]
budget_constraints <- bsd$initBudget[c("minPICbudget", "minLastStgBudget")]
budget_constraints <- c(budget_constraints, bsd$initBudget[grep("ratio", names(bsd$initBudget))])

runBatchBSD <- function(batchBSD, percentage){
  require(parallel)
  if (bsd$debug){
    require(here)
    on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data", "runBatch.rds"))
    })
  }

  runWithBSD <- function(bsd, percentage, returnBSD=F){
    return(runWithBudget(percentage, bsd, returnBSD))
  }

  if (bsd$debug){
    batchResults <- lapply(batchBSD, runWithBSD, percentage=percentage)
  } else{
    batchResults <- mclapply(batchBSD, runWithBSD, percentage=percentage,
                             mc.preschedule = F, mc.cores=bsd$nCores)
  }
  # batchResults is now a list of lists
  # Remove results where budget was not valid
  findBatchNA <- function(oneRunRes){
    !any(is.na(unlist(oneRunRes)))
  }
  batchResults <- batchResults[sapply(batchResults,
                                      function(v) findBatchNA(v))]

  if (bsd$debug) on.exit()
  return(batchResults)
}#END runBatchBSD

# For 6 cycles: 0.3874482 0.3870322 0.1624069 0.06311256
# For 12 cycles: 0.2520365 0.3929344 0.2805345 0.07449453
# BUT you want to use each of these budgets for 12 cycles
budg6Cyc <- c(0.3874482, 0.3870322, 0.1624069, 0.06311256)
budg12Cyc <- c(0.2520365, 0.3929344, 0.2805345, 0.07449453)
perc <- budg6Cyc
tstAgain <- runBatchBSD(burnedInPops, percentage=perc)
saveRDS(here::here("output", paste0("tstAgain", perc[1], ".rds")))
