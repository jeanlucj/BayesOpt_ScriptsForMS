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

# I think I really need to run both the average and the all together to see
# which comes out ahead.  They aren't so different for the 6 cycle optimum
# but they are pretty different for the 12 cycle optimum
# For 6 cycles analyzed all together: 0.4169595 0.4123836 0.15307223 0.01758467
# nEntries: 920 152 4, nBreedingProg: 844
# For 12 cycles analyzed all together: 0.2523531 0.3416312 0.38409934 0.02191636
# nEntries: 762 382 5, nBreedingProg: 511
# For 6 cycles: 0.3858742 0.4231775 0.1392178 0.05173058 [after 8 reps]
# For 12 cycles: 0.2620301 0.3966589 0.2810439 0.06026706 [after 8 reps]
# For 6 cycles: 0.3874482 0.3870322 0.1624069 0.06311256 [after 6 reps]
# For 12 cycles: 0.2520365 0.3929344 0.2805345 0.07449453 [after 6 reps]
# BUT you want to use each of these budgets for 12 cycles
budgAT6Cyc <- c(0.41512103, 0.41056531, 0.15239730, 0.02191636)
budgAT12Cyc <- c(0.2523531, 0.3416312, 0.38409934, 0.02191636) # Used
budg6Cyc <- c(0.3874482, 0.3870322, 0.1624069, 0.06311256) # 6 reps Used
budg12Cyc <- c(0.2520365, 0.3929344, 0.2805345, 0.07449453) # 6 reps
budg6Cyc <- c(0.3858742, 0.4231775, 0.1392178, 0.05173058) # 8 reps
budg12Cyc <- c(0.2620301, 0.3966589, 0.2810439, 0.06026706) # 8 reps
budgP1Cyc6 <- c(0.3235522, 0.4163751, 0.2357820, 0.02429063) # Used
budgP1Cyc12 <- c(0.3836027, 0.4179224, 0.1791765, 0.01929836) # Used
budg6Cyc <- c(0.3874482, 0.3870322, 0.1624069, 0.06311256)
budg12Cyc <- c(0.2520365, 0.3929344, 0.2805345, 0.07449453)
budgSml <- c(0.36805244, 0.33554188, 0.19840319, 0.09800249)
perc <- budg6Cyc
tstAgain <- runBatchBSD(burnedInPops, percentage=perc)
saveRDS(here::here("output", paste0("tstAgain", perc[1], ".rds")))
