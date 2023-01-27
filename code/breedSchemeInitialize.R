## Need to do this to get the budget constraints
require(here)
here::i_am("code/breedSchemeInitialize.R")

if (!exists("init_num")) init_num <- 1
print(paste("Initializing breeding scheme from Initialization", init_num))

## Check if there is a directory in output named init_num where the simulation
## parameters were stored
# Get subSeed from init_num
where_ <- regexpr("_", init_num)
if (where_ > 0) init_num <- substring(init_num, 1, where_ - 1)
# Check directories in output
fileOutput <- list.files(path=here::here("output"))
fileOutput <- fileOutput[-grep(".", fileOutput, fixed=T)]
if (init_num %in% fileOutput){
  subFolder <- paste0("output/", init_num)
} else{
  subFolder <- "data"
}

## ----Initialize program-----------------------------------------------
bsd <- BreedSimCost::initializeProgram(
  here::here(subFolder, "FounderCtrlFile.txt"),
  here::here(subFolder, "SchemeCtrlFile.txt"),
  here::here(subFolder, "CostsCtrlFile.txt"),
  here::here(subFolder, "OptimizationCtrlFile.txt")
  )

budget_constraints <- bsd$initBudget[c("minPICbudget", "minLastStgBudget")]
budget_constraints <- c(budget_constraints, bsd$initBudget[grep("ratio", names(bsd$initBudget))])
