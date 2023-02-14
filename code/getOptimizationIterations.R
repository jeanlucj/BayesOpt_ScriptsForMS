# Source this to retrieve stored optimization iterations
require(here)
require(tidyverse)
here::i_am("code/getOptimizationIterations.R")

if (!exists("init_num")) init_num <- 1
print(paste("Get optimization iterations from Initialization", init_num))

# init_num can include the optimization replication
# So something like 27823605_0
if (init_num %in% c("b06", "b12", "sml")){
  # Here, bring together optimizations from all founder populations
  simDF <- ifelse(init_num == "sml", "has499Iter.csv", "has99Iter.csv")
  simDF <- readr::read_csv(here::here("output", simDF))
  if (init_num != "sml"){
    if (init_num == "b06"){
      simDF <- simDF %>% dplyr::filter(nCycSet == 6)
    } else{
      simDF <- simDF %>% dplyr::filter(nCycSet == 12)
    }
  }
  btf <- simDF %>% dplyr::pull(fileBoTorch)
  budgets <- NULL
  gains <- NULL
  for (fileName in btf){
    initDat <- readRDS(here::here("output", fileName))
    where_ <- regexpr("_", fileName)
    if (where_ < 0){
      optimRepl <- 0
    } else{
      optimRepl <- as.numeric(substring(fileName, where_+1, where_+1))
    }
    budg <- initDat[[1]][[optimRepl+1]][[1]]
    if (ncol(budg) == 4) budg <- budg[, -4]
    budg <- budg[1:n_iter,]
    budgets <- rbind(budgets, as.matrix(budg))

    gns <- initDat[[2]][[optimRepl+1]][[1]]
    gns <- gns[1:n_iter]
    gains <- c(gains, gns)
  }
  # Adjust gains for their Founder effect
  founder <- rep(btf, each=n_iter)
  fitFound <- lm(gains ~ founder)
  gains <- gains - predict(fitFound)
} else{
  # Here, just get posteriors on ONE founder population (original version...)
  initDat <- readRDS(here::here("output", paste0("/BoTorchOut", init_num, ".rds")))
  where_ <- regexpr("_", init_num)
  if (where_ < 0){
    optimRepl <- 0
  } else{
    optimRepl <- as.numeric(substring(init_num, where_+1, where_+1))
  }
  budgets <- initDat[[1]][[optimRepl+1]][[1]]
  gains <- initDat[[2]][[optimRepl+1]][[1]]
  budgets <- as.matrix(budgets)
  if (ncol(budgets) == 4) budgets <- budgets[, -4]
  budgets <- budgets[1:n_iter,]
  gains <- gains[1:n_iter]
}

