# Process the outputs from getPosteriors.py
require(here)
here::i_am("code/processPosteriors.R")

if (!exists("init_num")) init_num <- 1
print(paste("Processing posterior predictions from Initialization", init_num))

posteriors <- list(bestBudget=bestBudget,
                   maxPredGain=maxPredGain,
                   postVarAtMax=postVarAtMax)

if (exists("predGains")){
  posteriors <- c(posteriors,
                  list(predBudgets=predBudgets,
                       predGains=predGains,
                       predVars=predVars)
  )
}

saveRDS(posteriors,
        file=here::here("output",
                        paste0("posteriors", init_num, "_", n_iter, ".rds")))
