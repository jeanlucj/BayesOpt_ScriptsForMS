# Source this to retrieve stored optimization iterations
require(here)
here::i_am("code/getOptimizationIterations.R")

# init_num can include the optimization replication
# So something like 27823605_0
if (!exists("init_num")) init_num <- 1
print(paste("Get optimization iterations from Initialization", init_num))

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
