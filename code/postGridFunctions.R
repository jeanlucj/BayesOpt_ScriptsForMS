require(here)
here::i_am("code/postGridFunctions.R")

# Script to make a ternary contour plot showing gain with different budgets
n_stages <- 3
# Make the budget grid to get posteriors from the Bayesian Optimization
# 1. Use the same constraints as in the BO.
# a) Make a rudimentary grid with a fairly fine step
# b) Eliminate grid points that do not satisfy using the outConstraint function
# outConstraint applies the same constraints as in OneOptimization.py
# c) Reduce the number of grid points by doing kmeans clustering
# d) Send the grid points to BoTorch to get their posterior means

# Make a grid across all budget parameters
# If justVDP, then picBudget is fixed and the grid is over the VDP budgets
# The picBudget parameter is only relevant if justVDP=T
# Remove budgets that do not meet bound or ratio constraints
# Reduce the number of budgets by clustering: one budget per cluster
# Because this is for BoTorch, delete the last budget
makeGridForTernary <- function(picBudget=0.3924, justVDP=F, nSteps=20){
  source(here::here("code", "breedSchemeInitialize.R"))

  # Set the lower and upper bounds
  # budget_constraints set in breedSchemeInitialize.R
  lower <- numeric(n_stages+1)
  lower[1] <- budget_constraints[1]
  lower[n_stages+1] <- budget_constraints[2]
  ratios <- budget_constraints[-(1:2)]
  for (stg in n_stages:2){
    lower[stg] <- ratios[stg]*lower[stg+1]
  }
  upper <- sapply(1:(n_stages+1), function(stg) 1 - sum(lower[-stg]))
  if (justVDP) lower[1] <- upper[1] <- picBudget

  # Function to determine if a budget fails bound and ratio constraints
  outConstraint <- function(budg){
    bad <- any(sapply(1:(n_stages+1), function(stg) budg[stg]<lower[stg])) |
      any(sapply(1:(n_stages+1), function(stg) budg[stg]>upper[stg]))
    for (stg in 1:n_stages){
      bad <- bad | budg[stg] / budg[stg+1] < ratios[stg]
    }
    return(bad)
  }

  # Make grid
  bsd$minPercentage <- lower
  bsd$maxPercentage <- upper
  bsd$percentageStep <- (upper - lower) / nSteps
  percGrid <- BreedSimCost::makeGrid(bsd, justVDP)
  percMat <- t(matrix(unlist(percGrid), nrow=length(percGrid[[1]])))
  if (justVDP){
    percMat <- percMat * (1 - picBudget)
    percMat <- cbind(picBudget, percMat)
  }

  # Eliminate budgets that do not satisfy constraints
  badBudg <- apply(percMat, 1, outConstraint)
  percMat <- percMat[!badBudg,]

  # To reduce the number of evaluation points, use kmeans clustering
  km <- kmeans(x=percMat, centers=200, iter.max=40)
  kmCtr <- km$centers

  # Send the grid points to BoTorch
  budgets <- kmCtr[,-4]
  saveRDS(budgets, here::here("data", "predPoints.rds"))
  return(budgets)
}

# Call python / BoTorch to get the posterior means
getPost <- function(n_iter, init_num, pred_points=0){
  print(n_iter)
  outFile <- paste0("pyOut", init_num, "_", n_iter, ".txt")
  callPyCommand <- paste("python getPosterior.py", init_num, n_stages,
                         n_iter, pred_points, ">", outFile)
  system(callPyCommand, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, timeout = 0)

  fileNotFound <- TRUE
  while(fileNotFound){
    tst <- try({
      dummy <- Sys.sleep(1)
      post <- readRDS(
        here::here("output",
                   paste0("posteriors", init_num, "_", n_iter, ".rds")))
    }, silent=T)
    fileNotFound <- "try-error" %in% class(tst)
    if (fileNotFound){
      cat(".")
    }
    else{
      cat("\n")
    }
  }
  return(post)
}

# Function to return values for the contour plot
# as a function of budget allocation
#' @param a, b, and c are all vectors
library(Ternary)
ternaryValFunc <- function(a, b, c, justVDP=F){
  knn <- 1
  xyz <- postGrid$predBudgets
  if (justVDP){
    xyz <- cbind(xyz[,2:3], 1 - rowSums(xyz[,1:3]))
    xyz <- xyz / rowSums(xyz)
  } else{
    xyz <- cbind(xyz[,1:2], 1 - rowSums(xyz[,1:2]))
  }
  abc <- cbind(a, b, c)
  oneABC <- function(abc){
    xyzDist <- apply(xyz, 1, function(v) sum(abs(v - abc)))
    return(mean(postGrid$predGains[order(xyzDist)[1:knn]]))
  }
  return(apply(abc, 1, oneABC))
}
