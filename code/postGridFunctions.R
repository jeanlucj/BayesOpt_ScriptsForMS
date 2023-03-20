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
makeGridForTernary <- function(picBudget=0.3924, justVDP=F, nSteps=20, doMax=T){
  source(here::here("code", "breedSchemeInitialize.R"))
  if (!exists("n_stages")) n_stages <- bsd$nStages
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
  if (doMax){
    bsd$minPercentage <- lower
    bsd$maxPercentage <- upper
    bsd$percentageStep <- (upper - lower) / nSteps
  } else{
    bsd$percentageStep <- (bsd$maxPercentage - bsd$minPercentage) / nSteps
  }
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
# Just for the record, BreedSimCost::makeGrid
makeGrid <- function (bsd, justVDP = F){
  if (justVDP) strt <- 2 else strt <- 1
  percList <- as.list(seq(from = bsd$minPercentage[strt],
                          to = bsd$maxPercentage[strt],
                          by = bsd$percentageStep[strt]))
  for (stage in setdiff(1 + 1:bsd$nStages, strt)){
    newPerc <- seq(from = bsd$minPercentage[stage],
                   to = bsd$maxPercentage[stage],
                   by = bsd$percentageStep[stage])
    percList <- mapply(c, rep(percList, each = length(newPerc)),
                       rep(newPerc, length(percList)), SIMPLIFY = F)
  }
  budgets <- lapply(percList, function(v) return(v/sum(v)))
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
#' @param knn is how many nearest neighbors to use
library(Ternary)
ternaryValFunc <- function(a, b, c, justVDP=F, knn=1, postGrid=postGrid){
  xyz <- postGrid$predBudgets
  if (justVDP){
    xyz <- cbind(xyz[,2:3], 1 - rowSums(xyz[,1:3]))
    xyz <- xyz / rowSums(xyz)
  } else{
    xyz <- cbind(xyz[,1:2], 1 - rowSums(xyz[,1:2]))
  }
  abc <- cbind(a, b, c)
  oneABC <- function(abc){
    # Manhattan distance to the prediction points in xyz
    xyzDist <- apply(xyz, 1, function(v) sum(abs(v - abc)))
    # Weighted mean to a number of nearest neighbors
    tv <- weighted.mean(postGrid$predGains[order(xyzDist)[1:knn]],
                        1 / xyzDist[order(xyzDist)[1:knn]])
    return(tv)
  }
  return(apply(abc, 1, oneABC))
}

# Function: take a file name, get the data, draw a ternary plot, make contour
# one tenth of an SD below the max.  See postGrid.R for the original version
contourMinusSD <- function(initnum=188690362, repl=0, plotName=NULL, knn=3,
                           allContours=F){
  require(Ternary)

  # Analyzing small test again
  smlTA <- getTestAgain(picBudgString="0.36805244")
  varCandMeanSD <- dplyr::last(smlTA$parmsSE$varCandMeanse)*20 # 5.107004

  if (is.null(plotName)){
    replTxt <- repl
    if (length(repl) > 1) replTxt <- NULL
    plotName <- paste0(initnum, "_", replTxt, "_smlTernary.pdf")
  }
  pdf(here::here("output", "_PubOut", "figures", plotName),
      width=6*length(repl), height=6)
  par(mfcol=c(1, length(repl)))
  par(mar = rep(0.2, 4))
  for (i in repl){
    # Construct name 188690362_0_1036sml06PIC.rds
    fileName <- paste0(initnum, "_", i, "_1536sml06PIC.rds")
    pg <- readRDS(here::here("output", "PostGrid", fileName))
    maxPost <- pg$bestBudget[1,]
    maxPost[3] <- 1 - sum(maxPost[1:2])

    # Make a new function that encapsulates the necessary parameters
    tvf <- function(a, b, c) return(ternaryValFunc(a, b, c,
                                    justVDP=F, knn=knn, postGrid=pg))
    values <- Ternary::TernaryPointValues(tvf, resolution = 36L)
    zlim <- round(range(values["z",]), 2)
    Ternary::TernaryPlot(atip=paste(zlim, collapse=" "),
                         alab = 'PIC', blab = 'SDN', clab = 'CET+PYT',
                         lab.cex=1.3, axis.cex=1.3)
    Ternary::ColourTernary(values)
    Ternary::AddToTernary(points, maxPost, pch = 16, cex = 2.2, col="red")
    contourLevel <- pg$maxPredGain - varCandMeanSD/5
    if (allContours){
      Ternary::TernaryContour(tvf, resolution = 36L)
    } else{
      Ternary::TernaryContour(tvf, resolution = 36L, levels=contourLevel,
                              labels="")
    }
  }
  dev.off()
}

# This is a bit nuts, but I have to make a special function for the mean
# across replications if I am to be able to draw the contour
has499Iter <- read_csv(here::here("data", "has499Iter.csv"))
allSubSeed <- has499Iter %>% pull(subSeed) %>% unique

allPG <- list()
for (initnum in allSubSeed){
  for (i in 0:1){
    # Construct name 188690362_0_1036sml06PIC.rds
    fileName <- paste0(initnum, "_", i, "_1536sml06PIC.rds")
    pg <- readRDS(here::here("output", "PostGrid", fileName))
    allPG <- c(allPG, list(pg))
  }
}
allMaxPreGain <- sapply(allPG, function(pg) return(pg$maxPredGain))
maxPredGain <- mean(allMaxPredGain)
allBestBudget <- sapply(allPG, function(pg) return(pg$bestBudget))
allBestBudget <- t(rbind(allBestBudget[1:2,], 1 - colSums(allBestBudget[1:2,])))
bestBudget <- colMeans(allBestBudget)

tvfMean <- function(a, b, c, apg){
  knn <- 3
  abc <- cbind(a, b, c)
  onePG <- function(postGrid, abc){
    xyz <- postGrid$predBudgets
    xyz <- cbind(xyz[,1:2], 1 - rowSums(xyz[,1:2]))
    oneABC <- function(abc){
      # Manhattan distance to the prediction points in xyz
      xyzDist <- apply(xyz, 1, function(v) sum(abs(v - abc)))
      # Weighted mean to a number of nearest neighbors
      tv <- weighted.mean(postGrid$predGains[order(xyzDist)[1:knn]],
                          1 / xyzDist[order(xyzDist)[1:knn]])
      return(tv)
    }
    return(apply(abc, 1, oneABC))
  }
  return(rowMeans(sapply(apg, FUN=onePG, abc=abc)))
}

# Analyzing small test again
smlTA <- getTestAgain(picBudgString="0.36805244")
varCandMeanSD <- dplyr::last(smlTA$parmsSE$varCandMeanse)*20 # 5.107004

plotName <- "___smlTernary.pdf"
pdf(here::here("output", "_PubOut", "figures", plotName),
    width=6, height=6)
par(mar = rep(0.2, 4))
values <- Ternary::TernaryPointValues(tvfMean, resolution = 36L, apg=allPG)
# WARNING: what's commented out below has to be defined
Ternary::TernaryPlot(alab = 'PIC', blab = 'SDN', clab = 'CET+PYT',
                     lab.cex=1.3, axis.cex=1.3)
Ternary::ColourTernary(values)
Ternary::AddToTernary(points, bestBudget, pch = 16, cex = 2.2, col="red")
contourLevel <- maxPredGain - (4:6)*(varCandMeanSD/20)
tvf <- function(a, b, c) return(tvfMean(a, b, c, apg=allPG))
Ternary::TernaryContour(tvf, resolution = 36L,
                        levels=contourLevel, labels=c("0.20", "0.25", "0.3"))
dev.off()
