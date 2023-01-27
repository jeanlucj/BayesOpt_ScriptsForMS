require(here)
here::i_am("code/postGrid.R")

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
                         n_iter, pred_points, ">", outFile, "&")
  system(callPyCommand, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, timeout = 0) # wait=TRUE isn't working here
  # because the callPyCommand has & at the end...

  fileNotFound <- TRUE
  while(fileNotFound){
    tst <- try({
      dummy <- Sys.sleep(1)
      post <- readRDS(paste0("../output/posteriors",
                             init_num, "_", n_iter, ".rds"))
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
###################################################
# Everything from here on out uses makeGridForTernary and getPost
###################################################
# Default plot margins
# par(mar=c(5.1, 4.1, 4.1, 2.1))
justVDP <- FALSE
makeGridForTernary(justVDP=justVDP)
# These are the random seeds that came out to initiate six optimizations
# For some reason for the long optimization, some initializations died
allInit <- c(151116416, 46109957, 725670485, 365248867, 632547493, 877775392)
longInit <- c(151116416, 725670485, 365248867, 877775392)
# Plot out the predictions
plotOut <- FALSE

bigShort <- TRUE
bigLong <- TRUE
if (bigShort){
  if (bigLong) ai <- longInit else ai <- allInit
  if (bigLong) fileName <- "gridPostBigLong" else
    fileName <- "gridPostBigShort"
  if (bigLong) folder <- "bigLongOut" else folder <- "bigShortOut"
  if (justVDP) vdpTxt <- "VDP" else vdpTxt <- ""

  for (n_iter in 36+0:15*50){
    # Generate the predictions from BoTorch
    for (init_num in ai){
      print(paste("***********", init_num, "***********"))
      postGrid <- getPost(n_iter, init_num, pred_points="predPoints.rds")
      saveRDS(postGrid, paste0("../output/", folder, "/", fileName, vdpTxt, init_num, "_", n_iter, ".rds"))
    }
  }

  # Just some tests
  # tst <- readRDS(paste0("../output/", folder, "/", fileName, vdpTxt, init_num, "_", n_iter, ".rds"))
  # t2 <- readRDS(paste0("../output/", folder, "/posteriors", vdpTxt, init_num, "_", n_iter, ".rds"))

  if (plotOut){
    for (init_num in ai){
      postGrid <- readRDS(paste0("../output/", folder, "/", fileName, init_num, "_786.rds"))
      if (!justVDP){
        maxPost <- postGrid$bestBudget[1,]
        maxPost[3] <- 1 - sum(maxPost[1:2])
        values <- TernaryPointValues(ternaryValFunc, resolution = 24L)
      } else{
        values <- TernaryPointValues(ternaryValFunc, resolution = 24L, justVDP=T)
      }

      plotName <- paste0("../output/contour", fileName, vdpTxt, init_num, ".pdf")

      pdf(plotName)
      par(mar = rep(0.2, 4))
      if (!justVDP){
        TernaryPlot(alab = 'PIC', blab = 'SDN', clab = 'CET+PYT')
      } else{
        TernaryPlot(alab = 'SDN', blab = 'CET', clab = 'PYT')
      }
      ColourTernary(values)
      # AddToTernary(points, maxPost, pch = 16, cex = 2.2, col="red")
      # TernaryContour(ternaryValFunc, resolution = 36L)
      dev.off()
    }
  }#END plotOut

  # Make plots that use the analytical model of the selectiongain package
  makeSelectionGainPlots <- FALSE
  if (makeSelectionGainPlots){
    source("breedSchemeInitialize.R")
    # Make selGainGrid
    postGrid <- readRDS(paste0("../output/gridPostBigShort",
                               vdpTxt, ai[1], ".rds"))
    getNEntries <- function(percentages){
      perc <- c(percentages, 1 - sum(percentages))
      return(budgetToScheme(perc, bsd)$nEntries)
    }
    nEntryMat <- apply(postGrid$predBudgets, 1, getNEntries) %>% t
    vdpCovMat <- calcVDPcovMat(bsd)
    postGrid$predGains <- apply(nEntryMat, 1, calcVDPGain,
                    nFinal=bsd$nToMarketingDept, cov=vdpCovMat)

    values <- TernaryPointValues(ternaryValFunc, resolution = 24L, justVDP=T)
    plotName <- paste0("../output/contourSelGainBigShort", vdpTxt, ".pdf")

    pdf(plotName)
    par(mar = rep(0.2, 4))
    if (!justVDP){
      TernaryPlot(alab = 'PIC', blab = 'SDN', clab = 'CET+PYT')
    } else{
      TernaryPlot(alab = 'SDN', blab = 'CET', clab = 'PYT')
    }
    ColourTernary(values)
    # AddToTernary(points, maxPost, pch = 16, cex = 2.2, col="red")
    # TernaryContour(ternaryValFunc, resolution = 36L)
    dev.off()
  }#END makeSelGainPlots

  getMeanPICbudget <- FALSE
  if (getMeanPICbudget){
    allMaxBS <- NULL
    for (init_num in allInit){
      postGrid <- readRDS(paste0("../output/bigShortOut/gridPostBigShort", init_num, "_786.rds"))
      maxPost <- postGrid$bestBudget
      allMaxBS <- rbind(allMaxBS, maxPost)
    }
    print("BigShort")
    print(c(colMeans(allMaxBS), 1 - sum(colMeans(allMaxBS))))
    print(c(apply(allMaxBS, 2, sd), sd(1 - rowSums(allMaxBS)))/sqrt(6))
    # meanPIC budget across optimizations is 0.3924

    allMaxBL <- NULL
    for (init_num in longInit){
      postGrid <- readRDS(paste0("../output/bigLongOut/gridPostBigLong", init_num, "_786.rds"))
      maxPost <- postGrid$bestBudget
      allMaxBL <- rbind(allMaxBL, maxPost)
    }
    print("BigLong")
    print(c(colMeans(allMaxBL), 1 - sum(colMeans(allMaxBL))))
    print(c(apply(allMaxBL, 2, sd), sd(1 - rowSums(allMaxBL)))/sqrt(4))
  }
}#END if bigShort

smlShortManyIter <- !bigShort
if (smlShortManyIter){
  for (init_num in 1:14){
    for (n_iter in 36+0:30*50){ # Oops: didn't save for n_iter
      print(paste("***********", init_num, n_iter, "***********"))
      postGrid <- getPost(n_iter, init_num, pred_points="predPoints.rds")
      saveRDS(postGrid,
              paste0("../output/gridPostSmlShortManyIter", init_num, ".rds"))
    }
  }
}

if (plotOut){
  if (smlShortManyIter){
    allPostSm <- readRDS("../output/allPostSmlShortManyIter.rds")
    colnames(allPostSm) <- c("initNum", "nIter", "PIC", "SDN", "CET", "maxPredGain", "postVarAtMax")
    allPostSm <- allPostSm %>% as_tibble(allPostSm)
    allPostSm <- allPostSm %>% filter(!(initNum == 7 | initNum == 14)) %>%
      dplyr::mutate(initNum=if_else(initNum > 6, initNum-1, initNum)) %>%
      dplyr::mutate(replication=(initNum-1) %% 2, realInit=(initNum-1) %/% 2) %>%
      dplyr::mutate(PYT=1 - (PIC+SDN+CET)) %>% relocate(PYT, .after=CET)

    for (init_num in 1:6){
      maxPost <- filter(allPostSm, initNum==init_num & nIter==1500) %>% transmute(PIC=PIC, SDN=SDN, CET_PYT=CET+PYT)

      postGrid <- readRDS(
        paste0("../output/gridPostSmlShortManyIter", init_num, ".rds"))
      values <- TernaryPointValues(ternaryValFunc, resolution = 24L)

      plotName <- paste0("../output/contourPostSmlShortManyIterMax", init_num, ".pdf")

      pdf(plotName)
      par(mar = rep(0.2, 4))
      TernaryPlot(alab = 'PIC', blab = 'SDN', clab = 'CET+PYT')
      ColourTernary(values)
      AddToTernary(points, maxPost, pch = 16, cex = 2.2, col="red")
      # TernaryContour(ternaryValFunc, resolution = 36L)
      dev.off()
    }
  }#END smlShortManyItr
}#END if plotOut

# Calculate distances between predictions
if (FALSE){
  allGrid <- tibble()
  # Assemble all of the datasets
  for (init_num in 1:14){
    for (n_iter in 36+0:30*50){
      pg <- readRDS(paste0("../output/smlShortManyIter/posteriors",
                           init_num, "_", n_iter,".rds"))
      allGrid <- bind_rows(allGrid,
                           tibble(init=init_num, n_iter=n_iter,
                                  predBudgets=pg$predBudgets,
                                  predGains=pg$predGains,
                                  predVars=pg$predVars))
    }
  }

  distSameInit <- function(n_iter){
    ni <- n_iter
    distSame <- function(i){
      si <- ifelse (i < 7, (i-1)%/%2*2 + i%%2, i%/%2*2 - i%%2) + 1
      agi <- filter(allGrid, init==i & n_iter==ni) %>% pull(predGains) %>%
        scale
      ags <- filter(allGrid, init==si & n_iter==ni) %>% pull(predGains) %>%
        scale
      return(sqrt(sum((agi - ags)^2)))
    }
    return(mean(sapply(c(1:6, 8:13), distSame)))
  }
  dsi <- sapply(36+0:30*50, distSameInit)
  plot(dsi)

  distDiffInit <- function(n_iter){
    ni <- n_iter
    distDiff <- function(i){
      di <- setdiff(c(1:6, 8:13),
              c(i, ifelse (i < 7, (i-1)%/%2*2 + i%%2, i%/%2*2 - i%%2) + 1))
      agi <- filter(allGrid, init==i & n_iter==ni) %>% pull(predGains) %>%
        scale
      dd <- function(d, agi){
        agd <- filter(allGrid, init==d & n_iter==ni) %>% pull(predGains) %>%
          scale
        return(sqrt(sum((agi - agd)^2)))
      }
      return(mean(sapply(di, dd, agi=agi)))
    }
    return(mean(sapply(c(1:6, 8:13), distDiff)))
  }
  ddi <- sapply(36+0:30*50, distDiffInit)
  plot(ddi)

  yLim <- range(c(dsi, ddi))
  pdf("../output/HowManyIterations.pdf")
  plot(36+0:30*50, dsi, pch=16, cex=1.3, ylim=yLim,
       xlab="Optimization iteration", ylab="Gain landscape difference",
       cex.axis=1.3, cex.lab=1.3)
  points(36+0:30*50, ddi, pch=16, cex=1.3, col="blue")
  legend(250, 15, c("Different founders", "Same founders"), pch=16, col=c("blue", "black"), bty="n", cex=1.5, pt.cex=1.5)
  dev.off()
}#END FALSE
