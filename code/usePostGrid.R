require(here)
here::i_am("code/usePostGrid.R")
# Functions to use here
source(here::here("code", "postGridFunctions.R"))
# Script to figure out what optimizations have been succesful...
source(here::here("code", "SortThroughOutputs.R"))

# Script to use getPosterior.py to find the maximum and the expected posterior
# gain at different iterations
n_stages <- 3

# Options are
# simSize "big" or "sml"
# nCycle "6" or "12"
# If simSize=="big" then nOptIter <- 99 else nOptIter <- 499
simSize <- "big"
justVDP <- FALSE

args <- commandArgs(trailingOnly=T)
whereSimSize <- grep("simSize", args)
if (length(whereSimSize) > 0) simSize <- args[whereSimSize[1]+1]
whereVDP <- grep("justVDP", args)
if (length(whereVDP) > 0) justVDP <- as.logical(args[whereVDP[1]+1])

if (simSize == "big"){
  simDF <- has99Iter
  nOptIter <- 99
  iterInterval <- 50
} else{
  simDF <- has499Iter
  nOptIter <- 499
  iterInterval <- 100
}
allInit_num <- dplyr::pull(simDF, fileBoTorch)
allInit_num <- gsub("BoTorchOut", "", allInit_num)
allInit_num <- gsub(".rds", "", allInit_num)

# This will save predPoints.rds
# NOTE: this function is slow
init_num <- simDF %>% filter(hasDirectory) %>% slice(1) %>% pull(subSeed)
makeGridForTernary(justVDP=justVDP)

if (simSize == "big") sizeTxt <- "big" else sizeTxt <- "sml"
if (justVDP) vdpTxt <- "VDP" else vdpTxt <- "PIC"

# Will put these files in a separate directory, so make it if it doesn't exist
outputFiles <- list.files(here::here("output"))
if (length(grep("PostGrid", outputFiles)) == 0){
  dir.create(here::here("output", "PostGrid"))
}

for (n_iter in 36+seq(from=0, to=(nOptIter+1)*3, by=iterInterval)){
  # Generate the predictions from BoTorch
  for (init_num in allInit_num){
    print(paste("***********", init_num, "***********"))
    postGrid <- getPost(n_iter, init_num, pred_points="predPoints.rds")

    cycTxt <- simDF %>% slice(grep(init_num, fileBoTorch)) %>% pull(nCycNBO)
    if (nchar(cycTxt) == 1) cycTxt <- paste0("0", cycTxt)
    saveRDS(postGrid, here::here("output", "PostGrid",
      paste0(init_num, "_", n_iter, sizeTxt, cycTxt, vdpTxt,  ".rds")))
  }
}
