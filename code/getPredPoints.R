# If there are a number of points at which you want to make a prediction
# get those points with this script
require(here)
here::i_am("code/getPredPoints.R")

if (!exists("predFile")) predFile <- "predFile.rds"
print(paste("Getting budgets for which to predict gain", predFile))

predBudgets <- readRDS(here::here("data", predFile))
print(predBudgets)
