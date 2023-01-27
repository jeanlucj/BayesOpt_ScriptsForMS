setwd("~/Documents/GitRepo/BayesOpt_Budgets/output")
allInit <- c(151116416, 46109957, 725670485, 365248867, 632547493, 877775392)
longInit <- c(151116416, 725670485, 365248867, 877775392)

allMaxBS <- NULL
for (init_num in allInit){
  postGrid <- readRDS(paste0("./bigShortOut/gridPostBigShort", init_num, "_786.rds"))
  maxPost <- postGrid$bestBudget
  allMaxBS <- rbind(allMaxBS, maxPost)
}
print("BigShort")
print(c(colMeans(allMaxBS), 1 - sum(colMeans(allMaxBS))))
# 0.39239415 0.38372097 0.19763981 0.02624507
print(c(apply(allMaxBS, 2, sd), sd(1 - rowSums(allMaxBS)))/sqrt(6))
# 0.03592147 0.03691799 0.04124302 0.00866042
# meanPIC budget across optimizations is 0.3924

allMaxBL <- NULL
for (init_num in longInit){
  postGrid <- readRDS(paste0("./bigLongOut/gridPostBigLong", init_num, "_786.rds"))
  maxPost <- postGrid$bestBudget
  allMaxBL <- rbind(allMaxBL, maxPost)
}
print("BigLong")
print(c(colMeans(allMaxBL), 1 - sum(colMeans(allMaxBL))))
# 0.32610319 0.38645675 0.25560986 0.03183021
print(c(apply(allMaxBL, 2, sd), sd(1 - rowSums(allMaxBL)))/sqrt(4))
#  0.02490447 0.06637956 0.08262343 0.01424555

# Significant differences?
for (i in 1:4){
  print(t.test(allMaxBS[, i], allMaxBL[, i])$p.value)
}
