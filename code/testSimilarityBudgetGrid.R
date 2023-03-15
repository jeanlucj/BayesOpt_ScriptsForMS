af <- list.files(here::here("output", "PostGrid", "Feb17"))
std <- readRDS(here::here("output", "PostGrid", "Feb17", af[1]))$predBudgets
#std <- std[order(std[,1], std[,2], std[,3]),]
allSC <- NULL
for (fn in af){
  tst <- readRDS(here::here("output", "PostGrid", "Feb17", fn))$predBudgets
  #tst <- tst[order(tst[,1], tst[,2], tst[,3]),]
  sc <- sum(cor(std[,1], tst[,1]),
            cor(std[,2], tst[,2]),
            cor(std[,3], tst[,3]))
  allSC <- c(allSC, sc)
  #if (sc < 2.90) print(fn)
}
hist(allSC)
table(round(allSC, 3))
afSubSeed <- strsplit(af, "_", fixed=T)
afSubSeed <- sapply(afSubSeed, function(v) v[1]) %>% unique
writeLines(afSubSeed, "bigSeeds.txt")

###############################################
# Small breeding programs here: Small Many Long
af <- list.files(here::here("output", "PostGrid"), pattern="sml")
std <- readRDS(here::here("output", "PostGrid", af[1]))$predBudgets
std <- std[order(std[,1], std[,2], std[,3]),]
allSC <- NULL
for (fn in af){
  tst <- readRDS(here::here("output", "PostGrid", fn))$predBudgets
  tst <- tst[order(tst[,1], tst[,2], tst[,3]),]
  sc <- sum(cor(std[,1], tst[,1]),
            cor(std[,2], tst[,2]),
            cor(std[,3], tst[,3]))
  allSC <- c(allSC, sc)
}
hist(allSC)
table(round(allSC, 3))
afSubSeed <- strsplit(af, "_", fixed=T)
afSubSeed <- sapply(afSubSeed, function(v) v[1]) %>% unique
writeLines(afSubSeed, here::here("output", "PostGrid", "smlSeeds.txt"))
