library(ggplot2)

bigMeans$nCyc <- as.factor(bigMeans$nCyc)
bigMeans$nPICcycles <- as.factor(bigMeans$nPICcycles)
bmL <- bigMeans %>% select(-meanGain, -meanVar, -contains("sd", ignore.case=F)) %>%
  tidyr::pivot_longer(cols=starts_with("mean"), names_to="stage",
                      values_to="allocation") %>%
  mutate(stage=gsub("mean", "", stage))
bmLsd <- bigMeans %>% select(-sdGain, -sdVar, -contains("mean")) %>%
  tidyr::pivot_longer(cols=starts_with("sd"), names_to="stageSD",
                      values_to="allocSD")
bmL <- bmL %>% bind_cols(allocSE=bmLsd$allocSD / sqrt(bmLsd$nOpt))

# I think I want to have separate graphs for each scheme
allocPlot <- function(keepCyc="06", keepPIC=1){
  p <- ggplot(bmL %>% filter(nCyc==keepCyc & nPICcycles==keepPIC),
              aes(x=nOptIter, y=allocation))
  return(p + geom_pointrange(aes(color=stage,
                          ymin=allocation-allocSE, ymax=allocation+allocSE)))
}

p <- allocPlot("06", 1)
p <- allocPlot("12", 1)
p <- allocPlot("06", 2)
p <- allocPlot("12", 2)

# Small Many Iter here
allSmlbyFnd <- allSml %>% group_by(subSeed, nOptIter) %>%
  summarize(meanPIC=mean(budgPIC),
            meanSPN=mean(budgSDN),
            meanCET=mean(budgCET),
            meanPYT=mean(budgPYT)) %>%
  tidyr::pivot_longer(cols=starts_with("mean"), names_to="stage",
                      names_prefix="mean", values_to="allocation")

smlMeans <- allSmlbyFnd %>% group_by(nOptIter, stage) %>%
  summarize(allocMean=mean(allocation), allocSE=sd(allocation)/4)

p <- ggplot(smlMeans,
            aes(x=nOptIter, y=allocMean))
p + geom_pointrange(aes(color=stage,
                               ymin=allocMean-allocSE, ymax=allocMean+allocSE))

# These plots show (I think) that after ~ 250 iterations there is relative
# stability across founders, with occasional perturbations due to the flat
# function surface close to the maximum
q <- ggplot(allSmlbyFnd %>% filter(stage=="SPN"),
            aes(x=nOptIter, y=allocation))
q + geom_line(aes(color=subSeed), show.legend=F)
q <- ggplot(allSmlbyFnd %>% filter(stage=="CET"),
            aes(x=nOptIter, y=allocation))
q + geom_line(aes(color=subSeed), show.legend=F)
q <- ggplot(allSmlbyFnd %>% filter(stage=="PYT"),
            aes(x=nOptIter, y=allocation))
q + geom_line(aes(color=subSeed), show.legend=F)
q <- ggplot(allSmlbyFnd %>% filter(stage=="PIC"),
            aes(x=nOptIter, y=allocation))
q + geom_line(aes(color=subSeed), show.legend=F)

# OK.  I need to figure out again how to make the overall surface plots
