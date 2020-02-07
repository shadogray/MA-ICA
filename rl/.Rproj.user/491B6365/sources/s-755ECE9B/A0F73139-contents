library(tidyverse)
options(keep.source = T)

corrFile <- 'correlations.rdata'
runCorrs <- NULL
allCorrs <- NULL
allX <- NULL

if (!file.exists(corrFile)) {
  tryCatchLog::tryCatchLog({
    dDirs <- list.dirs('.', recursive = F) %>% str_subset('ssteps\\d$')
    for (dDir in dDirs) {
      rdatas <- dir(path = dDir, pattern = 'finDemo.env_.*_Res.rdata', recursive = T)
      for (rdata in rdatas) {

        rdfile <- paste0(dDir,'/',rdata)
        if (file.exists(rdfile)) {
          print(paste('loading: ',rdfile))
          load(rdfile)
          if (!has_name(result, 'myEnvs')) {
            print(paste('old results format, ignoring', rdfile))
            next()
          }

          mev <- result$myEnv
          envRunCorr <- result$runCorr
          envRunCorr$numSigs <- mev$numSigs
          envRunCorr$numOrig <- mev$numOrig
          if (is.null(runCorrs)) {
            runCorrs <- envRunCorr
          } else {
            runCorrs <- runCorrs %>% bind_rows(envRunCorr)
          }
          print(runCorrs, n=Inf)

          for (run in 1:length(result$myEnvs)) {
            print(paste('evaluate run',run))

            mev <- result$myEnvs[[run]]

            envRcs <- as_tibble(result$runCorrs[[run]])
            envRcs$run <- run
            envRcs$numSigs <- mev$numSigs
            envRcs$numOrig <- mev$numOrig
            if (is.null(allCorrs)) {
              allCorrs <- envRcs
            } else {
              allCorrs <- allCorrs %>% bind_rows(envRcs)
            }

            X <- mev$X
            X$run <- run
            X$numSigs <- mev$numSigs
            X$numOrig <- mev$numOrig

            X.orig <- mev$X.orig %>% rowid_to_column('xid')
            X.orig$run <- run
            X.orig$numSigs <- mev$numSigs
            X.orig$numOrig <- mev$numOrig

            ica.worstX <- mev$ica.worstX %>% rowid_to_column('xid')
            ica.worstX$run <- run
            ica.worstX$numSigs <- mev$numSigs
            ica.worstX$numOrig <- mev$numOrig

            if (is.null(allX)) {
              allX <- X
              allX.orig <- X.orig
              allica.worstX <- ica.worstX
            } else {
              allX <- allX %>% bind_rows(X)
              allX.orig <- allX.orig %>% bind_rows(X.orig)
              allica.worstX <- allica.worstX %>% bind_rows(X.orig)
            }
          }
        } else {
          print(paste('Cannot load: ',rdfile))
        }
      }
    }

    save(runCorrs, allCorrs, allX, allX.orig, file = corrFile)
  })
} else {
  load(corrFile)
}

runCorrs.all <- runCorrs

for (i in unique(runCorrs.all$numOrig)) {
  print(paste0("Correlation for NumOrig=",i,":"))
  d <- runCorrs.all %>% mutate(AbsDiff=abs(numSigs-numOrig)) %>%
    filter(numOrig==i) %>%
    select(Si,AbsDiff) %>% mutate_all(funs(scale))
  d.lm <- lm(Si~., d)
  #print(d.lm)
  #print(summary(d.lm))
  d.lm.s <- summary(d.lm)
  print(d.lm.s$coefficients)
}

for (t in c('mad','chisq','z','t','iqr')) {
  runCorrs.all <- runCorrs.all %>% group_by(numOrig,numSigs) %>%
    mutate(!!paste0('Si.',t) := outliers::scores(Si, type=t)) %>%
    mutate(!!paste0('Si.',t,'.prob') := outliers::scores(Si, type=t, prob=T)) %>%
    ungroup()
}

#rCoutliers <- runCorrs.all %>% group_by(numSigs,numOrig) %>% summarize(Si = outliers::outlier(Si))
#runCorrs.outliers <- runCorrs.all %>% inner_join(rCoutliers) %>% arrange(numSigs,numOrig)
runCorrs.outliers <- runCorrs.all %>% filter(Si.mad.prob>.9) %>% select(run,numSigs,numOrig)
print('RunCorrelations - Outliers:')
print(runCorrs.outliers, n=Inf)

runCorrs <- runCorrs.all %>% anti_join(runCorrs.outliers)
runCorrs.mean <- runCorrs %>% group_by(numOrig,numSigs) %>%
  summarize(runs = n(),
            Si.med = median(Si), Si.mean = mean(Si), Si.var = var(Si),
            Si.q1 = quantile(Si, seq(0,1,.1))[1],
            Si.q10 = quantile(Si, seq(0,1,.1))[10]) %>%
  arrange(desc(Si.mean)) %>%
  mutate(Si = Si.mean) %>%
  rowid_to_column('rid')
print('RunCorrelations - Mean/Median/Variance:')
print(runCorrs.mean, n=Inf)

runCorrs.opt <- runCorrs.mean %>% filter(
  rid %in% (runCorrs.mean %>% group_by(numOrig) %>% summarise(mrid = min(rid)))$mrid) %>% ungroup()
print('RunCorrelations - Optimum:')
print(runCorrs.opt, n=Inf)

data.all <- runCorrs.all %>% # filter(nlf != 'NONE') %>%
  mutate(numSigs.f = as.factor(numSigs)) %>%
  mutate(numOrig.f = as.factor(numOrig))

data.outliers <- data.all %>% inner_join(runCorrs.outliers)

data <- runCorrs %>% # filter(nlf != 'NONE') %>%
  mutate(numSigs.f = as.factor(numSigs)) %>%
  mutate(numOrig.f = as.factor(numOrig))

data.sumr <- runCorrs %>% # filter(nlf != 'NONE') %>%
  select(numSigs,numOrig,Si) %>%
  group_by(numSigs,numOrig) %>% summarise(Si.mean = mean(Si), Si.max = max(Si)) %>% ungroup() %>%
  mutate(numSigs.f = as.factor(numSigs)) %>%
  mutate(numOrig.f = as.factor(numOrig))

data.sumr %>% ggplot(aes(numSigs,Si.max, color=numOrig.f))+
  geom_point()+
  geom_line()+
  geom_point(data=data.sumr %>% filter(numSigs==numOrig), aes(numSigs,Si.max, color=numOrig.f), size=3)+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Max. Correlation')
ggsave('AllCorrelation_Summary_Max.png', dpi = 'print', width = 15, height = 10, units = 'cm')

ggplot(data.sumr) + geom_point()

data.sumr %>% ggplot(aes(numSigs,Si.mean, color=numOrig.f))+
  geom_point()+
  geom_line()+
  geom_point(data=data.sumr %>% filter(numSigs==numOrig), aes(numSigs,Si.mean, color=numOrig.f), size=3)+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Mean Correlation')
ggsave('AllCorrelation_Summary_Mean.png', dpi = 'print', width = 15, height = 10, units = 'cm')

data.sumr %>% ggplot(aes(numSigs,Si.mean, color=numOrig.f))+
  geom_point()+
  geom_line()+
  geom_point(data=data.sumr %>% filter(numSigs==numOrig), aes(numSigs,Si.mean, color=numOrig.f), size=3)+
  geom_jitter(data=data.outliers, aes(numSigs,Si, color=numOrig.f), shape='asterisk', size=.7, width=.02)+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Mean Correlation')
ggsave('AllCorrelation_Summary_Mean_Outliers.png', dpi = 'print', width = 15, height = 10, units = 'cm')

data.all %>% ggplot()+
  geom_boxplot(data = data.all %>% mutate(numSigs = as.factor(numSigs)),
               aes(numSigs,Si, color=numOrig.f), position = position_dodge2(width=.5, ), varwidth = T)+
  coord_cartesian(ylim = c(.6,1))+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Correlation')
ggsave('AllCorrelation_Summary_Box.png', dpi = 'print', width = 15, height = 10, units = 'cm')

data.all %>% ggplot()+
  geom_smooth(aes(numSigs,Si, color=numOrig.f), span=.5)+
  #geom_line(data=data.sumr, aes(numSigs,Si.mean, color=numOrig.f))+
  geom_point(data=data.sumr %>% filter(numSigs==numOrig), aes(numSigs,Si.mean, color=numOrig.f), size=3)+
  coord_cartesian(ylim = c(.6,1))+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Correlation')
ggsave('AllCorrelation_Summary_Smooth.png', dpi = 'print', width = 15, height = 10, units = 'cm')

data.all %>% ggplot()+
  geom_jitter(aes(numSigs,Si, color=numOrig.f), width = .1)+
  geom_line(data=data.sumr, aes(numSigs,Si.mean, color=numOrig.f))+
  geom_point(data=data.sumr %>% filter(numSigs==numOrig), aes(numSigs,Si.mean, color=numOrig.f), size=3)+
  coord_cartesian(ylim = c(.6,1))+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Correlation')
ggsave('AllCorrelation_Summary.png', dpi = 'print', width = 15, height = 10, units = 'cm')

