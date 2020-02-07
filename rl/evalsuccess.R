library(tidyverse)

s <- tibble()

for (f in dir(pattern = 'successrate.*.rdata', recursive = T)) {
  load(f)
  print(succ)
  s <- bind_rows(s, succ)
}

s <- s %>% filter(!is.na(lwiter)) %>%
  mutate(confg=paste0(numOrig,'_',numSigs)) %>%
  mutate(numSigs.f = as.factor(numSigs)) %>%
  mutate(numOrig.f = as.factor(numOrig))


for (i in unique(s$numOrig)) {
  print(paste0("Correlation for numOrig: ",i))
  s.lm <- lm(lwc~., s %>% filter(numOrig==i) %>%
               mutate(AbsDiff=abs(numSigs-numOrig)) %>%
               select(AbsDiff,lwc) %>% mutate_all(funs(scale)))
  print(summary(s.lm)$coefficients)
}

for (i in unique(s$numOrig)) {
  print(paste0("Correlation for numOrig: ",i))
  s.lm <- lm(r.mean~., s %>% filter(numOrig==i) %>%
               mutate(AbsDiff=abs(numSigs-numOrig)) %>%
               select(AbsDiff,r.mean) %>% mutate_all(funs(scale)))
  print(summary(s.lm)$coefficients)
}

for (i in unique(s$numOrig)) {
  print(paste0("Correlation for numOrig: ",i))
  s.lm <- lm(r.centerMean~., s %>% filter(numOrig==i) %>%
               mutate(AbsDiff=abs(numSigs-numOrig)) %>%
               select(AbsDiff,r.centerMean) %>% mutate_all(funs(scale)))
  print(summary(s.lm)$coefficients)
}

for (i in unique(s$numOrig)) {
  print(paste0("Correlation for numOrig: ",i))
  s.lm <- lm(g.sumMean~., s %>% filter(numOrig==i) %>%
               mutate(AbsDiff=abs(numSigs-numOrig)) %>%
               select(AbsDiff,g.sumMean) %>% mutate_all(funs(scale)))
  print(summary(s.lm)$coefficients)
}

for (i in unique(s$numOrig)) {
  print(paste0("Correlation for numOrig: ",i))
  s.lm <- lm(g.sdMean~., s %>% filter(numOrig==i) %>%
               mutate(AbsDiff=abs(numSigs-numOrig)) %>%
               select(AbsDiff,g.sdMean) %>% mutate_all(funs(scale)))
  print(summary(s.lm)$coefficients)
}

ggplot(s)+geom_boxplot(aes(numSigs.f,(lwiter-lwtin)/20, color=numOrig.f))+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Iterations')
ggsave('success_Iterations.png')

ggplot(s)+geom_boxplot(aes(numSigs.f,r.mean, color=numOrig.f))+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Mean(R)')
ggsave('success_Mean_R.png')

ggplot(s)+geom_boxplot(aes(numSigs.f,r.centerMean, color=numOrig.f))+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Center Mean(R)')
ggsave('success_MeanCenter_R.png')

ggplot(s)+geom_boxplot(aes(numSigs.f,g.sumMean, color=numOrig.f))+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Mean(G)')
ggsave('success_Mean_G.png')

ggplot(s %>% filter(g.sdMean<100))+geom_boxplot(aes(numSigs.f,g.sdMean, color=numOrig.f))+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Std.Deviation(G)')
ggsave('success_Sd_G.png')

