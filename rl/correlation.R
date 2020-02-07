library(tidyverse)

k <- c(1.53,1.92,2.34,2.97,3.45,3.97)
m <-matrix(rep(k,6),6,6)
t(m)/k

# corr;ICA;numOrig;CorrOrig;comment
data <- read_delim('correlation.csv', delim = ';', locale = locale(decimal_mark = '.')) %>%
  extract(ICA, 'numSigs', regex = 'ICA(\\d)', remove = F) %>%
  mutate_at(vars(numSigs,corr,CorrOrig), as.numeric) %>%
  mutate(lag = grepl('Lag',comment)) %>%
  mutate(dif = grepl('Diff',comment)) %>%
  mutate_at(vars(ICA), as.factor) %>%
  mutate(numSigs.f = as.factor(numSigs)) %>%
  mutate(numOrig.f = as.factor(numOrig)) %>%
  mutate(nlf = as.factor(paste0(ifelse(dif,'DIFF','ORIG'),':', ifelse(lag,'LAG','NONE'))))

data %>% ggplot(aes(numSigs,corr, color=numOrig.f))+
  geom_point()+
  geom_line()+
  geom_point(data=data %>% filter(numSigs==numOrig), aes(numSigs,corr, color=numOrig.f), size=3)+
  facet_wrap(~nlf)+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Correlation')
ggsave('Correlation_Summary.png', dpi = 'print')

data.max <- data %>% # filter(nlf != 'NONE') %>%
  select(numSigs,numOrig,corr) %>%
  group_by(numSigs,numOrig) %>% summarise(corr.mean = mean(corr), corr.max = max(corr)) %>%
  mutate(numSigs.f = as.factor(numSigs)) %>%
  mutate(numOrig.f = as.factor(numOrig))

data.max %>% ggplot(aes(numSigs,corr.max, color=numOrig.f))+
  geom_point()+
  geom_line()+
  geom_point(data=data.max %>% filter(numSigs==numOrig), aes(numSigs,corr.max, color=numOrig.f), size=3)+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Max. Correlation')
ggsave('Correlation_Summary_Max.png', dpi = 'print')

data.none <- data %>% filter(!dif & !lag)
data.none %>% ggplot(aes(numSigs,corr, color=numOrig.f))+
  geom_point()+
  geom_line()+
  geom_point(data=data.none %>% filter(numSigs==numOrig), aes(numSigs,corr, color=numOrig.f), size=3)+
  coord_cartesian(ylim = c(.6,1))+
  #ggtitle('Correlation of Tracks with Target Signal')+
  scale_color_discrete(name = 'Original\nSignals')+
  xlab('Number of Signals of ICA')+
  ylab('Correlation')
ggsave('Correlation_Summary_None.png', dpi = 'print')

