library(tidyverse)
#library(timeSeries)
#library(corrr)
library(lubridate)
#library(directlabels)

source('finlearn.R')

sig.v <- NULL
orig.s <- NULL
test <- F
evalPlot <- T
evalPrint <- T
local <- F
args <- commandArgs(trailingOnly = T)
if (length(args) > 0) sig.v <- unlist(strsplit(args[1],','))
if (length(args) > 1) orig.s <- unlist(strsplit(args[2],','))
if (length(args) > 2) {
  test <- str_detect(args[3],'T')
  local <- str_detect(args[3],'L')
  if (str_detect(args[3],'[pP]')) evalPlot <- !str_detect(args[3],'p')
  if (str_detect(args[3],'[rR]')) evalPrint <- !str_detect(args[3],'r')
}

node <- Sys.info()['nodename']
nlog = T
deblog = F
if (deblog) {
  nlog = T
}
options(tibble.width = Inf)
options(width = 2000)
#options(error = function(){traceback(1);print();stop()})

tailRuns = 200
numIter = 100000
deltaSumR = 2*.001
epsilon = .0
alpha = .1
gamma = .5
ctrl.types = factor('Q','SARSA')

type = 'SARSA'
#args <- commandArgs(trailingOnly = T)
#if (length(args) > 0) type = args[1]
## load('test.rdata');res=result;env=result$myEnv;tracks=result$tracks
#######

getTrack <- function(env = myEnv) {
  env$SumR <- 0
  findRes <- findWay(env)
  b <- tibble(.rows = length(env$B1))
  for (n in names(env)[grep('B\\d',names(env))]) b <- bind_cols(b, env[n])
  findRes$b <- b
  track <- findRes$track
  if (!has_name(track,'t')) track <- track %>% rowid_to_column('t')
  findRes$track <- track
  findRes
}

printTrack <- function(res, env = res$myEnv, track = res$track) {
  print(track, n=nrow(track))
  print(paste('Type:',type,
              'Iterations:',env$Iterations,
              'FinalR:', track %>% filter(res.terminal) %>% select(res.R),
              'SumR:',toString(env$SumR),'W:',toString(env$W)))
}

plotResults <- function(res, env = res$myEnv, track = res$track) {
  print(paste0('processing result: ',env$numSigs,'/',env$numOrig,': plots'))

  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,': ')
  fileSfx <- paste0('S_density_ica')
  fileExt <- paste0('_',env$numSigs,'_',env$numOrig,'.png')

  logs <- bind_rows(map(res$myEnvs,'log'))
  comb <- gtools::combinations(env$xColNum-1,2)
  for (i in 1:nrow(comb)) {
    si <- paste0('SVal.S',comb[i,])
    #ggplot(env$S %>% select(si,n,sig.p) %>% group_by_at(c(si,'sig.p')) %>% summarize(n=sum(n))) +
    #  geom_contour(aes_string(si[1],si[2], z='n'))+
    #  facet_wrap(~sig.p)+theme(strip.text.x = element_text(size = 5)) +
    #ggplot(env$S %>% select(si,n) %>% group_by_at(si) %>% summarize(n=sum(n))) +
    #  geom_contour(aes_string(si[1],si[2], z='n'))+
    #  theme(strip.text.x = element_text(size = 5)) +
    #  ggtitle(paste0(titleSfx,'State-Evaluation Density - States ',si[2],'~',si[1]))
    lab.si <- paste0('S',comb[i,])
    ggplot(logs %>% select(si) %>% mutate_all(list(~round(.))), mapping=aes_string(si[1],si[2])) +
      geom_density_2d() +
      theme(strip.text.x = element_text(size = 5)) +
      #ggtitle(paste0(titleSfx,'State-Evaluation Density - States ',si[2],'~',si[1]))
      xlab(lab.si[1]) + ylab(lab.si[2])

    ggsave(paste0(fileSfx,'_',i,fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  }

  b <- as_tibble(as.matrix(res$b)*env$baseX)
  w.s <- colMeans(env$log %>% tail(tailRuns) %>% select(matches('W\\.S.')))
  W <- env$W / sum(abs(env$W))
  fs <- env$FS %>% rowid_to_column('t')
  for (n in env$xColNames) fs[paste0(n,'.w')] <- w.s[paste0('W.',n)]*env$X[n]
  fs$S.w <- rowSums(fs %>% select(matches('^S(\\d|C)\\.w$')))
  
  fileSfx <- paste0('Signals_ica')

  ggplot(bind_cols(b,FS=env$FS$Si) %>% rowid_to_column('t') %>%
           gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Base Timeseries (Bi) and Target Sum (FS)')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_BaseSignals',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  fs <- fs %>% mutate(Si = Si-mean(Si)) %>% mutate(Si = env$baseX*Si/max(abs(Si))) %>%  
    mutate(S.w = S.w-mean(S.w)) %>% mutate(S.w = env$baseX*S.w/max(abs(S.w)))
  x.x <- env$X %>% select(env$xColNames)
  FS.scaled <- env$FS$Si-mean(env$FS$Si) 
  FS.scaled <- env$baseX * FS.scaled /max(abs(FS.scaled)) 
  
  ggplot(bind_cols(x.x, FS = FS.scaled) %>% rowid_to_column('t') %>%
           gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal))+
    #ggtitle(paste0(titleSfx,'Target Signal (FS), ICA Signals')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_ICASignals',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  x.xWorst <- env$ica.worst
  ggplot(bind_cols(x.xWorst, FS = FS.scaled) %>% rowid_to_column('t') %>%
           gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal))+
    #ggtitle(paste0(titleSfx,'Target Signal (FS), Worst ICA Signals')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_WorstICASignals',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(fs %>% gather(Signal,Amplitude,matches('^Si$|^S.?\\.w$'),-t))+
    geom_line(aes(t,Amplitude,color=Signal))+
    #ggtitle(paste0(titleSfx,'Target Signal (FS), Weighted ICA Signals (W*Si)')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,S.w), show.legend = F, linetype = 'dotted') +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_WeightedICASignals',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')


  ggplot(as_tibble(env$X %>% select(env$xColNames)) %>% rowid_to_column('t') %>%
           gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Raw ICA Signals (Si)')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,S.w), show.legend = F, color='red') +
    geom_line(data=fs, mapping=aes(t,S.w), show.legend = F, linetype = 'dotted') +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_RawSignals',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(as_tibble(env$X.orig %>% select(-matches('S10'))) %>% rowid_to_column('t') %>% gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Artificial TimeSeries generated from BaseSignals Bi')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_ArtificialSeries',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  #ggplot(b %>% rowid_to_column('t') %>% gather(Signal,Amplitude,-t))+geom_line(aes(t,Amplitude,color=Signal))
  #ggplot(select(env$FS,Si) %>% rowid_to_column('t') %>% gather(Signal,Amplitude,-t))+geom_line(aes(t,Amplitude,color=Signal))
  #ggplot(select(env$FS,S) %>% rowid_to_column('t') %>% gather(Signal,Amplitude,-t))+geom_line(aes(t,Amplitude,color=Signal))
  #ggplot(select(env$FS,p) %>% rowid_to_column('t') %>% gather(Signal,Amplitude,-t))+geom_line(aes(t,Amplitude,color=Signal))
  #ggplot(env$A %>% filter(n>0)) + geom_density_2d(aes(sid,aid))
  #ggplot(env$A %>% filter(n>0) %>% group_by(sid,aid) %>% summarise(n = sum(n))) + geom_point(aes(sid,aid,size=n))


  fileSfx <- paste0('Tracks_ica')
  track <- track %>% inner_join(fs, suffix = c('','.FS'), by='t')
  track$F.Si <- env$FS$Si[1:nrow(track)]

  track <- track %>% mutate_at(vars(matches(regex('^res\\..*R\\d?$', ignore_case = F))),
                      list(~./env$xColNum))

  ggplot(track %>% select(t,matches('^S.i$')) %>% gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Resulting State-Track following generated Policy'), ) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=track, mapping = aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_StateTrack_Si',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(track %>% select(t,matches('^res\\.R\\d$')) %>% gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Resulting Reward vector following generated Policy')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Reward_Vector',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(track %>% select(t,matches('^res\\.SumR\\d$')) %>% gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Resulting Sum(Reward) vector following generated Policy')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Sum_Reward_Vector',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(track %>% select(t,matches('^res\\..*R$')) %>% gather(Signal,Amplitude,-t))+
    geom_line(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Resulting Reward following generated Policy')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Reward',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  #ggplot(track %>% select(t,matches('S\\di1')) %>% gather(Signal,Amplitude,-t))+geom_line(aes(t,Amplitude,color=Signal))

  #ggplot(env$S %>% select(matches('S\\d'),P,n) %>% gather(Signal,Amplitude,-P,-n)) +
  #ggplot(env$S %>% select(matches('S\\d'),n,P)) + geom_contour(aes_string(x='S1',y='S2', z='n'))+facet_wrap(~P)
}

run_many <- function(env) {
  for (i in 1:9) {
    assign(paste0('res',i), getTrack(env))
  }
  res <- gather(res1$track %>% select(matches('S\\di')), Signal, Amplitude)
  for (i in 2:9) {
    res <- bind_cols(res, gather(get(paste0('res',i))$track %>% select(matches('S\\di')), Signal, Amplitude) %>% select(Amplitude))
  }
  ggplot(res %>% rowid_to_column('t') %>% gather(k,v,-t,-Signal)) + geom_line(aes(t,v,color=k))
}

evalLog <- function(res, env = res$myEnv, tag = '') {
  print(paste0('processing result: ',env$numSigs,'/',env$numOrig,': logs'))

  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,': ')
  fileSfx <- paste0('eval',tag)
  fileExt <- paste0('_',env$numSigs,'_',env$numOrig,'.png')
  tailNum <- tailRuns*(env$tSteps-1)

  idxs <- str_extract(env$xColNames, '.$')
  fs <- env$FS %>% rowid_to_column('t')
  l <- env$log %>% inner_join(fs, suffix = c('','.FS'), by='t')
  l$Qpo <- rowSums(l %>% select(matches('Qpo.S\\d')))
  for (s in env$xColNames) {
    l <- l %>% mutate(!!paste0('dqw.',s) := !!as.name(paste0('dQ.',s)) * !!as.name(paste0('W.',s)))
  }
  ggplot(l %>% gather(Weight,Value,matches('^[W]\\.*'),-t))+geom_smooth(aes(i,Value,color=Weight))
  ggsave(paste0(fileSfx,'_i_W',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  ggplot(l %>% gather(Weight,Value,matches('^[W]\\.*'),-t))+geom_smooth(aes(t,Value,color=Weight))
  ggsave(paste0(fileSfx,'_W',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  ggplot(tail(l,tailRuns) %>% gather(Weight,Value,matches('^[W]\\.*'),-t))+geom_line(aes(t,Value,color=Weight))
  ggsave(paste0(fileSfx,'_W_500',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  l.w <- l %>% select(matches('^SVal\\.S.$|^W\\.S.$|^[ti]$')) %>% mutate(S.w = 0)
  for (idx in idxs) l.w <- l.w %>%
    mutate(!!paste0('SVal.S',idx,'.w') := !!as.name(paste0('W.S',idx)) * !!as.name(paste0('SVal.S',idx)))
  l.w$SVal.S.w <- rowSums(select(l.w, matches('SVal\\.S.\\.w')))
  l.w <- l.w %>% mutate(SVal.S.w - min(SVal.S.w)) %>% mutate(SVal.S.w = env$baseX*SVal.S.w/max(abs(SVal.S.w)))

  ggplot(l.w %>% select(matches('^SVal\\.S.$|^t$')) %>% gather(Signal,Amplitude,-t))+
    geom_smooth(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'Steps SVals')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Steps',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(l.w %>% select(matches('^SVal\\.S.?.w$|^t$')) %>% gather(Signal,Amplitude,-t))+
    geom_smooth(aes(t,Amplitude,color=Signal)) +
    #ggtitle(paste0(titleSfx,'W Weighted Steps (SVals.i * W)')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_WeightedSteps',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')


  #### Evaluate Variance W
  log <- env$log
  w.max <- tibble(.rows=1)
  W.diff <- tibble(.rows = length(log$W.S1)-1)
  for (n in grep('W\\.S.',names(log), value = T)) {
    W.diff[n] <- diff(pull(log,n))
    print(paste('Diff: ', n, 'val:',max(abs(pull(W.diff,n))),'idx:',which.max(abs(pull(W.diff,n)))))
    w.max[n] <- which.max(abs(pull(W.diff,n)))
  }

  wmx <- unlist(w.max)
  wmx <- c(wmx-1, wmx, wmx+1) %>% sort()
  log[unlist(wmx),] %>% select(t,time,i,t,sid,epsilon,matches('dQ\\.|W'))

  ggplot(W.diff %>% tail(5000) %>% gather(Signal,Diff)) + geom_histogram(aes(Diff))+
    scale_y_log10() + facet_wrap(~Signal)+
    #ggtitle(paste0(titleSfx,'W Histogram')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_W_Histogram',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')



  #ggplot(l %>% gather(Signal,Amplitude,matches('^[Z]\\.*'),-t,-n))+geom_smooth(aes(n,Amplitude,color=Signal))+
  #  facet_wrap(~t)+theme(strip.text.x = element_text(size = 5))
  #ggsave(paste0(fileSfx,'_Z_t',fileExt))
  #ggplot(l %>% gather(Signal,Amplitude,matches('d[Q]\\..*'),-sid))+geom_smooth(aes(sid,Amplitude,color=Signal))
  #ggsave(paste0(fileSfx,'_dQ',fileExt))
  #ggplot(l %>% gather(Signal,Amplitude,matches('d[Q]\\..*'),-n,-t))+geom_smooth(aes(n,Amplitude,color=Signal))+
  #  facet_wrap(~t)+theme(strip.text.x = element_text(size = 5))
  #ggsave(paste0(fileSfx,'_dQ_t',fileExt))
  #ggplot(l %>% gather(Signal,Amplitude,matches('dqw\\..*'),-t))+geom_smooth(aes(t,Amplitude,color=Signal))
  #ggsave(paste0(fileSfx,'_dqw',fileExt))
  #ggplot(l %>% gather(Signal,Amplitude,matches('dqw\\..*'),-n,-t))+geom_smooth(aes(n,Amplitude,color=Signal))+
  #  facet_wrap(~t)+theme(strip.text.x = element_text(size = 5))
  #ggsave(paste0(fileSfx,'_dqw_t',fileExt))

  #comb <- gtools::combinations(env$xColNum-1,2)
  #for (i in 1:nrow(comb)) {
  #  si <- paste0('S',comb[i,])
  #  #ggplot(l %>% select(si,matches('^(sid|Qpo|sig.p)$')) %>%
  #  #         group_by_at(c(si,'sig.p')) %>% summarise(Qpo = mean(Qpo)))+
  #  #  geom_contour(aes_string(x=si[1],y=si[2],z='Qpo'))+
  #  #  facet_wrap(~sig.p)+theme(strip.text.x = element_text(size = 5))+
  #  ggplot(l %>% select(si,matches('^(sid|Qpo)$')) %>%
  #           group_by_at(si) %>% summarise(Qpo = mean(Qpo)))+
  #    geom_contour(aes_string(x=si[1],y=si[2],z='Qpo'))+
  #    theme(strip.text.x = element_text(size = 5))+
  #    ggtitle(paste0(titleSfx,'State-Evaluation Sum(Qpo.i) - States ',si[2],'~',si[1]))
  #  ggsave(paste0(fileSfx,'_SumQ_',i,fileExt))
  #}
}

evalTracks <- function(res, env = res$myEnv, tracks = res$tracks) {
  print(paste0('processing result: ',env$numSigs,'/',env$numOrig,': tracks'))

  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,': ')
  fileSfx <- paste0('track')
  fileExt <- paste0('_',env$numSigs,'_',env$numOrig,'.png')
  tailNum <- tailRuns*(env$tSteps-1)

  idxs <- str_extract(env$xColNames, '.$')
  fs <- env$FS %>% rowid_to_column('t') %>% rename(Si.FS=Si)
  tracks.all <- tracks %>%
    mutate_at(vars(matches(regex('^res\\..*R\\d?$', ignore_case = F))), list(~./env$xColNum)) %>%
    mutate_at(vars(i), as.factor) %>%
    mutate(t = stid) %>%
    inner_join(fs, suffix = c('','.FS'), by='t') %>%
    left_join(env$X %>% set_names(paste0('X.',names(.))), by=c('t'='X.xid'))
  tracks <- tracks.all %>% tail(tailNum)

  trs.all <- tracks.all %>%
    select(matches('^SVal\\.S.$|^S.?i?$|^W\\.S.$|^S.i.w$|^t$|^i$|^X\\.S.$')) %>%
    rowid_to_column('id') %>%
    inner_join(fs, suffix = c('','.FS'), by='t') %>%
    mutate(Si.w = 0, i = as.factor(i))
  trs <- trs.all %>% tail(tailNum)

  for (idx in idxs) trs <-
    mutate(trs, !!as.name(paste0('S',idx,'.w')) := !!as.name(paste0('W.S',idx))*!!as.name(paste0('X.S',idx)))
  for (idx in idxs) trs <-
    mutate(trs, Si.w := Si.w + !!as.name(paste0('W.S',idx))*!!as.name(paste0('X.S',idx)))
  trs <- trs %>% mutate(Si.w = Si.w-min(Si.w)) %>% mutate(Si.w = env$baseX * Si.w / max(abs(Si.w)))

  for (idx in idxs) tracks <-
    mutate(tracks, !!as.name(paste0('S',idx,'i.w')) := !!as.name(paste0('W.S',idx))*!!as.name(paste0('X.S',idx)))
  tracks$Si.w <- rowSums(tracks %>% select(matches('S.i\\.w')))
  tracks <- tracks %>% group_by(i) %>% mutate(Si.w = Si.w - min(Si.w)) %>% ungroup()
  tracks <- tracks %>% mutate(Si.w = env$baseX * Si.w / max(abs(Si.w)))

  tsi <- tracks %>% select(t,i,matches('^SVal\\.S.$')) %>% gather(Signal,Amplitude,-t,-i)
  sumR <- tracks %>% filter(t==19) %>%
    rename(SumR = res.SumR) %>%
    select(i,SumR)
  topSumR <- sumR %>% arrange(desc(SumR)) %>% head(5)

  tsi.all <- tracks.all %>% select(t,i,matches('^SVal\\.S.$')) %>% gather(Signal,Amplitude,-t,-i)
  tsi <- tsi.all %>% left_join(sumR, by='i')
  topTsi <- tsi %>% filter(i %in% topSumR$i)
  topTracks <- tracks %>% filter(i %in% topSumR$i)

  tsi.w <- tracks %>% select(t,i,matches('^S.?i.w$')) %>% gather(Signal,Amplitude,-t,-i)
  tsi.w <- tsi.w %>% left_join(sumR, by='i')
  topTsi.w <- tsi.w %>% filter(i %in% topSumR$i)

  ggplot()+
    geom_point(data=trs %>% tail(19*10), mapping=aes(t,Si.w, color='Tracks'), linetype='dotted', color='black', size=1) +
    geom_smooth(data=trs, mapping=aes(t,Si.w, color='FinalTracks'), linetype='dotted', color='red', size=2, span=.2) +
    geom_smooth(data=tracks, mapping=aes(t,Si.w, color='AllTracks'), linetype='dashed', color='green', size=.5, span=.2)

  ggplot()+
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')

  ggplot(trs)+geom_line(aes(t,S,color=i))+
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    #geom_line(aes(t,Si.w), linetype='dotted', color='red', size=2)+
    geom_line(aes(t,Si.FS), linetype='dashed')+
    theme(legend.position = 'none')
  ggsave(paste0(fileSfx,'_WeightedTrack_FS',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(tail(tracks,50*(env$tSteps-1)) %>% gather(Signal,Amplitude,matches('^(S.i\\.w)$'),-t,-i))+
    geom_smooth(aes(t,Amplitude,color=Signal), span=.2) #+
    #geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2)+
    #geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_WeightedTracks_Tail',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(topTracks %>% gather(Signal,Amplitude,matches('^(S.?i\\.w)|Si\\.FS$'),-t,-i))+
    geom_smooth(aes(t,Amplitude,color=Signal), span=.2) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_WeightedTracks_Top',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')


  ggplot(tsi)+geom_smooth(aes(t,Amplitude,color=Signal)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_Tracks_smooth',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  ggplot(tsi %>% mutate(ti=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,Amplitude,color=ti), position = position_jitter(.3,.3))+
    theme(legend.position = 'none') +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
    #+ ggtitle(paste0(titleSfx,'Tracks'))
  ggsave(paste0(fileSfx,'_Tracks',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  ggplot(topTsi %>% mutate(Signal=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,Amplitude,color=Signal), position = position_jitter(.3,.3)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
    #+ ggtitle(paste0(titleSfx,'most successful Tracks'))
  ggsave(paste0(fileSfx,'_Tracks_Top5',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(tsi.w)+geom_smooth(aes(t,Amplitude,color=Signal)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
    #+ ggtitle(paste0(titleSfx,'Weighted Tracks'))
  ggsave(paste0(fileSfx,'_WeightedTracks_smooth',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  ggplot(tsi.w %>% mutate(Signal=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,Amplitude,color=Signal), position = position_jitter(.3,.3))+
    theme(legend.position = 'none') +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
    #+ ggtitle(paste0(titleSfx,'Weighted Tracks'))
  ggsave(paste0(fileSfx,'_WeightedTracks',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  ggplot(topTsi.w %>% mutate(Signal=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,Amplitude,color=Signal), position = position_jitter(.3,.3)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
    #+ ggtitle(paste0(titleSfx,'most successful Weighted Tracks'))
  ggsave(paste0(fileSfx,'_WeightedTracks_Top5',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  trs <- tracks
  tr.stat <- tracks %>% filter(t > 5)
  r.mean <- mean(tracks$act.R)
  r.center <- tr.stat %>% group_by(t) %>% summarize(R.center = mean(act.R))
  r.sd <- inner_join(tr.stat, r.center, by='t', suffix=c('','.c')) %>%
    group_by(t) %>% summarize(R.sd = sd(act.R-R.center))
  r.csd <- inner_join(r.center, r.sd, by='t', suffix=c('','.sd'))
  r.centerMean <- mean(r.center$R.center %>% head(-5))
  r.y <- max(r.mean,r.centerMean)

  ggplot(trs)+geom_smooth(aes(t,act.R), span=.2) +
    geom_hline(yintercept = r.mean, linetype='dashed') +
    geom_text(aes(0,r.y), label=paste('mean(R):',format(r.mean,digits=3)), vjust=-2.5) +
    geom_hline(yintercept = r.centerMean, linetype='dashed', color='red') +
    geom_text(aes(0,r.y), label=paste('RCenter:',format(r.centerMean,digits=3)), colour='red', vjust=-1) +
    geom_line(data=r.csd, mapping=aes(t,R.center+R.sd), linetype='dotted', size=.5) +
    geom_line(data=r.csd, mapping=aes(t,R.center-R.sd), linetype='dotted', size=.5)
    #+ ggtitle(paste0(titleSfx,'Reward'))
  ggsave(paste0(fileSfx,'_Reward_smooth',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(trs %>% filter(stid==19) %>% mutate(i = as.numeric(as.character(i))))+
    geom_line(aes(i,res.SumR)) +
    geom_smooth(aes(i,res.SumR), span=.5)
    #+ ggtitle(paste0(titleSfx,'final SumReward'))
  ggsave(paste0(fileSfx,'_SumReward_Tail',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

  ggplot(trs %>% filter(stid==19) %>% mutate(i = as.numeric(as.character(i))) %>%
           gather(w,Value,matches('W.S.'),-i))+
    geom_line(aes(i,Value,color=w)) +
    geom_smooth(aes(i,Value,color=w), linetype='dashed', span=.3, size=.5)
    #+ ggtitle(paste0(titleSfx,'final W'))
  ggsave(paste0(fileSfx,'_W_Tail',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')

}

runOptim <- function(myXtsFs) {
  #ofs <- xts.fs[5:nrow(xts.fs),]
  ofs <- myXtsFs[5:(nrow(myXtsFs)-5),c('Si','S.w')]
  ofs <- scale(ofs, center=T, scale=T)

  doLag <- function(d.t) {
    dt<-round(d.t);
    #print(paste('d.t=',d.t,'dt=',dt))
    mean((ofs$Si[1:(nrow(ofs)-dt)]-ofs$S.w[(dt+1):nrow(ofs)])^2, na.rm=T)
  }
  #o <- optim(c(0), doLag, method = 'L-BFGS-B', lower=0, upper=5, control=list(trace=3, num))
  #o <- optimx::optimx(c(.3), doLag, lower=0, upper=1, control=list(trace=3, maxit=10000, ndeps=seq(0,1,.01)))

  o <- map_dbl(seq(0,nrow(myXtsFs)/10), doLag)
  which.min(o)
}

printTracks <- function(res, env = res$myEnv, tracks = res$tracks) {
  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,': ')
  fileSfx <- paste0('tracks_',env$numSigs,'_',env$numOrig)
  tailNum <- tailRuns*(env$tSteps-1)

  tracks.all <- tracks %>% mutate(t = stid) %>%
    mutate_at(vars(matches(regex('^res\\..*R\\d?$', ignore_case = F))),
                               list(~./env$xColNum))%>%
    left_join(env$X %>% set_names(paste0('X.',names(.))), by=c('t'='X.xid'))
  tracks <- tracks.all %>% tail(tailNum)

  tsi.all <- tracks.all %>% select(t,i,matches('^SVal\\.S.$')) %>% gather(si,v,-t,-i)
  tsi <- tsi.all %>% tail(tailNum)

  sumR <- tracks %>% filter(t==19) %>%
    rename(SumR = res.SumR) %>%
    select(i,SumR)
  w.s <- colMeans(env$log %>% tail(tailRuns) %>% select(matches('W\\.S.')))
  W <- env$W
  fs <- env$FS %>% rowid_to_column('t')
  for (n in env$xColNames) fs[paste0(n,'.w')] <- w.s[paste0('W.',n)]*env$X[n]
  fs$S.w <- rowSums(fs %>% select(matches('^S(\\d|C)\\.w$')))

  tr.stat <- tail(tracks, tailNum) %>% filter(t > 4)
  r.mean <- mean(tail(tracks,tailNum)$act.R)
  r.center <- tr.stat %>% group_by(t) %>% summarize(R.center = mean(act.R))
  r.sd <- inner_join(tr.stat, r.center, by='t', suffix=c('','.c')) %>%
    group_by(t) %>% summarize(R.sd = sd(act.R-R.center))
  r.csd <- inner_join(r.center, r.sd, by='t', suffix=c('','.sd'))

  numVar <- ncol(fs)
  d <- as_datetime(1000*fs$t)
  xts.fs <- xts::xts(as.matrix(fs), order.by = d)
  xts.fs <- xts::merge.xts(xts.fs, as_datetime(1000*seq(1,max(fs$t),1/10)), fill=NA)
  xts.fs <- zoo::na.spline(xts.fs[,1:numVar])
  xts.fs <- scale(xts.fs)
  png(paste0(fileSfx,'_timeline.png'))
  print(plot(xts.fs[,c('S.w','Si')], main = paste(titleSfx, 'orig.')))
  dev.off()

  xts.fs.diff <- scale(diff(xts.fs))
  png(paste0(fileSfx,'_timeline_diff.png'))
  print(plot(xts.fs.diff[,c('S.w','Si')], main = paste(titleSfx, 'Lag:', format(lag, digits=2), 'differentiated')))
  dev.off()

  lag <- runOptim(xts.fs)
  xts.fs.lag <- xts.fs
  xts.fs.lag$Si <- xts::lag.xts(xts.fs.lag$Si, lag)
  xts.fs.lag <- scale(xts.fs.lag)
  png(paste0(fileSfx,'_timeline_lag.png'))
  print(plot(xts.fs.lag[,c('S.w','Si')], main = paste(titleSfx, 'Lag:',format(lag, digits=2))))
  dev.off()

  xts.fs.lag.diff <- xts.fs.lag
  xts.fs.lag.diff$Si <- xts::lag.xts(xts.fs$Si, lag)
  xts.fs.lag.diff <- diff(xts.fs.lag.diff)
  xts.fs.lag.diff <- scale(xts.fs.lag.diff)
  png(paste0(fileSfx,'_timeline_lag_diff.png'))
  print(plot(xts.fs.lag.diff[,c('S.w','Si')], main = paste(titleSfx, 'Lag:', format(lag, digits=2), 'differentiated')))
  dev.off()

  print(paste(titleSfx,'Track: ',titleSfx,', W=',toString(w.s)))
  print(w.s)
  print(paste(titleSfx,'Summary of SumR:'))
  print(summary(sumR$SumR))
  print(paste(titleSfx,'Track: Reward=',toString(r.center$R.center)))
  print(paste(titleSfx,'Track: Reward: Sum=',toString(sum(r.center$R.center)),
              ', Mean=',toString(mean(r.center$R.center)),
              ', SD=',toString(sd(r.sd$R.sd))))
  print(paste(titleSfx,'Track: Reward Summary=',toString(summary(r.center$R.center))))
  print(paste(titleSfx,'Track: Reward Summary:'))
  print(summary(r.center$R.center))

  print(paste(titleSfx,'Summary of FS - Sum(W*Xi)=',toString(summary(fs$Si - fs$S.w))))
  print(paste(titleSfx,'Summary of FS - Sum(W*Xi):'))
  print(summary(fs$Si - fs$S.w))

  snames <- na.exclude(names(xts.fs) %>% str_extract('^Si$|^S\\d?\\.w$'))

  print(paste(titleSfx,'Correlation of Signals:'))
  print(suppressWarnings(corrr::correlate(fs[,snames])))

  print(paste(titleSfx,'Correlation of Diff-Signals:'))
  print(suppressWarnings(corrr::correlate(xts.fs.diff[,snames], quiet = T)))

  print(paste(titleSfx,'Correlation of Signals - Lag =',lag,':'))
  print(suppressWarnings(corrr::correlate(xts.fs.lag[,snames], quiet = T)))

  print(paste(titleSfx,'Correlation of Diff-Signals - Lag =',lag,':'))
  print(suppressWarnings(corrr::correlate(xts.fs.lag.diff[,snames], quiet = T)))

  #print(sumR %>% arrange(desc(SumR)), n=20)
}

corrTrackPerRun <- function(res, env = res$myEnv, run = 1) {
  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,'[',run,']',': ')
  fileSfx <- paste0('tracks_',env$numSigs,'_',env$numOrig,'_',run)
  tailNum <- tailRuns*(env$tSteps-1)

  tracks <- env$tracks
  tracks.all <- tracks %>% mutate(t = stid) %>%
    mutate_at(vars(matches(regex('^res\\..*R\\d?$', ignore_case = F))),
              list(~./env$xColNum))%>%
    left_join(env$X %>% set_names(paste0('X.',names(.))), by=c('t'='X.xid'))
  tracks <- tracks.all %>% tail(tailNum)

  tsi.all <- tracks.all %>% select(t,i,matches('^SVal\\.S.$')) %>% gather(si,v,-t,-i)
  tsi <- tsi.all %>% tail(tailNum)

  sumR <- tracks %>% filter(t==19) %>%
    rename(SumR = res.SumR) %>%
    select(i,SumR)
  w.s <- colMeans(env$log %>% tail(tailRuns) %>% select(matches('W\\.S.')))
  W <- env$W
  fs <- env$FS %>% rowid_to_column('t')
  for (n in env$xColNames) fs[paste0(n,'.w')] <- w.s[paste0('W.',n)]*env$X[n]
  fs$S.w <- rowSums(fs %>% select(matches('^S(\\d|C)\\.w$')))

  snames <- na.exclude(names(fs) %>% str_extract('^Si$|^S\\d?\\.w$'))

  print(paste(titleSfx,'Correlation of Signals:'))
  print(suppressWarnings(corrr::correlate(fs[,snames])))

  fscorr <- corrr::correlate(fs[,snames])
  return(fscorr)
}


icns <- c(2:10)
#icns <- c(3)
#baseDir <- list('v3sig3'='storage','v3sig4'='storage','v3sig2'='tfrlnx.fruehbeck.at')
#baseDir <- list('v4sig2'='tfrlnx.fruehbeck.at','v4sig3'='storage','v4sig4'='storage')
#baseDir <- list('v5sig2'='tfrlnx.fruehbeck.at','v5sig3'='storage','v5sig4'='storage')
#baseDir <- list('v6sig2'='tfrlnx.fruehbeck.at','v6sig3'='storage','v6sig4'='storage')
#baseDir <- list('v6sig2'='storage','v6sig3'='storage')
#baseDir <- list('v6sig3'='monster')
baseDir <- list('v6sig2'='monster','v6sig3'='monster','v6sig4'='monster')

if (!is.null(sig.v)) {
  #bd <- str_detect(names(baseDir),paste0('v\\dsig',sig.v))
  baseDir <- sig.v
}
if (!is.null(orig.s)) icns <- c(as.numeric(orig.s))

dn <- NA
icn <- NA

if (evalPlot) {
  print(paste('Running evalPlot args=',toString(args)))
  for (dn in baseDir) {
    currDir <- getwd()
    if (local) dn <- '.'
    print(paste('switched to:',getwd()))
    tryCatch({
      host <- baseDir[dn]

      for (icn in icns) {
        file <- paste0('finDemo.env_',host,'_ica_',icn,'.rdata')
        if (!file.exists(file)) host <- node
        file <- paste0('finDemo.env_',host,'_ica_',icn,'.rdata')

        if (!file.exists(file)) {
          print(paste('file not existing:',file))
          next()
        }

        print(paste('loading:',file))
        try({
          loaded <- F
          for (loadTries in 1:40) {
            tryCatch({ load(file); loaded <- T; },
                     error=function(e) { print(paste('cannot load',file,e)); Sys.sleep(5) })
            if (loaded) {
              print(paste('loaded:',file))
              break;
            }
          }
          if (!loaded) stop(paste('failed to load file:',file))
          if (test) next

          assign(paste0('x',icn),myEnv)

          #myEnv <- get(paste0('x',icn))
          result <- getTrack(myEnv)
          result$tracks <- tracks
          result$myEnvs <- myEnvs
          result$runCorrs <- list()

          for (run in 1:length(result$myEnvs)) {
            runCorr <- corrTrackPerRun(result, env = result$myEnvs[[run]], run = run)
            runCorr$run <- run
            result$runCorrs[[run]] <- runCorr
            runTib <- runCorr %>% filter(rowname == 'S.w')
            if (run == 1) {
              result$runCorr <- runTib
            } else {
              result$runCorr <- bind_rows(result$runCorr, runTib)
            }
          }

          if (icn < 5) assign(paste0(dn,'r',icn),result)
          save(result, file=paste0('finDemo.env_',host,'_ica_',icn,'_Res.rdata'))

          plotResults(result)
          evalLog(result)
          for (i in 1:length(result$myEnvs)) evalLog(result, result$myEnvs[[i]], paste0('_I',i,'_'))
          evalTracks(result)
        })
      }
    },
      error = function(err) { print(paste('Error:',dn,'/',icn,':',toString(err))) },
      finally = setwd(currDir))
  }
}

if (evalPrint) {
  print(paste('Running evalPrint args=',toString(args)))
  for (dn in baseDir) {
    host <- baseDir[dn]
    for (icn in icns) {
      try({
        if (local) dn <- '.'
        file <- paste0(dn,'/','finDemo.env_',host,'_ica_',icn,'_Res.rdata')
        if (!file.exists(file)) {
          host <- node
          file <- paste0(dn,'/','finDemo.env_',host,'_ica_',icn,'_Res.rdata')
        }

        print(paste('loading:',file))
        load(file)
        print(paste('loaded:',file))

        printTracks(result)

      })
    }
  }
}

