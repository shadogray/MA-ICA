library(tidyverse)
#library(timeSeries)
#library(corrr)
library(lubridate)

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
  result <- findWay(env)
  b <- tibble(.rows = length(env$B1))
  for (n in names(env)[grep('B\\d',names(env))]) b <- bind_cols(b, env[n])
  if (!has_name(result$track,'t')) result$track <- result$track %>% rowid_to_column('t')
  result
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
  fileSfx <- paste0('S_density_ica_',env$numSigs,'_',env$numOrig)

  b <- as_tibble(as.matrix(res$b)*env$baseX)
  w.s <- colMeans(env$log %>% tail(200) %>% select(matches('W\\.S.')))
  W <- env$W / sum(abs(env$W))
  fs <- env$FS %>% rowid_to_column('t')
  for (n in env$xColNames) fs[paste0(n,'.w')] <- w.s[paste0('W.',n)]*env$X[n]
  fs$S.w <- rowSums(fs %>% select(matches('^S(\\d|C)\\.w$')))

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
    ggplot(env$log %>% select(si) %>% mutate_all(list(~round(.))), mapping=aes_string(si[1],si[2])) +
      geom_density_2d() +
      theme(strip.text.x = element_text(size = 5)) +
      ggtitle(paste0(titleSfx,'State-Evaluation Density - States ',si[2],'~',si[1]))
    ggsave(paste0(fileSfx,'_',i,'.png'))
  }

  fileSfx <- paste0('Signals_ica_',env$numSigs,'_',env$numOrig)

  ggplot(bind_cols(b,FS=env$FS$Si) %>% rowid_to_column('tx') %>%
           gather(signal,value,-tx))+
    geom_line(aes(tx,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Base Timeseries (Bi) and Target Sum (FS)')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_BaseSignals','.png'))

  x.x <- env$X %>% select(env$xColNames)
  ggplot(bind_cols(x.x, FS = env$FS$Si) %>% rowid_to_column('tx') %>%
           gather(signal,value,-tx))+
    geom_line(aes(tx,value,color=signal))+
    ggtitle(paste0(titleSfx,'Target Signal (FS), ICA Signals')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_ICASignals','.png'))

  ggplot(fs %>% gather(signal,value,matches('^Si$|^S.?\\.w$'),-t))+
    geom_line(aes(t,value,color=signal))+
    ggtitle(paste0(titleSfx,'Target Signal (FS), Weighted ICA Signals (W*Si)')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_WeightedICASignals','.png'))


  ggplot(as_tibble(env$X %>% select(env$xColNames)) %>% rowid_to_column('tx') %>%
           gather(signal,value,-tx))+
    geom_line(aes(tx,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Raw ICA Signals (Si)')) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=fs, mapping=aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_RawSignals','.png'))

  ggplot(as_tibble(env$X.orig) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+
    geom_line(aes(tx,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Artificial TimeSeries generated from BaseSignals Bi')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_ArtificialSeries','.png'))

  #ggplot(b %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  #ggplot(select(env$FS,Si) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  #ggplot(select(env$FS,S) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  #ggplot(select(env$FS,p) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  #ggplot(env$A %>% filter(n>0)) + geom_density_2d(aes(sid,aid))
  #ggplot(env$A %>% filter(n>0) %>% group_by(sid,aid) %>% summarise(n = sum(n))) + geom_point(aes(sid,aid,size=n))


  fileSfx <- paste0('Tracks_ica_',env$numSigs,'_',env$numOrig)
  track <- track %>% inner_join(fs, suffix = c('','.FS'), by='t')
  track$F.Si <- env$FS$Si[1:nrow(track)]

  track <- track %>% mutate_at(vars(matches(regex('^res\\..*R\\d?$', ignore_case = F))),
                      list(~./env$xColNum))

  ggplot(track %>% select(t,matches('^S.i$')) %>% gather(signal,value,-t))+
    geom_line(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Resulting State-Track following generated Policy'), ) +
    theme(plot.title = element_text(size = 10)) +
    geom_line(data=track, mapping = aes(t,Si), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_StateTrack_Si','.png'))

  ggplot(track %>% select(t,matches('^res\\.R\\d$')) %>% gather(signal,value,-t))+
    geom_line(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Resulting Reward vector following generated Policy')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Reward_Vector','.png'))

  ggplot(track %>% select(t,matches('^res\\.SumR\\d$')) %>% gather(signal,value,-t))+
    geom_line(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Resulting Sum(Reward) vector following generated Policy')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Sum_Reward_Vector','.png'))

  ggplot(track %>% select(t,matches('^res\\..*R$')) %>% gather(signal,value,-t))+
    geom_line(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Resulting Reward following generated Policy')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Reward','.png'))

  #ggplot(track %>% select(t,matches('S\\di1')) %>% gather(signal,value,-t))+geom_line(aes(t,value,color=signal))

  #ggplot(env$S %>% select(matches('S\\d'),P,n) %>% gather(key,value,-P,-n)) +
  #ggplot(env$S %>% select(matches('S\\d'),n,P)) + geom_contour(aes_string(x='S1',y='S2', z='n'))+facet_wrap(~P)
}

run_many <- function(env) {
  for (i in 1:9) {
    assign(paste0('res',i), getTrack(env))
  }
  res <- gather(res1$track %>% select(matches('S\\di')), key, value)
  for (i in 2:9) {
    res <- bind_cols(res, gather(get(paste0('res',i))$track %>% select(matches('S\\di')), key, value) %>% select(value))
  }
  ggplot(res %>% rowid_to_column('t') %>% gather(k,v,-t,-key)) + geom_line(aes(t,v,color=k))
}

evalLog <- function(res, env = res$myEnv) {
  print(paste0('processing result: ',env$numSigs,'/',env$numOrig,': logs'))

  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,': ')
  fileSfx <- paste0('eval_',env$numSigs,'_',env$numOrig)

  idxs <- str_extract(env$xColNames, '.$')
  fs <- env$FS %>% rowid_to_column('t')
  l <- env$log %>% inner_join(fs, suffix = c('','.FS'), by='t')
  l$Qpo <- rowSums(l %>% select(matches('Qpo.S\\d')))
  for (s in env$xColNames) {
    l <- l %>% mutate(!!paste0('dqw.',s) := !!as.name(paste0('dQ.',s)) * !!as.name(paste0('W.',s)))
  }
  ggplot(l %>% gather(key,val,matches('^[W]\\.*'),-lid))+geom_smooth(aes(lid,val,color=key))
  ggsave(paste0(fileSfx,'_W','.png'))
  ggplot(tail(l,500) %>% gather(key,val,matches('^[W]\\.*'),-lid))+geom_line(aes(lid,val,color=key))
  ggsave(paste0(fileSfx,'_W_500','.png'))

  l.w <- l %>% select(matches('^SVal\\.S.$|^W\\.S.$|^[ti]$')) %>% mutate(S.w = 0)
  for (idx in idxs) l.w <- l.w %>%
    mutate(!!paste0('SVal.S',idx,'.w') := !!as.name(paste0('W.S',idx)) * !!as.name(paste0('SVal.S',idx)))
  l.w$SVal.S.w <- rowSums(select(l.w, matches('SVal\\.S.\\.w')))

  ggplot(l.w %>% select(matches('^SVal\\.S.$|^t$')) %>% gather(signal,value,-t))+
    geom_smooth(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Steps SVals')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Steps','.png'))

  ggplot(l.w %>% select(matches('^SVal\\.S.?.w$|^t$')) %>% gather(signal,value,-t))+
    geom_smooth(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'W weighted Steps (SVals.i * W)')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_WeightedSteps','.png'))

  #ggplot(l %>% gather(key,val,matches('^[Z]\\.*'),-t,-n))+geom_smooth(aes(n,val,color=key))+
  #  facet_wrap(~t)+theme(strip.text.x = element_text(size = 5))
  #ggsave(paste0(fileSfx,'_Z_t','.png'))
  #ggplot(l %>% gather(key,val,matches('d[Q]\\..*'),-sid))+geom_smooth(aes(sid,val,color=key))
  #ggsave(paste0(fileSfx,'_dQ','.png'))
  #ggplot(l %>% gather(key,val,matches('d[Q]\\..*'),-n,-t))+geom_smooth(aes(n,val,color=key))+
  #  facet_wrap(~t)+theme(strip.text.x = element_text(size = 5))
  #ggsave(paste0(fileSfx,'_dQ_t','.png'))
  #ggplot(l %>% gather(key,val,matches('dqw\\..*'),-lid))+geom_smooth(aes(lid,val,color=key))
  #ggsave(paste0(fileSfx,'_dqw','.png'))
  #ggplot(l %>% gather(key,val,matches('dqw\\..*'),-n,-t))+geom_smooth(aes(n,val,color=key))+
  #  facet_wrap(~t)+theme(strip.text.x = element_text(size = 5))
  #ggsave(paste0(fileSfx,'_dqw_t','.png'))

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
  #  ggsave(paste0(fileSfx,'_SumQ_',i,'.png'))
  #}
}

evalTracks <- function(res, env = res$myEnv, tracks = res$tracks) {
  print(paste0('processing result: ',env$numSigs,'/',env$numOrig,': tracks'))

  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,': ')
  fileSfx <- paste0('track_',env$numSigs,'_',env$numOrig)

  idxs <- str_extract(env$xColNames, '.$')
  fs <- env$FS %>% rowid_to_column('t') %>% rename(Si.FS=Si)
  tracks <- tracks %>%
    mutate_at(vars(matches(regex('^res\\..*R\\d?$', ignore_case = F))), list(~./env$xColNum)) %>%
    mutate_at(vars(i), as.factor) %>%
    mutate(t = stid) %>%
    inner_join(fs, suffix = c('','.FS'), by='t') %>%
    left_join(env$X %>% set_names(paste0('X.',names(.))), by=c('t'='X.xid'))

  trs <- select(tracks,matches('^SVal\\.S.$|^S.?i?$|^W\\.S.$|^S.i.w$|^t$|^i$|^X\\.S.$')) %>%
    rowid_to_column('id') %>%
    inner_join(fs, suffix = c('','.FS'), by='t') %>%
    mutate(Si.w = 0, i = as.factor(i))
  for (idx in idxs) trs <-
    mutate(trs, !!as.name(paste0('S',idx,'.w')) := !!as.name(paste0('W.S',idx))*!!as.name(paste0('X.S',idx)))
  for (idx in idxs) trs <-
    mutate(trs, Si.w := Si.w + !!as.name(paste0('W.S',idx))*!!as.name(paste0('X.S',idx)))

  for (idx in idxs) tracks <-
    mutate(tracks, !!as.name(paste0('S',idx,'i.w')) := !!as.name(paste0('W.S',idx))*!!as.name(paste0('X.S',idx)))
  tracks$Si.w <- rowSums(tracks %>% select(matches('S.i\\.w')))

  tsi <- tracks %>% select(t,i,matches('^SVal\\.S.$')) %>% gather(si,v,-t,-i)
  sumR <- tracks %>% filter(t==19) %>%
    rename(SumR = res.SumR) %>%
    select(i,SumR)
  topSumR <- sumR %>% arrange(desc(SumR)) %>% head(5)
  tsi <- tsi %>% left_join(sumR, by='i')
  topTsi <- tsi %>% filter(i %in% topSumR$i)
  topTracks <- tracks %>% filter(i %in% topSumR$i)

  tsi.w <- tracks %>% select(t,i,matches('^S.?i.w$')) %>% gather(si,v,-t,-i)
  tsi.w <- tsi.w %>% left_join(sumR, by='i')
  topTsi.w <- tsi.w %>% filter(i %in% topSumR$i)

  ggplot()+
    geom_point(data=trs %>% tail(19*10), mapping=aes(t,Si.w), linetype='dotted', color='black', size=1) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2) +
    geom_smooth(data=tracks, mapping=aes(t,Si.w), linetype='dashed', color='green', size=.5, span=.2)

  ggplot()+
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')

  ggplot(tail(trs,50*(env$tSteps-1)))+geom_line(aes(t,S,color=i))+
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    #geom_line(aes(t,Si.w), linetype='dotted', color='red', size=2)+
    geom_line(aes(t,Si.FS), linetype='dashed')+
    theme(legend.position = 'none')
  ggsave(paste0(fileSfx,'_WeightedTrack_FS','.png'))

  ggplot(tail(tracks,50*(env$tSteps-1)) %>% gather(k,s,matches('^(S.i\\.w)$'),-t,-i))+
    geom_smooth(aes(t,s,color=k), span=.2) #+
    #geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2)+
    #geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_WeightedTracks_Tail','.png'))

  ggplot(topTracks %>% gather(k,s,matches('^(S.?i\\.w)|Si\\.FS$'),-t,-i))+
    geom_smooth(aes(t,s,color=k), span=.2) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_WeightedTracks_Top','.png'))


  ggplot(tsi)+geom_smooth(aes(t,v,color=si)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')
  ggsave(paste0(fileSfx,'_Tracks_smooth','.png'))
  ggplot(tsi %>% mutate(ti=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,v,color=ti), position = position_jitter(.3,.3))+
    theme(legend.position = 'none') +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')+
    ggtitle(paste0(titleSfx,'Tracks'))
  ggsave(paste0(fileSfx,'_Tracks','.png'))
  ggplot(topTsi %>% mutate(ti=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,v,color=ti), position = position_jitter(.3,.3)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')+
    ggtitle(paste0(titleSfx,'most successful Tracks'))
  ggsave(paste0(fileSfx,'_Tracks_Top5','.png'))

  ggplot(tsi.w)+geom_smooth(aes(t,v,color=si)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')+
    ggtitle(paste0(titleSfx,'Weighted Tracks'))
  ggsave(paste0(fileSfx,'_WeightedTracks_smooth','.png'))
  ggplot(tsi.w %>% mutate(ti=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,v,color=ti), position = position_jitter(.3,.3))+
    theme(legend.position = 'none') +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')+
    ggtitle(paste0(titleSfx,'Weighted Tracks'))
  ggsave(paste0(fileSfx,'_WeightedTracks','.png'))
  ggplot(topTsi.w %>% mutate(ti=as.factor(paste0(si,'_',i)))) +
    geom_line(aes(t,v,color=ti), position = position_jitter(.3,.3)) +
    geom_smooth(data=trs, mapping=aes(t,Si.w), linetype='dotted', color='red', size=2, span=.2)+
    geom_line(data=tracks, mapping=aes(t,Si.FS), show.legend = F, linetype = 'dashed')+
    ggtitle(paste0(titleSfx,'most successful Weighted Tracks'))
  ggsave(paste0(fileSfx,'_WeightedTracks_Top5','.png'))

  discNum <- round(min(.1*nrow(tracks),200)/(env$tSteps-1))*(env$tSteps-1)
  tailNum <- 500*(env$tSteps-1)
  trs <- tail(tracks, -discNum)

  tr.stat <- tail(tracks, tailNum) %>% filter(t > 4)
  r.mean <- mean(tail(tracks, tailNum)$act.R)
  r.center <- tr.stat %>% group_by(t) %>% summarize(R.center = mean(act.R))
  r.sd <- inner_join(tr.stat, r.center, by='t', suffix=c('','.c')) %>%
    group_by(t) %>% summarize(R.sd = sd(act.R-R.center))
  r.csd <- inner_join(r.center, r.sd, by='t', suffix=c('','.sd'))

  ggplot(trs)+geom_smooth(aes(t,act.R), span=.2) +
    geom_hline(yintercept = r.mean, linetype='dashed') +
    geom_hline(yintercept = mean(r.center$R.center), linetype='dashed', color='red') +
    geom_line(data=r.csd, mapping=aes(t,R.center+R.sd), linetype='dotted', size=.5) +
    geom_line(data=r.csd, mapping=aes(t,R.center-R.sd), linetype='dotted', size=.5) +
    ggtitle(paste0(titleSfx,'Reward'))
  ggsave(paste0(fileSfx,'_Reward_smooth','.png'))

  ggplot(trs %>% filter(stid==19) %>% mutate(i = as.numeric(as.character(i))))+
    geom_line(aes(i,res.SumR)) +
    geom_smooth(aes(i,res.SumR), span=.5) +
    ggtitle(paste0(titleSfx,'final SumReward'))
  ggsave(paste0(fileSfx,'_SumReward_Tail','.png'))

  ggplot(trs %>% filter(stid==19) %>% mutate(i = as.numeric(as.character(i))))+
    geom_smooth(aes(i,res.SumR), span=.2)+
    ggtitle(paste0(titleSfx,'SumReward'))
  ggsave(paste0(fileSfx,'_SumReward','.png'))

  ggplot(trs %>% filter(stid==19) %>% mutate(i = as.numeric(as.character(i))) %>%
           gather(w,val,matches('W.S.'),-i))+
    geom_line(aes(i,val,color=w)) +
    geom_smooth(aes(i,val,color=w), linetype='dashed', span=.3, size=.5) +
    ggtitle(paste0(titleSfx,'final W'))
  ggsave(paste0(fileSfx,'_W_Tail','.png'))

}

runOptim <- function(xts.fs) {
  ofs <- xts.fs[5:nrow(xts.fs),]

  doLag <- function(d.t) {
    dt<-round(d.t);
    sum((ofs$Si[1:(nrow(ofs)-dt)]-ofs$S.w[(dt+1):nrow(ofs)])^2, na.rm=T)
  }
  o <- optim(c(0), doLag, method = 'Brent', lower=0, upper=5)
  o$par
}

printTracks <- function(res, env = res$myEnv, tracks = res$tracks) {
  titleSfx <- paste0('ICA',env$numSigs,'/',env$numOrig,': ')
  fileSfx <- paste0('xts_',env$numSigs,'_',env$numOrig)

  tracks <- tracks %>% mutate(t = stid) %>%
    mutate_at(vars(matches(regex('^res\\..*R\\d?$', ignore_case = F))),
                               list(~./env$xColNum))%>%
    left_join(env$X %>% set_names(paste0('X.',names(.))), by=c('t'='X.xid'))

  tsi <- tracks %>% select(t,i,matches('^SVal\\.S.$')) %>% gather(si,v,-t,-i)
  sumR <- tracks %>% filter(t==19) %>%
    rename(SumR = res.SumR) %>%
    select(i,SumR)
  w.s <- colMeans(env$log %>% tail(200) %>% select(matches('W\\.S.')))
  W <- env$W
  fs <- env$FS %>% rowid_to_column('t')
  for (n in env$xColNames) fs[paste0(n,'.w')] <- w.s[paste0('W.',n)]*env$X[n]
  fs$S.w <- rowSums(fs %>% select(matches('^S(\\d|C)\\.w$')))

  tailNum <- 500*(env$tSteps-1)
  tr.stat <- tail(tracks, tailNum) %>% filter(t > 4)
  r.mean <- mean(tail(tracks, tailNum)$act.R)
  r.center <- tr.stat %>% group_by(t) %>% summarize(R.center = mean(act.R))
  r.sd <- inner_join(tr.stat, r.center, by='t', suffix=c('','.c')) %>%
    group_by(t) %>% summarize(R.sd = sd(act.R-R.center))
  r.csd <- inner_join(r.center, r.sd, by='t', suffix=c('','.sd'))

  numVar <- ncol(fs)
  d <- as_datetime(1000*fs$t)
  xts.fs <- xts::xts(as.matrix(fs), order.by = d)
  xts.fs <- xts::merge.xts(xts.fs, as_datetime(1000*seq(1,max(fs$t),1/4)), fill=NA)
  xts.fs <- zoo::na.spline(xts.fs[,1:numVar])
  png(paste0(fileSfx,'_timeline.png'))
  plot(xts.fs[,c('S.w','Si')])
  dev.off()

  lag <- runOptim(xts.fs)
  xts.fs.lag <- xts.fs
  xts.fs.lag$Si <- xts::lag.xts(xts.fs.lag$Si, lag)
  png(paste0(fileSfx,'_timeline_lag.png'))
  plot(xts.fs.lag[,c('S.w','Si')])
  dev.off()

  xts.fs.lag.diff <- diff(xts.fs)
  png(paste0(fileSfx,'_timeline_lag_diff.png'))
  plot(xts.fs.lag.diff[,c('S.w','Si')])
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
  print(suppressWarnings(corrr::correlate(fs[5:(nrow(fs)-5),snames])))

  print(paste(titleSfx,'Correlation of Signals - Lag =',lag,':'))
  print(suppressWarnings(corrr::correlate(xts.fs.lag[5:(nrow(fs)-5),snames], quiet = T)))

  print(paste(titleSfx,'Correlation of Diff-Signals - Lag =',lag,':'))
  print(suppressWarnings(corrr::correlate(xts.fs.lag.diff[5:(nrow(fs)-5),snames], quiet = T)))

  #print(sumR %>% arrange(desc(SumR)), n=20)
}

icns <- c(2,3,4,5)
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

        print(paste('loading:',file))
        try({
          load(file)
          print(paste('loaded:',file))
          if (test) next

          assign(paste0('x',icn),myEnv)

          #myEnv <- get(paste0('x',icn))
          result <- getTrack(myEnv)
          result$tracks <- tracks
          if (icn < 5) assign(paste0(dn,'r',icn),result)
          save(result, file=paste0('finDemo.env_',host,'_ica_',icn,'_Res.rdata'))

          plotResults(result)
          evalLog(result)
          evalTracks(result)
        })
      }
    },
      error = function(err) { print(paste('Error:',dn,'/',icn,':',toString(err))) },
      finally = setwd(currDir))
  }
}

if (evalPrint) {
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

