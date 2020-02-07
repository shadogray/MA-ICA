library(tidyverse)

source('worlds.R')
source('finlearn.R')

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
args <- commandArgs(trailingOnly = T)
if (length(args) > 0) type = args[1]

#######

getTrack <- function(env = myEnv) {
  env$SumR <- 0
  result <- findWay(env)
  result$b <- tibble(B1 = env$B1, B2 = env$B2)
  if (!is.null(env$B3)) result$b <- bind_cols(result$b, B3=env$B3)
  if (!is.null(env$B4)) result$b <- bind_cols(result$b, B4=env$B4)
  result$track <- result$track %>% rowid_to_column('t')
  result
}

printTrack <- function(result = result, env = result$myEnv, track = result$track) {
  print(track, n=nrow(track))
  print(paste('Type:',type,
              'Iterations:',env$Iterations,
              'FinalR:', track %>% filter(res.terminal) %>% select(res.R),
              'SumR:',toString(env$SumR),'W:',toString(env$W)))
}

plotResults <- function(result = result, env = result$myEnv, track = result$track) {

  numOrig <- ifelse(is.null(env$numOrig), 3, env$numOrig)
  if (is.null(env$numBase)) numOrig <- 2
  b <- as_tibble(as.matrix(result$b)*env$baseX)
  W <- env$W / sum(abs(env$W))

  titleSfx <- paste0('ICA',env$numSigs,'/',numOrig,': ')
  fileSfx <- paste0('S_density_ica_',env$numSigs,'_',numOrig)
  comb <- gtools::combinations(env$xColNum,2)
  for (i in 1:nrow(comb)) {
    si <- paste0('S',comb[i,])
    ggplot(env$S %>% select(si,n,P) %>% group_by_at(c(si,'P')) %>% summarize(n=sum(n))) +
      geom_contour(aes_string(si[1],si[2], z='n'))+facet_wrap(~P) +
      ggtitle(paste0(titleSfx,'State-Evaluation Density - States ',si[2],'~',si[1]))
    ggsave(paste0(fileSfx,'_',i,'.png'))
  }

  ggplot(bind_cols(b,FS=env$FS$Si) %>%
           rowid_to_column('tx') %>% gather(signal,value,-tx))+
    geom_line(aes(tx,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Base Timeseries (Bi) and Target Sum (FS)')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_BaseSignals','.png'))

  x.x <- as_tibble(t(apply(as.matrix(env$X[env$xColNames]),1,`*`,1)) * env$FS$p) %>%
    mutate(S = rowSums(.))
  ggplot(bind_cols(x.x, FS = env$FS$Si) %>%
    rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))+
    ggtitle(paste0(titleSfx,'Target Signal (FS), Price-weighted ICA Signals (Si * p(t)) and Sum of Si (S)')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Signals','.png'))

  x.Wx <- as_tibble(t(apply(as.matrix(env$X[env$xColNames]),1,`*`,as.numeric(W))) * env$FS$p) %>%
    mutate(S = rowSums(.))
  ggplot(bind_cols(x.Wx, FS = env$FS$Si) %>%
           rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Price-and-W weighted ICA Signals (Si * W * p(t))')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_WeightedSignals','.png'))

  ggplot(as_tibble(env$X %>% select(env$xColNames)) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+
    geom_line(aes(tx,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Raw ICA Signals (Si)')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_RawSignals','.png'))

  ggplot(as_tibble(env$X.orig) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Artificial TimeSeries generated from BaseSignals Bi')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_ArtificialSeries','.png'))

  ggplot(b %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  ggplot(select(env$FS,Si) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  ggplot(select(env$FS,S) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  ggplot(select(env$FS,p) %>% rowid_to_column('tx') %>% gather(signal,value,-tx))+geom_line(aes(tx,value,color=signal))
  #ggplot(env$A %>% filter(n>0)) + geom_density_2d(aes(sid,aid))
  #ggplot(env$A %>% filter(n>0) %>% group_by(sid,aid) %>% summarise(n = sum(n))) + geom_point(aes(sid,aid,size=n))


  track$F.Si <- env$FS$Si[1:nrow(track)]
  for (n in env$vColNames) track[paste0(n,'.p')] <- track[n] * env$FS$p
  track$S.p <- rowSums(track %>% select(matches('^S\\di$'))) * env$FS$p[1:nrow(track)]

  ggplot(track %>% select(t,matches('^S\\di.p$')) %>% gather(signal,value,-t))+
    geom_line(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Resulting State-Track following generated Policy'), ) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_StateTrack_Si','.png'))

  ggplot(track %>% select(t,matches('^S.*i\\.p$'),F.Si) %>% gather(signal,value,-t))+
    geom_line(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Resulting State-Track following generated Policy'), ) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_StateTrack','.png'))

  ggplot(track %>% select(t,matches('^res\\..*R$')) %>% gather(signal,value,-t))+geom_line(aes(t,value,color=signal)) +
    ggtitle(paste0(titleSfx,'Resulting Reward following generated Policy')) +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(fileSfx,'_Reward','.png'))

  #ggplot(track %>% select(t,matches('S\\di1')) %>% gather(signal,value,-t))+geom_line(aes(t,value,color=signal))

  #ggplot(env$S %>% select(matches('S\\d'),P,n) %>% gather(key,value,-P,-n)) +
  #ggplot(env$S %>% select(matches('S\\d'),n,P)) + geom_contour(aes_string(x='S1',y='S2', z='n'))+facet_wrap(~P)
}

run_many <- function(env = myEnv) {
  for (i in 1:9) {
    assign(paste0('res',i), getTrack(env))
  }
  res <- gather(res1$track %>% select(matches('S\\di')), key, value)
  for (i in 2:9) {
    res <- bind_cols(res, gather(get(paste0('res',i))$track %>% select(matches('S\\di')), key, value) %>% select(value))
  }
  ggplot(res %>% rowid_to_column('t') %>% gather(k,v,-t,-key)) + geom_line(aes(t,v,color=k))
}

for (icn in c(2,3,4,5)) {
  load(paste0('finDemo.env_storage_ica_',icn,'.rdata'))
  result <- getTrack()
  result$tracks <- tracks
  printTrack(result)
  plotResults(result)
}
