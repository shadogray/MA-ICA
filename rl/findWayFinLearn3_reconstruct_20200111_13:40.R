library(tidyverse)
library(lubridate)
source('finlearn.R')

sig.v <- NULL
orig.s <- NULL
test <- F
evalPlot <- T
evalPrint <- F
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
load('ssteps33.rdata');res=result;env=result$myEnv;tracks=result$tracks

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
  ggplot(l %>% gather(Weight,Value,matches('^[W]\\.*'),-t))+geom_line(aes(i,Value,color=Weight))
  ggsave(paste0(fileSfx,'_i_W',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  ggplot(l %>% gather(Weight,Value,matches('^[W]\\.*'),-t))+geom_smooth(aes(i,Value,color=Weight))
  ggsave(paste0(fileSfx,'_i_W_smooth',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
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

}

evalLogSuccess <- function(res, env = res$myEnv, tag = '') {
  
  # detect approximation success
  lw <- env$log %>% select(matches('^W.*'),lid)
  lwtm <- lw %>% select(matches('^W')) %>% tail(200) %>% summarize_all(mean)
  lwt <- matrix(nrow = nrow(lw),ncol = ncol(lw)-1)
  for (i in 1:(ncol(lw)-1)) {
    lwtma <- abs(as.numeric(lwtm))
    lwfloor <- lwtma-.1
    lwceil <- lwtma+.1
    lwt[,i] <- between(abs(pull(lw,i)), lwfloor[i], lwceil[i])
  }
  lwtin <- apply(lwt,1,all)

  #ggplot(lw %>% filter(lwtin) %>% gather(k,v,-lid))+geom_line(aes(lid,v,color=k))
  #ggplot(lw %>% gather(k,v,-lid))+geom_line(aes(lid,v,color=k))
  lwc <- sum(lwtin)/length(lwtin)
  
  res$lwc <- lwc
  res
}


icns <- c(2:10)
baseDir <- list('v6sig2'='monster','v6sig3'='monster','v6sig4'='monster')

if (!is.null(sig.v)) {
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
      
      succ <- tibble(numOrig=numeric(),numSigs=numeric(),iter=numeric(),lwc=numeric())
      
      for (icn in icns) {
        file <- paste0('finDemo.env_',host,'_ica_',icn,'_Res.rdata')
        if (!file.exists(file)) host <- node
        file <- paste0('finDemo.env_',host,'_ica_',icn,'_Res.rdata')
        
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
          
          for (i in 1:length(result$myEnvs)) {
            #evalLog(result, result$myEnvs[[i]], paste0('_I',i,'_'))
            res <- evalLogSuccess(result, result$myEnvs[[i]], paste0('_I',i,'_'))
            e <- res$myEnv
            succ <- add_row(succ, numOrig=e$numOrig, numSigs=e$numSigs, iter=i, lwc=res$lwc)
          }
          
          save(succ, file = 'successrate.rdata')
        })
      }
    },
    error = function(err) { print(paste('Error:',dn,'/',icn,':',toString(err))) },
    finally = setwd(currDir))
  }
}

