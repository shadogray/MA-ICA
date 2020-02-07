library(tidyverse)
library(lubridate)
source('finlearn.R')

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
  print(xts.fs)
  png('timeline.png')
  print(xts::plot.xts(xts.fs))
  dev.off()
  
}

load('ssteps44.rdata')
printTracks(res = result, env = result$myEnv, tracks = result$tracks)