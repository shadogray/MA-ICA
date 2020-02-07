library(tidyverse)
library(lubridate)
source('finlearn.R')

sig.v <- NULL
orig.s <- NULL
node <- Sys.info()['nodename']

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
## load('ssteps33.rdata');res=result;env=result$myEnv;tracks=result$tracks
#######

evalLogSuccess <- function(res, env, iter, tag) {
  
  fileExt <- paste0('_',env$numSigs,'_',env$numOrig,tag,'.png')
  
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
  res$lwtin <- sum(lwtin)
  res$lwiter <- nrow(lw)
  
  
  trs <- env$tracks %>% tail(20000)
  tr.stat <- trs %>% filter(t > 5)
  r.mean <- mean(trs$act.R)
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
  ggsave(paste0('Reward_Tracks_smooth',fileExt), dpi = 'print', width = 15, height = 10, units = 'cm')
  
  res$r.mean <- r.mean
  res$r.centerMean <- r.centerMean
  
  res
}


icns <- c(2:10)
if (T) {
  print(paste('Running evalPlot:',getwd()))
  tryCatch({
    succ <- tibble(numOrig=numeric(),numSigs=numeric(),iter=numeric(),
                   lwc=numeric(), lwtin=numeric(), lwiter=numeric(),
                   r.mean=numeric(), r.centerMean=numeric())
    
    for (icn in icns) {
      file <- paste0('finDemo.env_monster_ica_',icn,'_Res.rdata')

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

        for (i in 1:length(result$myEnvs)) {
          #evalLog(result, result$myEnvs[[i]], i, paste0('_I',i))
          res <- evalLogSuccess(result, result$myEnvs[[i]], i, paste0('_I',i))
          e <- res$myEnv
          succ <- add_row(succ, numOrig=e$numOrig, numSigs=e$numSigs, iter=i, 
                          lwc=res$lwc, lwtin=res$lwtin, lwiter=res$lwiter,
                          r.mean=res$r.mean, r.centerMean=res$r.centerMean)
        }
        
        save(succ, file = paste0('successrate_',icn,'.rdata'))
      })
    }
  },
  error = function(err) { print(paste('Error:',icn,':',toString(err))) },
  finally = {})
}

