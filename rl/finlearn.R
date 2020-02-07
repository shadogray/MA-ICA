source('FinanceDemo.R')

chooseAction <- function(myS, myEnv = finDemo, epsilon = myEnv$currEpsilon) {
  if (deblog) { print(paste('chooseAction:', toString(list(myS)))) }
  greedy <- runif(1) > epsilon

  topQs <- possibleActions(myEnv, myS)
  if (is.null(topQs)) return(NULL)

  if (greedy) {
    topQ <- topQs %>% arrange(desc(Q)) %>% slice(1)
  } else {
    topQ <- topQs[round(runif(1,1,nrow(topQs))),]
  }
  if (nlog) print(paste('chosen:', toString(list(topQ)), ', S=',toString(list(myS))))
  return (topQ)
}

takeAction <- function(A, myEnv = finDemo) {
  if (F&&deblog) { print(paste('takeAction:', toString(list(A)))) }
  resEnv <- act(myEnv, A$sid, A$aid)
  if (nlog) print(paste0('takeAction Res: resA=', toString(list(resEnv$resA)), ', logs=',toString(dim(resEnv$log))))
  return(resEnv)
}

control <- function(currAct, R, nextAoS, myEnv = finDemo, type = 'Q') {
  if (deblog) { print(paste('control: ', toString(list(currAct, R, nextAoS)))) }

  myEnv <- updateQ(myEnv, currAct, R, nextAoS, alpha = alpha, gamma = gamma)

  W <- as.numeric(myEnv$W)
  if (any(is.na(W) | is.nan(W) | W == Inf | W == -Inf)) { stop(paste('invalid W:', toString(W))) }
  return (myEnv)
}

learn <- function(myEnv = finDemo, numIter = 100, type = 'SARSA') {

  saveFileName <- paste0('finDemo.env_',Sys.info()['nodename'],ifelse(myEnv$ica,paste0('_ica_',myEnv$numSigs),''),'.rdata')
  lastSumR <- 0
  meanDiff <- 0
  iStart <- 1
  tracks <- tibble()
  myEnvs <- list()
  myEnv$run <- 0

  if (file.exists(saveFileName)) {
      load(saveFileName)
      print(paste('loaded:',saveFileName))
      #iStart <- myEnv$i
      myEnv <- init(myEnv)
      numIter <- ifelse(iStart < numIter, numIter, iStart + numIter)
  }
  myEnv$run <- myEnv$run + 1

  print(paste0('Environment Size: ',object.size(myEnv)))
  print(paste0('States Size     : ',object.size(myEnv$S)))
  print(paste0('Actions Size    : ',object.size(myEnv$A)))
  print(paste0('Environments    : ',length(myEnvs)))

  for (i in iStart:numIter) {
    print(paste('run:',myEnv$run,'iter:',i,'time:',Sys.time(),'type:',type,'logs:',nrow(myEnv$log)))
    myEnv$Iterations <- i
    myEnv$i <- i

    if (type == 'SARSA') {

      myEnv <- getStart(myEnv)
      currState <- getState(myEnv)
      currAct <- chooseAction(currState, myEnv)
      myEnv$chosenAction <- currAct

      while(TRUE) {
        myEnv <- takeAction(currAct, myEnv)
        resA <- myEnv$resA
        R <- resA$R

        if (resA$A$terminal || (has_name(resA$A,'failed') && resA$A$failed)) {
          myEnv <- control(resA$A, R, NULL, type = type, myEnv = myEnv)
          break
        }

        nextState <- resA$S

        if (is.na(nextState$sid)) {
          R <- -100
          myEnv <- control(resA$A, R, NULL, type = type, myEnv = myEnv)

          nextAct <- chooseAction(currState, myEnv)
          if (is.null(nextAct)) { print(paste0('Failed: S=',toString(currState),', SValues=',toString(myEnv$SValues))); break }

          myEnv$chosenAction <- nextAct
          currAct <- nextAct
          next
        } else {
          nextAct <- chooseAction(nextState, myEnv)
          if (is.null(nextAct)) { print(paste0('Failed: S=',toString(nextState),', SValues=',toString(myEnv$SValues))); break }

          myEnv$chosenAction <- nextAct
        }

        myEnv <- control(resA$A, R, nextAct, type = type, myEnv = myEnv)

        if (deblog) print(tibble('while', bind_rows(getStateAction(myEnv, currState$sid, currAct$aid),
                          getStateAction(myEnv, nextState$sid, nextAct$aid))), n=Inf, width=Inf)

        currState <- nextState
        currAct <- nextAct
      }

    }
    if (nlog) print(paste('W=',toString(myEnv$W)))

    if (i%%10 == 0) {
      res <- findWay(myEnv)
      tracks <- bind_rows(tracks, res$track %>% mutate(i = myEnv$i))
      myEnv$tracks <- tracks
    }

    if (i%%(100/myEnv$numSigs) == 0) {
      myEnvs[[myEnv$run]] <- myEnv
      save(myEnv, tracks, myEnvs, file = saveFileName)
    }


    meanDiff <- (meanDiff*(i-1) + abs(lastSumR - myEnv$SumR)) / i
    if (myEnv$SumR > 0 && abs(meanDiff/myEnv$SumR) < deltaSumR) {
      print(paste('END: run:',myEnv$run,'SumR:',toString(myEnv$SumR),'diffSumR:',toString(lastSumR - myEnv$SumR),'meanDiff:',toString(meanDiff)))
      myEnvs[[myEnv$run]] <- myEnv
      save(myEnv, tracks, myEnvs, file = saveFileName)
      #break
    } else {
      if (nlog) print(paste('SumR:',toString(myEnv$SumR),'diffSumR:',toString(lastSumR - myEnv$SumR),'meanDiff:',toString(meanDiff)))
      lastSumR <- myEnv$SumR
    }
  }

  print(paste('END: run:',myEnv$run,'SumR:',toString(myEnv$SumR)))
  myEnv
}

findWay <- function(myEnv = finDemo) {

  myEnv <- getStart(myEnv)
  currState <- getState(myEnv)
  track <- tibble()
  n.t <- microbenchmark::get_nanotime()

  for (i in 1:myEnv$tSteps+20) {

    nextAct <- chooseAction(currState, myEnv, epsilon = 0)
    if (is.null(nextAct)) { print(paste0('Failed: S=',toString(currState),', SValues=',toString(myEnv$SValues))); next }
    if (deblog) n.t <- nanoTime(paste('nextAct',myEnv$numSigs), n.t)

    if (deblog) { print('Track:'); print(track) }

    myEnv <- takeAction(nextAct, myEnv)
    if (deblog) n.t <- nanoTime(paste('takeAction',myEnv$numSigs), n.t)
    if (is.null(myEnv$resA)) {
      break # we are done?
    }

    resA <- myEnv$resA
    nextState <- resA$S

    for (r in 1:length(resA$R.r)) resA$A[paste0('R',r)] <- resA$R.r[r]
    resA$A$R <- sum(resA$R.r)
    resA$A$SumR <- sum(myEnv$SumR)
    for (r in 1:length(myEnv$SumR)) resA$A[paste0('SumR',r)] <- myEnv$SumR[r]

    if (deblog) n.t <- nanoTime(paste('move',myEnv$numSigs), n.t)

    track <- track %>% bind_rows(
      bind_cols(currState, myEnv$SValues %>% set_names(paste0('SVal.',names(.)))) %>%
        bind_cols(nextAct %>% set_names(paste0('act.',names(.)))) %>%
        bind_cols(resA$A %>% set_names(paste0('res.',names(.)))) %>%
        bind_cols(nextState %>% set_names(paste0('next.',names(.)))) %>%
        bind_cols(as.data.frame(t(myEnv$W)) %>% set_names(paste0('W.',myEnv$xColNames))) %>%
        mutate(time=microbenchmark::get_nanotime(), t=myEnv$t)
    )

    if (resA$A$terminal) {
      print(paste0('END[',myEnv$t,']: resA=',toString(resA)))
      if (nlog) print(currState %>% bind_rows(nextAct, nextState))
      break
    }
    if (identical(currState, nextState)) {
      print('same state found:  ')
      print(currState %>% bind_rows(nextAct, nextState))
      #break # we are done?
    }

    if (deblog) print(currState %>% bind_rows(nextState))
    currState <- nextState
  }
  list(track = track %>% rowid_to_column('stid'), myEnv = myEnv)
}
