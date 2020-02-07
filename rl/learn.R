chooseAction <- function(myS, myEnv = envTheCliff, myEpsilon = epsilon) {
  if (deblog) print(paste('myS=',myS$Sy,',',myS$Sy))
  greedy <- runif(1) > myEpsilon
  topQs <- possibleActions(myEnv, myS)
  if (greedy) {
    topQ <- topQs %>% arrange(desc(Q)) %>% slice(1)
  } else {
    topQ <- topQs[round(runif(1,1,nrow(topQs))),]
  }
  return (topQ)
}

takeAction <- function(A, myEnv = envTheCliff) {
  if (deblog) print(A)
  res <- act(myEnv, A$sid, A$aid)
  if (deblog) print(res)
  return(res)
}

control <- function(currAct, R, nextAoS, myEnv = envTheCliff, type = 'Q') {
  if (deblog) print(paste('curr=',currAct$Sy,',',currAct$Sy,'next=',nextAct$Sy,',',nextAct$Sy))
  if (type == 'Q') {
    # the max Q for all actions in next state
    maxNextQ <- possibleActions(myEnv, nextAoS$sid) %>% arrange(desc(Q)) %>% slice(1)
    q <- currAct$Q + alpha * (R + gamma * maxNextQ$Q - currAct$Q)
  } else if (type == 'SARSA') {
    q <- currAct$Q + alpha * (R + gamma * nextAoS$Q - currAct$Q)
  }
  myEnv <- updateQ(myEnv, currAct$aid, q)
  if (deblog) print(currAct %>% bind_rows(nextAoS, getStateAction(myEnv, currAct$sid, currAct$aid)))
  if (deblog) print(paste('q=',q))
  return (myEnv)
}

learn <- function(myEnv = envTheCliff) {

  lastSumR <- 0
  meanDiff <- 0

  for (i in 1:numIter) {
    myEnv$Iterations <- i

    if (type == 'Q') {
      currState <- getStart(myEnv)
      while(TRUE) {
        nextAct <- chooseAction(currState, myEnv = myEnv)
        nextState <- takeAction(nextAct, myEnv = myEnv)
        if (is.null(nextState)) break

        R <- nextState$R
        myEnv$SumR <- myEnv$SumR + R
        myEnv <- control(nextAct, R, nextState, type = type, myEnv = myEnv)

        if (deblog) print(bind_rows(getStateAction(myEnv, currState$sid, nextAct$aid),nextState))

        currState <- nextState
      }

    } else if (type == 'SARSA') {

      currState <- getStart(myEnv)
      currAct <- chooseAction(currState, myEnv = myEnv)
      while(TRUE) {
        nextState <- takeAction(currAct, myEnv = myEnv)
        if (is.null(nextState)) break

        R <- nextState$R
        myEnv$SumR <- myEnv$SumR + R
        nextAct <- chooseAction(nextState, myEnv = myEnv)

        myEnv <- control(currAct, R, nextAct, type = type, myEnv = myEnv)

        if (deblog)
          print(bind_rows(getStateAction(myEnv, currState$sid, currAct$aid),
                          getStateAction(myEnv, nextState$sid, nextAct$aid)))

        currState <- nextState
        currAct <- nextAct
      }

    }
    if (deblog) print(myEnv$Q %>% filter(Q != 0) %>% arrange(desc(Q)))

    meanDiff <- (meanDiff*(i-1) + abs(lastSumR - myEnv$SumR)) / i
    if (abs(meanDiff/myEnv$SumR) < deltaSumR) {
      print(paste('SumR:',myEnv$SumR,'diffSumR:',lastSumR - myEnv$SumR,'meanDiff:',meanDiff))
      break
    } else {
      if (deblog) print(paste('SumR:',myEnv$SumR,'diffSumR:',lastSumR - myEnv$SumR,'meanDiff:',meanDiff))
      lastSumR <- myEnv$SumR
    }
  }

  myEnv
}

findWay <- function(myEnv = envTheCliff) {
  currState <- getStart(myEnv)
  track <- currState

  while (TRUE) {

    nextAct <- chooseAction(currState, myEnv = myEnv, myEpsilon = 0)
    track <- track %>% bind_rows(bind_cols(currState, nextAct))

    if (deblog) print(track)

    nextState <- takeAction(nextAct, myEnv = myEnv)

    if (is.null(nextState)) {
      break # we are done?
    }
    if (identical(currState, nextState)) {
      print('no useful state found:  ')
      print(currState %>% bind_rows(nextAct, nextState))
      break # we are done?
    }

    myEnv$SumR <- myEnv$SumR + nextState$R

    if (deblog) print(currState %>% bind_rows(nextState))
    currState <- nextState
  }
  list(track = track, myEnv = myEnv)
}
