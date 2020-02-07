library(tidyverse)

source('worlds_OK.R')

deblog = F
epsilon = .1
alpha = .1
gamma = 1
ctrl.types = factor('Q','SARSA')

chooseAction <- function(myS, myEnv = envTheCliff, myEpsilon = epsilon) {
  #if (deblog) print(paste('myS=',myS$Sy,',',myS$Sy))
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
  #if (deblog) print(paste('curr=',currAct$Sy,',',currAct$Sy,'next=',nextAct$Sy,',',nextAct$Sy))
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
    if (deblog) print(currState %>% bind_rows(nextState))
    currState <- nextState
  }
  print(track, n=nrow(track))
}

type = 'SARSA' # 'Q' #'SARSA'
envTheCliff <- theCliff(1:12,1:5)
envTheCliff$offLimR = -1
envTheCliff <- init(envTheCliff)

stTime <- Sys.time()
for (i in 1:500) {
  if (type == 'Q') {
    currState <- getStart(envTheCliff)
    while(TRUE) {
      nextAct <- chooseAction(currState)
      nextState <- takeAction(nextAct)
      if (is.null(nextState)) break

      R <- nextState$R
      envTheCliff <- control(nextAct, R, nextState, type = type)

      if (deblog) print(bind_rows(getStateAction(envTheCliff, currState$sid, nextAct$aid),nextState))

      envTheCliff$SumR <- envTheCliff$SumR + R
      currState <- nextState
    }

  } else if (type == 'SARSA') {

    currState <- getStart(envTheCliff)
    currAct <- chooseAction(currState)
    while(TRUE) {
      nextState <- takeAction(currAct)
      if (is.null(nextState)) break

      R <- nextState$R
      nextAct <- chooseAction(nextState)

      envTheCliff <- control(currAct, R, nextAct, type = type)

      if (deblog)
        print(bind_rows(getStateAction(envTheCliff, currState$sid, currAct$aid),
                        getStateAction(envTheCliff, nextState$sid, nextAct$aid)))

      envTheCliff$SumR <- envTheCliff$SumR + R
      currState <- nextState
      currAct <- nextAct
    }
  }
  if (deblog) print(envTheCliff$Q %>% filter(Q != 0) %>% arrange(desc(Q)))
}
print(paste('execution: ',Sys.time()-stTime,'SumR:',envTheCliff$SumR))

track <- findWay()

