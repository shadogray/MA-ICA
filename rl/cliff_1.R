library(tidyverse)

X <- 1:12
Y <- 1:5
minX = min(X)
maxX = max(X)
minY = min(Y)
maxY = max(Y)
A <- factor(c('up','down','left','right'))
Q <- as_tibble(expand.grid(Sx = X, Sy = Y, A = A)) %>%
#  filter(!(Sx==minX & A=='down')) %>% filter(!(Sx==maxX & A=='up')) %>%
#  filter(!(Sy==minY & A=='left')) %>% filter(!(Sx==maxY & A=='right')) %>%
  mutate(Q = 0) %>% rowid_to_column('id')

epsilon = .1
alpha = .1
gamma = 1
ctrl.types = factor('Q','SARSA')
type = 'SARSA' # 'Q' #'SARSA'

envTheCliff <- function(actQ, myQ = Q) {
  myAct <- actQ
  if (myAct$A == 'right') myAct$Sx <- myAct$Sx + 1
  else if (myAct$A == 'left') myAct$Sx <- myAct$Sx - 1
  else if (myAct$A == 'up') myAct$Sy <- myAct$Sy + 1
  else if (myAct$A == 'down') myAct$Sy <- myAct$Sy - 1

  myAct$R <- -1 # on all normal transitions
  if (myAct$Sx == maxX && myAct$Sy == minY) {
    return (NULL) # end of episode
  }
  if (myAct$Sx > minX && myAct$Sy == minY) { # the Cliff!!
    myAct$R <- -100
    myAct$Sx <- minX
    myAct$Sy <- minY
  } else if (myAct$Sx < minX || myAct$Sx > maxX || myAct$Sy < minY || myAct$Sy > maxY) {
    myAct$R <- -2 # add penalty for Off-Track
    myAct$Sx <- actQ$Sx
    myAct$Sy <- actQ$Sy
  }
  return(myAct)
}

chooseAction <- function(myS, myQ = Q, myEpsilon = epsilon) {
  #print(paste('myS=',myS$Sy,',',myS$Sy))
  greedy <- runif(1) > myEpsilon
  topQs <- myQ %>% filter(Sx == myS$Sx, Sy == myS$Sy)
  if (greedy) {
     topQ <- topQs %>% arrange(desc(Q)) %>% slice(1)
  } else {
    topQ <- topQs[round(runif(1,1,nrow(topQs))),]
  }
  return (topQ)
}

takeAction <- function(action, env = envTheCliff) {
  #print(paste('action=',action$Sy,',',action$Sy))
  res <- env(action)
  #print(paste('res=',res$Sy,',',res$Sy,',',res$A))
  return(res)
}

control <- function(currAct, nextAct, myQ = Q, type = 'Q') {
  #print(paste('curr=',currAct$Sy,',',currAct$Sy,'next=',nextAct$Sy,',',nextAct$Sy))
  if (type == 'Q') {
    # the max Q for all actions in next state
    maxNextQ <- myQ %>% filter(Sx == nextAct$Sx, Sy == nextAct$Sy) %>%
      arrange(desc(Q)) %>% slice(1)
    q <- currAct$Q + alpha * (currAct$R + gamma * maxNextQ$Q - currAct$Q)
  } else if (type == 'SARSA') {
    # the max Q for all actions in next state
    currQ <- myQ %>% filter(Sx == currAct$Sx, Sy == currAct$Sy, A == currAct$A) %>% slice(1)
    nextQ <- myQ %>% filter(Sx == nextAct$Sx, Sy == nextAct$Sy, A == nextAct$A) %>% slice(1)
    q <- currAct$Q + alpha * (currAct$R + gamma * nextQ$Q - currAct$Q)
  }
  myQ <- myQ %>% mutate(Q = ifelse(id == currAct$id, q, Q))
  #print(myQ %>% filter(id==currAct$id))
  return (myQ)
}

Q$Q <- 0

for (i in 1:500) {
  if (type == 'Q') {
    currS <- tibble(Sx=minX, Sy=minY)
    while(TRUE) {
      act <- chooseAction(currS)
      nextAct <- takeAction(act)
      if (is.null(nextAct)) break

      act$R <- nextAct$R
      Q <- control(act, nextAct = nextAct, type = type)

      #print(currS %>% bind_rows(act,nextAct))
      currS <- nextAct
    }
    #print(currS %>% bind_rows(act,nextAct))

  } else if (type == 'SARSA') {

    currState <- Q %>% filter(Sx==minX, Sy==minY) %>% head(1)
    currState <- chooseAction(currState)
    while(TRUE) {
      actResult <- takeAction(currState)
      if (is.null(actResult)) break

      nextState <- actResult
      nextState <- chooseAction(nextState)

      currState$R <- actResult$R
      Q <- control(currState, nextAct = nextState, type = type)

      #print(currState %>% bind_rows(actResult,nextState))

      currState <- nextState
    }
    print(currAct %>% bind_rows(actAct,nextAct))
  }
  print(Q %>% filter(Q != 0) %>% arrange(desc(Q)))
}

findWay <- function(myQ = Q) {
  currS <- tibble(Sx=minX, Sy=minY)
  while (TRUE) {
    act <- chooseAction(currS, myEpsilon = 0)
    nextAct <- takeAction(act)
    print(act %>% bind_rows(nextAct))
    if (is.null(nextAct)) break
    currS <- nextAct
  }
}
findWay()
