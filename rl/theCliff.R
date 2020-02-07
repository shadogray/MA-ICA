##############################################################
######################### The Cliff ##########################

theCliff <- function(X, Y,
                     A = factor(c('up','down','left','right')),
                     offLimR = -1, cliffR = -100, ...) {
  new_theCliff(X, Y, A, offLimR, cliffR, ...)
}

new_theCliff <- function(X, Y, A, offLimR, cliffR, ..., class = character()) {
  minX = min(X)
  maxX = max(X)
  minY = min(Y)
  maxY = max(Y)
  S = as_tibble(expand.grid(Sx = X, Sy = Y)) %>% rowid_to_column('sid')
  #    filter(!trim | !(Sx==minX & A=='down')) %>%
  #    filter(!trim | !(Sx==maxX & A=='up')) %>%
  #    filter(!trim | !(Sy==minY & A=='left')) %>%
  #    filter(!trim | !(Sx==maxY & A=='right')) %>%
  Q = as_tibble(expand.grid(sid = S$sid, A = A)) %>%
    mutate(Q = 0) %>% rowid_to_column('aid')
  new_world(
    X = X, Y = Y, A = A,
    minX = minX, maxX = maxX, minY = minY, maxY = maxY,
    S = S, Q = Q, SumR = .0, offLimR = offLimR, cliffR = cliffR,
    class = c(class, 'theCliff')
  )
}

#' initialize environment (e.g. reset Q to zero)
init.theCliff <- function(x, ...) {
  x$Q$Q <- .0
  x$SumR <- .0
  x
}

#' get the starting state of this environment
getStart.theCliff <- function(x, ...) {
  x$S %>% filter(Sx == x$minX, Sy == x$minY)
}

#' find the real state for the possibly updated location values (Sx,Sy,..) of state S
#' @param S the possbly updated input state
#' @return the real state element matching the example state
getState.theCliff <- function(x, S, ...) {
  if (is.numeric(S)) return (x$S[S,])
  return (x$S %>% filter(Sx == S$Sx, Sy == S$Sy))
}

#' get the real State/Action for the possibly updated S and A
#' @param state S as tibble/list
#' @param action A as tibble/list or factor
#' @return a merged object of S and A
#' @seealso getState
getStateAction.theCliff <- function(x, S, A, ...) {
  if (is.numeric(S))
    s <- x$S[S,]
  else
    s <- getState(x, S)
  if (is.numeric(A)) {
    a <- x$Q[A,]
  } else {
    mya <- A
    if (!(is.factor(A) || is.character(A))) mya <- A$A
    a <- x$Q %>% filter(sid == s$sid, A == mya)
  }
  bind_cols(s, a)
}

#' get possible actions (tibble) for the state S
#' @param state S as tibble/list or numeric (sid)
#' @return tibble of actions
#'
possibleActions.theCliff <- function(x, S) {
  if (is.numeric(S)) {
    id <- S
  } else {
    S <- getState(x, S)
    id <- S$sid
  }
  x$Q %>% filter(sid == id)
}

#' Updates state S according action A (as tibble/list or as factor/character)
#' @param state S as tibble/list
#' @param action A as tibble/list or factor
#' @return updated S does not call getState
evalAction.theCliff <- function(x, S, A, ...) {
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$Q[A,]
  a <- A
  if (!(is.factor(A)||is.character(A))) a <- A$A
  if (a == 'right') {
    S$Sx <- S$Sx + 1
  } else if (a == 'left') {
    S$Sx <- S$Sx - 1
  } else if (a == 'up') {
    S$Sy <- S$Sy + 1
  } else if (a == 'down') {
    S$Sy <- S$Sy - 1
  }
  return (S)
}

#' act out the action A in state S
#' @param state S as tibble/list
#' @param action A as tibble/list or factor/character
#' @return the new state S' augmented by reward R
#' @seealso evalAction
#' @seealso getState
act.theCliff <- function(x, S, A, ...) {
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$Q[A,]
  Sin <- S
  Ain <- A
  S <- evalAction(x, S, A)

  R <- -1 # on all normal transitions
  if (S$Sx == x$maxX && S$Sy == x$minY) {
    return (NULL) # end of episode
  }
  if (S$Sx > x$minX && S$Sy == x$minY) { # the Cliff!!
    R <- x$cliffR
    S$Sx <- x$minX
    S$Sy <- x$minY
  } else if (S$Sx < x$minX || S$Sx > x$maxX ||
             S$Sy < x$minY || S$Sy > x$maxY) {
    R <- x$offLimR # add penalty for Off-Track
    S$Sx <- Sin$Sx
    S$Sy <- Sin$Sy
  }

  S <- getState(x, S)
  S$R <- R
  return(S)
}

#' check if acting out action A in state S will leave the world's context
#' @param state S as tibble/list
#' @param action A as tibble/list or factor/character
#' @return true if next state would be outside world's context
#' @seealso evalAction
offLimits.theCliff <- function(x, S, A, ...) {
  S <- evalAction(x, S, A)
  is.null(S) || nrow(S) == 0
}

#' update Q for this action A
#' @param action A as tibble/list (match is done by A$aid) or as aid
#' @return the environment
updateQ.theCliff <- function(x, A, q, ...) {
  if (is.numeric(A))
    aid <- A
  else
    aid <- A$aid
  x$Q$Q[aid] <- q
  x
}

#' test core functions of this world
test.theCliff <- function(x, ...) {
  envTheCliff <- theCliff(1:12, 1:10)
  print(sloop::s3_dispatch(possibleActions(envTheCliff, list())))
  startState <- getStart(envTheCliff)
  assertthat::assert_that(startState$Sx==1 && startState$Sy==1)
  startState <- getState(envTheCliff, startState)
  assertthat::assert_that(startState$Sx==1 && startState$Sy==1)
  nextState <- startState
  nextState$Sy <- nextState$Sy + 1
  upState <- getState(envTheCliff, nextState)
  assertthat::assert_that(upState$Sx==nextState$Sx && upState$Sy==nextState$Sy)

  poss <- possibleActions(envTheCliff, startState)
  print(poss)
  assertthat::assert_that(identical(poss, possibleActions(envTheCliff, startState$sid)))
  assertthat::assert_that(nrow(poss) == length(envTheCliff$A))

  updState <- evalAction(envTheCliff, startState, poss[poss$A=='up',])
  assertthat::assert_that(updState$Sx==1 && updState$Sy==2)
  updState <- evalAction(envTheCliff, startState$sid, poss[poss$A=='up',]$aid)
  assertthat::assert_that(updState$Sx==1 && updState$Sy==2)

  deadAct <- act(envTheCliff, startState, poss[poss$A=='right',])
  assertthat::assert_that(deadAct$R == -100)
  deadAct <- act(envTheCliff, startState$sid, poss[poss$A=='right',]$aid)
  assertthat::assert_that(deadAct$R == -100)
  deadAct <- act(envTheCliff, startState, 'right')
  assertthat::assert_that(deadAct$R == -100)

  SA <- getStateAction(envTheCliff, startState, 'up')
  assertthat::assert_that(SA$sid == startState$sid && SA$A == 'up')
  SA <- getStateAction(envTheCliff, startState$sid, SA$aid)
  assertthat::assert_that(SA$sid == startState$sid && SA$A == 'up')

  env2 <- updateQ(envTheCliff, SA, 33)
  assertthat::assert_that(env2$Q$Q[env2$Q$aid == SA$aid] == 33)
  env2 <- updateQ(envTheCliff, SA$aid, 44)
  assertthat::assert_that(env2$Q$Q[env2$Q$aid == SA$aid] == 44)
}

######################### The Cliff ##########################
##############################################################