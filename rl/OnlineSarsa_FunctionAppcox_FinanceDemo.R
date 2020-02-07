##############################################################
source('worlds.R')
######################### Finance Demo ##########################
sSteps = 20; tSteps = 30; baseX = 100; X = NULL; deltaX = .1; deflM = .05; txFee = .1; ica=T; numSigs=2; numBase = 7; numOrig = 4

finDemo <- function(sSteps = 10, tSteps = 20, baseX = 100, alpha = .1, gamma = .9, X = NULL, deltaX = .1,
                    deflM = .05, txFee = .1, ica = T, numSigs = 2, numBase = 7, numOrig = 3, ...) {

  x <- new_finDemo(sSteps = sSteps, tSteps = tSteps, baseX = baseX, alpha = alpha, gamma = gamma,
                   X = X, deltaX = deltaX, deflM = deflM, txFee = txFee, ...)
  x$ica <- ica
  x$numSigs <- numSigs
  x$numBase <- numBase
  x$numOrig <- numOrig

  x$B1 <- .5*(1+sin(2*pi*1.5*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-.05,.05))+runif(x$tSteps,-.05,.05))
  x$B2 <- .5*(1+cos(2*pi*2.5*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-.05,.05))+runif(x$tSteps,-.05,.05))
  if (x$numOrig > 2) {
    x$B3 <- .5*(1+cos(pi/5+2*pi*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-.05,.05))+runif(x$tSteps,-.05,.05))
  }
  if (x$numOrig > 3) {
    x$B4 <- .5*(1+sin(pi/4+2*pi*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-.05,.05))+runif(x$tSteps,-.05,.05))
  }
  if (is.null(x$X)) {
    if (numOrig == 4) {
      x$X <- tibble(S1 = .5*(x$B1+x$B2), S2 = 1+.5*(x$B2-x$B3), S3 = .7*(x$B3-x$B4), S4 = .5*(x$B4-x$B1))
    } else if (numOrig == 3) {
      x$X <- tibble(S1 = .5*(x$B1+x$B2), S2 = 1+.5*(x$B2-x$B3), S3 = .7*(x$B3-x$B1))
    } else {
      x$X <- tibble(S1 = .5*(x$B1+x$B2), S2 = 1+.5*(x$B1-x$B2))
    }
    for (i in (x$numOrig+1):x$numBase) x$X <-
        bind_cols(x$X, !!paste0('S',i) :=
                    i/(2*5)*x$B1 * (1+runif(x$tSteps,-.01,.01)) +
                    (6-i)/(2*5)*x$B2 * (1+runif(x$tSteps,-.01,.01))
                  )
    x$X <- .95*x$baseX * scale(x$X) # v5:/ ncol(x$X) # X just a scaling factor?
    x$X.orig <- x$X
  }
  #ggplot(as_tibble(x$X) %>% rowid_to_column('t') %>%  gather(key,value,-t))+geom_line(aes(t,value, color=key))

  if (x$ica) {
    x$xpca <- prcomp(x$X.orig)
    x$X <- as.matrix(x$X.orig) %*% x$xpca$rotation[,1:x$numSigs]
    x$xica <- fastICA::fastICA(x$X, n.comp = x$numSigs, tol=.01)
    #x$X <- sweep(x$xica$S,2, apply(x$xica$S,2,min))
    ica.s <- apply(x$xica$S, 2, function(s) { s <- s-min(s); s/max(s) * x$baseX*.98 })
    x$X <- as.data.frame(ica.s)
    names(x$X) <- paste0('S',1:ncol(x$X))
    #x$X <- as_tibble(x$X / (max(x$X)-min(x$X)) * .98*x$baseX) # v5: / ncol(x$X))
  }
  #ggplot(as_tibble(x$X) %>% rowid_to_column('t') %>%  gather(key,value,-t))+geom_line(aes(t,value, color=key))

  x$xColNames <- names(x$X)[str_detect(names(x$X),'S\\d')]
  x$xColNum <- length(x$xColNames)
  x$vColNames <- paste0(x$xColNames,'i')
  x$allColNames <- c(x$xColNames,'S')
  x$itv <- x$sSteps+1

  #valSteps <- baseX * xColNum * c(0,10^((1:sSteps)/sSteps)-1)/10
  x$valSteps <- x$baseX * c(0,(1:x$sSteps)/x$sSteps)
  x$SIntv <- tibble(sid=1:x$itv)
  for (n in x$xColNames) x$SIntv <- bind_cols(x$SIntv, !!n := x$valSteps)
  x$SIntv <- x$SIntv %>% select(-sid)

  x$Si <- tibble()
  for (n in x$xColNames) x$Si <- crossing(x$Si, !!n := x$SIntv[n])
  names(x$Si)<-paste0(names(x$Si),'i')
  x$Si <- rowid_to_column(x$Si, 'sid')

  x$S <- tibble()
  for (n in x$xColNames) x$S <- crossing(x$S, !!n := 1:x$itv)
  x$S <-  x$S %>% rowid_to_column('sid')
  #x$S <- crossing(x$S, sig.p = c(-1,1))
  x$S <- inner_join(x$S, x$Si, 'sid') %>% mutate(n=0)
  #x$S['spid'] <- 1:nrow(x$S)
  x$S <- x$S %>% select(-sid) %>% rowid_to_column('sid')

  x$X$Si <- rowSums(x$X[x$xColNames])
  for (n in x$xColNames) x$X[paste0(n,'.i')] <- findInterval(pull(x$X,n), pull(x$SIntv))
  x$X$S <- findInterval(pull(x$X,'Si'), pull(x$SIntv))
  x$X <- rowid_to_column(x$X, 'xid')

  x$XRev <- x$X %>% arrange(Si)
  #for (n in xColNames) { XRev[paste0(n,'.i')] <- findInterval(pull(XRev,n),pull(SIntv,1)) }

  x$SValues <- rep(0, x$xColNum) # v3: x$X$S[1] + x$SIntv[2,]
  x$SValues.init <- x$SValues

  x$Acts <- tibble()
  for (n in x$xColNames) x$Acts <- crossing(x$Acts, !!n := c(-1, 0, 1))

  x$A <- crossing(x$Acts, sid = x$S$sid) %>% rowid_to_column('aid') %>% mutate(t=0,n=0)

  x$sRows <- nrow(x$S)
  x$aRows <- nrow(x$A)
  #q <- crossing(sid=x$S$sid,aid=x$A$aid) %>% mutate(qid = toQid(sid,aid,x$aRows))
  x$Q <- tibble(aid = x$A$aid)
  for (n in x$xColNames) x$Q[n] <- 0
  #x$Q <- bind_cols(x$Q, sid=q$sid,aid=q$aid)

  x$W <- rep(1, x$xColNum)
  x$Z <- rep(0, x$xColNum)
  x$Qold <- 0

  x$brokeR <- -10000
  x$xCols <- which(names(x$X) %in% x$xColNames)
  x$sCols <- which(names(x$S) %in% x$xColNames)
  x$vCols <- which(names(x$S) %in% paste0(x$xColNames,'i'))
  x$offLimR <- -1
  x$A$dV <- 0
  x$SumR <- .0
  x$numSigs <- numSigs
  x$FS <- mapState(x) %>% mutate(p = Si/x$baseX, dp = c(0,diff(Si))) %>% mutate(sig.p = sign(dp))
  x$R <- 0
  x$t <- 1
  x$trace <- c()
  x$n.t <- microbenchmark::get_nanotime()
  x$logR.cols <- paste0('R.',x$xColNames)
  x$logQpo.cols <- paste0('Qpo.',x$xColNames)
  x$logdQ.cols <- paste0('dQ.',x$xColNames)
  x$logdQn.cols <- paste0('dQn.',x$xColNames)
  x$logQa.cols <- paste0('Qa.',x$xColNames)
  x$logQan.cols <- paste0('Qan.',x$xColNames)
  x$logDelta.cols <- paste0('delta.',x$xColNames)
  x$logW.cols <- paste0('W.',x$xColNames)
  x$logZ.cols <- paste0('Z.',x$xColNames)
  x$log <- tibble(lid=integer(), time=numeric(), n=integer(), t=integer(), sid=integer(), aid=integer(), Qpo=numeric())
  for (n in c(x$logRcols,x$logQcols,x$logdQcols,x$logdQncols,x$logWcols,x$logZcols)) x$log[n] <- numeric()
  x
}

new_finDemo <- function(..., class = character()) {
  new_world(..., class = c(class, 'finDemo')
  )
}

#' initialize environment (e.g. reset Q to zero)
init.finDemo <- function(x, ...) {
  x <- getStart(x)
  x$W <- (x$W * 0) + 1 # /x$baseX/3
  x
}

#' get the starting state of this environment
getStart.finDemo <- function(x, ...) {
  x$t = 1
  x$SValues <- x$SValues.init
  x$SumR <- .0
  x$Z <- x$Z * 0
  x$Qold <- 0
  x$trace <- c()
  x$n.t <- microbenchmark::get_nanotime()
  x
}

#' find the real state for the possibly updated location values (Sx,Sy,..) of state S
#' @param S the possbly updated input state
#' @return the real state element matching the example state
getState.finDemo <- function(x, SValues = x$SValues, ...) {
  #if (deblog) n.t <- nanoTime(paste('getState',x$numSigs), x$n.t)
  #sIdx <- 0
  #for (col in 1:x$xColNum) { sIdx <- (sIdx * x$itv) + findInterval(pull(SValues,col), pull(x$SIntv,col))-1 }

  #if (deblog) x$n.t <- nanoTime(paste('getState',x$numSigs), n.t)
  #return (x$S[x$S$sid == sIdx+1,][1,])
  sid = toSid(t(findInterval(SValues, pull(x$SIntv))),1, x$itv, x$sCols)
  x$S[sid,]
}

mapState.finDemo <- function(x, X = NULL, ...) {
  if (!is.null(X) && is.numeric(X)) {
    si <- x$XRev[findInterval(X,pull(x$XRev,Si)),x$vCols]
    names(si)<-x$xColNames
    si
  } else {
    b <- x$B1+x$B2
    if (!is.null(x$B3)) b <- b + x$B3
    if (!is.null(x$B4)) b <- b + x$B4
    tibble(Si = b * x$baseX/x$numOrig) %>% mutate(S = findInterval(Si, x$valSteps))
  }
}

toQid <- function(sid, aid, aRows) {
  aRows*(sid-1)+aid
}

#' calculate S-index by Si and sig.p
#' validation: sum(x$S$sid-toSid(x$S[,x$xCols],x$S$P))
toSid <- function(s, sig.p, itv, sCols) {
  #m <- as.matrix(cbind(s,p)-1)
  s <- replace(s, s>itv, itv)
  s <- replace(s, s<1, 1)
  m <- as.matrix(s)-1
  m[m<0] <- NA
  m %*% itv^((length(sCols)-1):0) + 1
  #(m %*% itv^((length(sCols)-1):0))*2+ifelse(sig.p<0,0,1) + 1
}

move <- function(x, S, A, t = x$t) {
  if (deblog) n.t <- nanoTime(paste('move',x$numSigs), x$n.t)
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$A[A,]
  s.idx <- sweep(A[x$xColNames],2,-as.numeric(S[x$sCols]))
  #s.p <- right_join(x$S[x$S$P == as.numeric(x$FS$S[t]),], s.idx, by=x$xColNames)
  #s.p <- x$S[toSid(s.idx, as.numeric(x$FS$S[t]), x$itv, x$sCols),] # P all intervals
  s.p <- x$S[toSid(s.idx, x$FS$sig.p[t], x$itv, x$sCols),]
  if (deblog) print(paste0('move[',x$i,']: S=',toString(S),' A=',toString(list(A[1,])),', next=',toString(s.p$sid)))
  if (deblog) x$n.t <- nanoTime(paste('move',x$numSigs), n.t)
  s.p
}

estimateR <- function (x, S, S.next = NULL, t = x$t) {
  if (deblog) n.t <- nanoTime(paste('estimateR',x$numSigs), x$n.t)

  #p <- x$FS$p[t]
  if (is.null(S.next)) {
    #return(rowSums(S[,x$vCols]) * x$FS$p[t])
    return(-(abs(rowSums(S[,x$vCols]) - x$FS$Si[t])))
  }

  xn <- x$vCols
  #p.t <- ifelse(t > 1, x$FS$p[t-1], p)
  #dp <- p - p.t
  #dS.tp <- sweep(S.next[,xn], 2, as.numeric(S[xn])) * p.t
  #s <- S.next[,xn] * dp - dS.tp
  #s <- rowSums(s)
  s <- -abs(S.next[,xn] - x$FS$Si[t+1])
  if (deblog) x$n.t <- nanoTime(paste('estimateR',x$numSigs), n.t)
  s
}

#' R: reward of the State/Action pair
#' t: time t in State.next (after action A)
#' Difference of QuasiOpt policy (F(t) * delta.price - delta.F * price(t-1))
#' not weighted by W?
#' and the actual value of state: sum(S.next * delta.price - delta.S * price(t-1))
calculateR <- function(x, S, S.next = NULL, t = x$t) {
  if (deblog) n.t <- nanoTime(paste('calculateR',x$numSigs), x$n.t)

  if (is.null(S.next)) {
    #return(rowSums(S[,x$vCols] * as.numeric(x$W)) * x$FS[t,]$p - x$FS$Si[t])
    #return(rowSums(S[,x$vCols]) * x$FS[t,]$p)
    #return(S[,x$vCols] * x$FS[t,]$p)
    return(-abs(S[,x$vCols] - x$FS$Si[t]))
  }

  xn <- x$vCols
  #p.t <- x$FS$p[t-1]
  #dp <- x$FS$p[t] - p.t
  #v4: the to-be-expected target
  #target.s <- x$FS$Si[t]*dp - (x$FS$Si[t]-x$FS$Si[t-1])*x$FS$p[t-1]

  #dS.tp <- sweep(S.next[,xn], 2, as.numeric(S[xn])) * p.t
  #s <- S.next[,xn] * dp - dS.tp
  #s <- t(apply(s, 1, `*`, as.numeric(x$W)))
  #s <- rowSums(s) - target.s
  #v4: s <- s - target.s
  s <- -abs(S.next[,xn]  - x$FS$Si[t])
  if (deblog) x$n.t <- nanoTime(paste('calculateR',x$numSigs), n.t)
  s
}

calculateDx <- function(x, S, A, S.moves, alpha = x$alpha, t = x$t) {
  if (deblog) n.t <- nanoTime(paste('calculateDx',x$numSigs), x$n.t)
  R.estimate <- estimateR(x, S, S.moves, t)
  R.estimate[is.na(R.estimate)] <- -x$baseX
  #dX <- S.moves[x$vCols]/x$xColNum * R.estimate
  Q <- x$Q[A$aid,x$xCols]
  Q[is.na(Q)] <-  0
  #dx <- (dX/(1+rowSums(abs(dX))) * alpha / (1+rowSums(abs(Q)))) + Q
  #dx <- dX + Q
  dx <- R.estimate
  if (deblog) { print(paste0('calculateDx[',x$i,']:')); print(Q); print(dx) }
  if (deblog) x$n.t <- nanoTime(paste('calculateDx',x$numSigs), n.t)
  dx
}

#' get the real State/Action for the possibly updated S and A
#' @param state S as tibble/list
#' @param action A as tibble/list or factor
#' @return a merged object of S and A
#' @seealso getState
getStateAction.finDemo <- function(x, S, A, ...) {
  if (is.null(x) || is.null(S) || is.null(A)) stop('missing input')
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$A[A,]
  A$t <- x$t
  bind_cols(S, A)
}

#' get possible actions (tibble) for the state S
#' @param state S as tibble/list or numeric (sid)
#' @return tibble of actions
#'
possibleActions.finDemo <- function(x, S, SValues = x$SValues) {
  if (deblog) n.t <- nanoTime(paste('possible',x$numSigs), x$n.t)
  if (is.null(x) || is.null(S)) stop('missing input')
  if (is.numeric(S)) {
    S <- x$S[S,]
  }

  #xA <- x$A[x$A$sid == S$sid,] %>% filter(tx == x$t) %>% mutate(z=ifelse(n==0,1,n)/sum(n))
  xA <- x$A[x$A$sid == S$sid,]

  #dX <- x$X[x$t+1,] - x$X[x$t,] # calculate Delta-X
  S.moves <- move(x, S, xA)
  xA <- xA[!is.na(S.moves$sid),]
  S.moves <- S.moves[!is.na(S.moves$sid),]
  #xA$spid <- S.moves$spid

  dX <- calculateDx(x, S, xA, S.moves)
  xA <- bind_cols(xA,dX)

  W <- x$W # rep(1,x$xColNum) # x$W
  # remember: apply-transpsed-rowSums!
  Q <- colSums(apply(dX,1,`*`,W)) # the predicted cost/value of a move
  xA$Q <- Q

  #if (any(abs(Q)>x$baseX*10)) browser()
  if (nlog) {
    print(paste0('possible[',x$i,']: t=',x$t,', S=',toString(list(S)),', Q=',toString(Q)))
  }
  #if (deblog) print(paste0('possible[',x$i,']: t=',x$t,', W=',toString(W)))

  # by: Value increase by Delta-X, sell: Value decreased by Delta-X:
  po <- xA %>% arrange(desc(Q))

  if (deblog) print(tibble(m=paste0('possible[',x$i,']:'), po), width=Inf)
  if (deblog) x$n.t <- nanoTime(paste('possible',x$numSigs), n.t)

  po
}

#' get the index for the currently worst performing share to sell
#' @param action A as tibble/numeric
getReplacement <- function(x, A) {
  if (is.numeric(A)) A <- x$A[A,]
  S <- x$S[A$sid,]
  v <- S[x$vCols] # current volumes
  dX <- enframe(x$X[x$t+1,] - x$X[x$t,], name='Xx', value='dX') # current value changes
  va <- gather(v,'Sx','v') %>% rowid_to_column('idx') %>% bind_cols(dX) %>%
    filter(idx!=A$A & Sx > 0) %>% mutate(dv = v*dX) %>% arrange(dv) %>%  head(1)
  if (nlog) print(paste0('getReplacement[',x$i,']: ', toString(list(va))))
  va
}

#' Updates state S according action A (as tibble/list or as factor/character)
#' @param state S as tibble/list
#' @param action A as tibble/list
#' @return updated real (S,A,SValues) object
evalAction.finDemo <- function(x, S, A, SValues = x$SValues, ...) {
  if (deblog) n.t <- nanoTime(paste('evalAction',x$numSigs), x$n.t)
  if (is.null(x) || is.null(S) || is.null(A)) stop('missing input')
  if (is.factor(A)||is.character(A)) { print(A); stop('A cannot be primitive') }
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$A[A,]

  AIn <- A
  dMFee = 0

  R <- rep(0,x$xColNum)
  SValues <- S[x$vCols]

  if (all(A[x$xCols] == 0)) { # no action
    #dX <- x$X[x$t+1,] - x$X[x$t,] # calculate Delta-X
    #A$R <- sum((1 - x$X[,x$t]) * S[x$vCols])
    R <- calculateR(x, S, S, t=x$t+1)
  } else { # try to trade
    #vRepl <- getReplacement(x,AIn)
    S.move <- move(x, S, A)
    if (is.null(S.move) || nrow(S.move) == 0 || is.na(S.move$sid)) {
      A <- NULL
      R <- calculateR(x, S, S, t=x$t+1)
    } else {
      S.next <- S.move
      R <- calculateR(x,S.next, t=x$t+1)
      S <- S.next
      SValues <- S[x$vCols]
    }
  }

  if (!is.null(A)) {
    print(paste0('ACTION[',x$i,']: t=',x$t,', A=',toString(list(A))))
    print(paste0('ACTION[',x$i,']: t=',x$t,', S=',toString(list(S)),', SValues=',toString(SValues),'/',toString(x$SValues)))
  } else {
    if (nlog) print(paste0('DENIED[',x$i,']: t=',x$t,', S=',toString(list(S))))
    if (nlog) print(paste0('DENIED[',x$i,']: t=',x$t,', A=',toString(list(AIn)),', SValues=',toString(SValues),'/',toString(x$SValues)))
  }
  if (deblog) x$n.t <- nanoTime(paste('evalAction',x$numSigs), n.t)
  return (list(S = S, A = A, R = R, SValues = SValues))
}

#' act out the action A in state S
#' @param state S as tibble/list
#' @param action A as tibble/list or factor/character
#' @return the new state S' augmented by reward R
#' @seealso evalAction
#' @seealso getState
act.finDemo <- function(x, S, A, SValues = x$SValues, ...) {
  if (deblog) n.t <- nanoTime(paste('act',x$numSigs), x$n.t)
  if (is.null(x) || is.null(S) || is.null(A)) stop('missing input')
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$A[A,]

  terminal <- F
  x$resA <- NULL
  SIn <- S
  SValIn <- SValues
  x$trace[x$t] <- SIn$sid

  res <- evalAction(x, S, A, SValIn)
  S <- res$S
  SValues <- res$SValues

  if (F&&sum(SValues) <= 1) { # BROKE!! end of episode
    print(paste0('Broke[',x$i,']! t=',x$t,', S=',toString(list(S)),', SValues=',toString(SValues)));
    R <- rep(-1000, x$xColNum)
    terminal = T
  } else {
    R <- res$R
    if (!is.null(res$A)) {
      A <- res$A
      A$sid <- res$S$sid
    }
  }

  x$t <- x$t + 1
  if (x$t >= x$tSteps) { # end of episode!
    terminal = T
    # v2: R <- calculateR(x, SIn) + x$SumR
    R <- x$SumR
    R.r <- R
  } else {
    x$SumR <- R + x$SumR
    # absolute delayed reward
    R.r <- R
    R <- R * 0
    if (sum(SValues) == 0) { # penalize no investment
      R <- rep(-x$baseX, x$xColNum)
      A$failed <- T
      #R.r <- R -- remember real R
    }
  }

  if (terminal) {
    print(paste0('END[',x$i,']: t=',x$t,', R=',toString(R),', SValues=',toString(x$SValues),', S=',toString(list(S))))
  }

  A$in.sid <- SIn$sid
  A$sid <- S$sid
  #A$spid <- S$spid
  A$terminal <- terminal
  A$t <- x$t
  A$n <- A$n+1

  S$n <- S$n+1
  x$S$n[S$sid] <- S$n

  x$A$t[A$aid] <- x$t
  #x$A$spid[A$aid] <- S$spid
  x$A$n[A$aid] <- A$n
  x$resA <- list(A=A, R=R, SValues = SValues, R.r = R.r)
  if (!terminal & nlog) {
    print(paste0('act[',x$i,']: t=',x$t,', resA=',toString(list(x$resA)),', SIn=',toString(list(SIn))));
    print(paste0('act[',x$i,']: t=',x$t,', R=',toString(R.r),', SumR=',toString(x$SumR),', SValues=',toString(x$resA$SValues),', S=',toString(list(S))))
  }
  if (deblog) x$n.t <- nanoTime(paste('act',x$numSigs), n.t)
  return (x)
}

#' check if acting out action A in state S will leave the world's context
#' @param state S as tibble/list
#' @param action A as tibble/list or factor/character
#' @return true if next state would be outside world's context
#' @seealso evalAction
offLimits.finDemo <- function(x, S, A, ...) {
  if (is.null(x) || is.null(S) || is.null(A)) stop('missing input')
  S <- evalAction(x, S, A)
  is.null(S) || nrow(S) == 0
}

#' update W for this action A
#' @param action A as tibble/list (match is done by A$aid) or as aid
#' @param return R from last action
#' @return the environment
updateQ.finDemo <- function(x, A, R, Anext = NULL, alpha = x$alpha, gamma = x$gamma, lambda = .8, ...) {
  if (deblog) n.t <- nanoTime(paste('updateQ',x$numSigs), x$n.t)
  if (deblog) {
    print(paste0('updateQ[',x$i,']: t=',x$t,', A=',toString(list(A)),', R=',toString(R)));
    if(!is.null(Anext)) print(paste0('updateQ[',x$i,']: t=',x$t,', Anext=',toString(Anext)));
  }
  if (is.null(x) || is.null(A) || is.null(R)) stop('missing input')
  if (is.numeric(A)) A <- x$A[A,]
  if (is.numeric(Anext)) Anext <- x$A[Anext,]

  aid <- A$aid
  S <- x$S[A$sid,]
  Rc <- R #* replace(A[x$xCols],A[x$xCols]==0,1) / 2

  # Update Q
  currQ <- x$Q[A$aid,x$xCols]
  currQ.org <- currQ
  if (!is.null(Anext)) {
    currQ <- currQ + alpha * (Rc + gamma * x$Q[Anext$aid,x$xCols] - currQ)
  } else {
    currQ <- currQ + alpha * (Rc - currQ)
  }
  if (nlog) print(paste0('updateQ[',x$i,']: x$Q=',toString(as.numeric(currQ.org)),', currQ=',toString(currQ)))
  x$Q[A$aid,x$xCols] <- currQ

#  if (is.null(Anext) && x$t>1) {
#    q <- x$Q[A$sid,x$xCols]
#    for (t in (x$t-1):1) {
#      q <- x$gamma * q
#      x$Q[x$trace[t],x$xCols] <- x$Q[x$trace[t],x$xCols] + q
#    }
#  }

  # Update W
  W <- x$W # [A$aid,]
  #dQ <- currQ # / x$baseX
  dQ <- estimateR(x,S)
  Qa <- W * dQ
  Qa.next <- NA
  dQ.next <- NA
  delta <- NA

  if (is.null(x$resA) || !has_name(x$resA$A,'failed') || !x$resA$A$failed) {
    if (is.null(Anext)) {
      dW <- (sum(Rc - Qa)) * dQ
      W <- W + alpha * dW # / max(1,abs(as.numeric(dW)))
      if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dQ=',toString(dQ),', Qa=',toString(Qa),', dW=',toString(dW)))
    } else {
      dQ.next <- calculateDx(x, S, Anext, x$S[Anext$sid,]) # / x$baseX
      Qa.next <- W * dQ.next # dQ.next # predicted next day value

      delta <- Rc + gamma * (Qa.next - Qa)

      dQ.z <- dQ # / max(1,abs(as.numeric(dQ)))
      x$Z <- lambda * x$Z + (1 - alpha*gamma*lambda*sum(x$Z*dQ.z)) * dQ.z
      #x$Z <- x$Z / sum(.001+abs(x$Z))

      dw.dqz <- as.numeric(sum(delta + Qa - x$Qold) * x$Z)
      dw.qz <- as.numeric(sum(Qa - x$Qold) * dQ.z)
      W <- W + alpha * dw.dqz/max(1,abs(dw.dqz))  - alpha * dw.qz #/max(1,abs(dw.qz)) 

      if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dQ=',toString(dQ),', Qa=',toString(Qa),', delta=',toString(delta),', dQ.z=',toString(dQ.z)))
      if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dQ.n=',toString(dQ.next),', Qa.n=',toString(Qa.next),', Z=',toString(x$Z),', Qold=',toString(x$Qold)))
      if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dw.dqz=',toString(dw.dqz),'/',toString(dw.dqz/max(1,abs(dw.dqz))),', dw.qz=',toString(dw.qz),'/',toString(dw.qz/max(1,abs(dw.qz)))))
      x$Qold <- Qa.next
    }

    if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', x$W=',toString(x$W),', W=',toString(W)))
    #W <- W / max(replace(abs(as.numeric(W)), 0, 1))
  } else {
    if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,' failed: x$W=',toString(x$W),', W=',toString(W)))
  }

  #if (any(abs(W)>x$baseX)) browser()
  
  po <- x$chosenAction
  lid <- nrow(x$log)+1
  x$log <- add_row(x$log, lid = lid, time = microbenchmark::get_nanotime())
  x$log$n[lid] <- x$i
  x$log$t[lid] <- x$t
  x$log$sid[lid] <- A$sid
  x$log$aid[lid] <- A$aid
  x$log[lid,x$logR.cols] <- x$resA$R.r
  x$log[lid,x$logdQ.cols] <- dQ
  x$log[lid,x$logdQn.cols] <- dQ.next
  x$log[lid,x$logQa.cols] <- Qa
  x$log[lid,x$logQan.cols] <- Qa.next
  x$log[lid,x$logDelta.cols] <- delta
  x$log[lid,x$logW.cols] <- W
  x$log[lid,x$logZ.cols] <- x$Z
  x$log$Qpo <- po$Q
  x$log[lid,x$logQpo.cols] <- po[x$vColNames]

  if (nlog) {
    print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', A=',toString(list(A))))
    print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', Rc=',toString(Rc),', Qa=',toString(Qa),', Qa.next=',toString(Qa.next)))
  }

  wd <- as.numeric(W)
  if (any(is.na(wd) | is.nan(wd) | wd == Inf | wd == -Inf)) {
    print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', A=',toString(A),', x$W=',toString(x$W)));
    print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', R=',R,', Rc=',toString(Rc)));
    print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', Qa=',toString(Qa),', W=',toString(W)));
    stop(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,' invalid W: ', toString(W)))
  }
  if (any(x$SValues < 0)) {
    print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', SValues:',toString(SValues))); stop('Invalid SValues: going SHORT!')
  }

  x$W <- as.numeric(W)
  x$A$t[aid] <- x$t
  x$A$n[aid] <- x$A$n[aid] + 1
  if (deblog) x$n.t <- nanoTime(paste('updateQ',x$numSigs), n.t)
  return (x)
}

nanoTime <- function(action = 'undef', n.t.start = n.t) {
  n.t.next <- microbenchmark::get_nanotime()
  print(paste(action,'time:',(n.t.next - n.t.start)/10^6))
  n.t.next
}

#' test core functions of this world
test.finDemo <- function(x, ...) {
  deblog.orig <- deblog
  deblog <- F
  nlog.orig <- nlog
  nlog <- F
  numSigs <- 2
  envfinDemo <- finDemo(numSigs = numSigs)

  #assertthat::are_equal(envfinDemo$S[200,]$sid,
  #                      as.numeric(convertBase(as.numeric(envfinDemo$S[200,1+1:envfinDemo$xColNum]),envfinDemo$sSteps+1)))
  #assertthat::are_equal(envfinDemo$S[200,]$sid, indexS(envfinDemo, envfinDemo$S[200,]))
  #assertthat::are_equal(envfinDemo$S$sid, indexS(envfinDemo, envfinDemo$S))
  #microbenchmark::microbenchmark(indexS(envfinDemo, envfinDemo$S[1,]),times=100)
  #microbenchmark::microbenchmark(indexS(envfinDemo, envfinDemo$S),times=100)

  print(sloop::s3_dispatch(possibleActions(envfinDemo, tibble())))
  startEnv <- getStart(envfinDemo)
  assertthat::assert_that(startEnv$t == 1 && startEnv$SValues[1] > 0)
  startState <- getState(envfinDemo, startEnv$SValues)
  assertthat::assert_that(startState$S1 > 1 && startState[2] > 1)
  nextSValues <- startEnv$SValues
  nextSValues[1] <- startEnv$baseX/2
  upState <- getState(envfinDemo, nextSValues)
  assertthat::assert_that(upState$S1 < startState$S1)
  microbenchmark::microbenchmark(getState(envfinDemo, startEnv$SValues))

  poss <- possibleActions(envfinDemo, startState)
  print(poss)
  assertthat::assert_that(identical(poss, possibleActions(envfinDemo, startState$sid)))
  #assertthat::assert_that(nrow(poss) == length(envfinDemo$A))
  microbenchmark::microbenchmark(possibleActions(envfinDemo, startState),times=10)

  envfinDemo <- finDemo(numSigs = numSigs)
  myAct <- poss %>% tail(1)
  updState <- evalAction(envfinDemo, myAct$sid, myAct$aid)
  assertthat::assert_that(updState$SValues[1] > 0)
  updState <- evalAction(envfinDemo, myAct$sid, myAct$aid)
  assertthat::assert_that(updState$SValues[1] > 0)
  microbenchmark::microbenchmark(evalAction(envfinDemo, myAct$sid, myAct$aid))

  envfinDemo <- finDemo(numSigs = numSigs)
  myAct <- poss %>% tail(1)
  actEnv <- act(envfinDemo, myAct$sid, myAct$aid)
  print(actEnv$resA)
  assertthat::assert_that(actEnv$resA$R != 0 && actEnv$t == 2)
  actEnv2 <- act(actEnv, actEnv$resA$sid, actEnv$resA$aid)
  print(actEnv2$resA)
  assertthat::assert_that(actEnv2$resA$R != 0 && actEnv2$t == 3)
  microbenchmark::microbenchmark(act(actEnv, myAct$sid, myAct))

  envfinDemo <- finDemo(numSigs = numSigs)
  myAct <- max(poss$aid)
  deadActEnv <- act(envfinDemo, envfinDemo$S[1,], myAct)
  assertthat::assert_that(deadActEnv$resA$terminal)

  envfinDemo <- finDemo(numSigs = numSigs)
  startState <- getState(envfinDemo, envfinDemo$SValues)
  SA <- getStateAction(envfinDemo, startState, envfinDemo$A[1,])
  assertthat::assert_that(SA$sid == startState$sid && SA$aid == envfinDemo$A[1,]$aid)
  SA <- getStateAction(envfinDemo, startState$sid, SA$aid)
  assertthat::assert_that(SA$sid == startState$sid && SA$aid == envfinDemo$A[SA$aid,]$aid)
  microbenchmark::microbenchmark(getStateAction(envfinDemo, startState, envfinDemo$A[1,]))

  envfinDemo <- finDemo(numSigs = numSigs)
  startEnv <- getStart(envfinDemo)
  startState <- getState(startEnv, startEnv$SValues)
  myAct <- envfinDemo$A %>% filter(sid == startState$sid) %>% head(1)
  nextAct <- envfinDemo$A %>% filter(sid == (move(envfinDemo, startState, myAct)$sid)) %>% head(1)
  env2 <- updateQ(envfinDemo, myAct$aid, .33, nextAct$aid)
  assertthat::assert_that(all(env2$W != 1))
  env3 <- updateQ(env2, myAct$aid, 44)
  assertthat::assert_that(!identical(env3$W, env2$W))
  microbenchmark::microbenchmark(updateQ(envfinDemo, SA$aid, 44))

  deblog <- deblog.orig
  nlog <- nlog.orig
}

######################### The Cliff ##########################
##############################################################