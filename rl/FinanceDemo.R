##############################################################
source('worlds.R')
source('utils.R')
######################### Finance Demo #######################
sSteps = 20; tSteps = 20; baseX = 100; X = NULL; deltaX = .1;
ica=T; numSigs=2; numBase = 10; numOrig = 4
fErr.v0 = baseX; fErr.base = baseX; fErr.a = 10; fs.pred='S'
ica.ratio = .5; ica.maxRetries = 500; ica.high = 10

finDemo <- function(sSteps = 10, tSteps = 20, baseX = 100,
                    epsilon = .1, alpha = .1, gamma = .9, X = NULL, deltaX = .1,
                    ica = T, numSigs = 2, numBase = 10, numOrig = 3, ica.ratio = 1e-3, ica.maxRetries = 500, ica.high = 10,
                    fErr.v0 = baseX, fErr.base = baseX, fErr.a = 10, fs.pred = 'S', ...) {

  x <- new_finDemo(sSteps = sSteps, tSteps = tSteps, baseX = baseX,
                   epsilon = epsilon, alpha = alpha, gamma = gamma,
                   X = X, deltaX = deltaX, ...)
  x$run <- 1
  x$numSigs <- numSigs
  x$ica <- ica
  x$numBase <- numBase
  x$numOrig <- numOrig
  x$fErr.v0 <- fErr.v0
  x$fErr.base <- fErr.base
  x$fErr.a <- fErr.a
  x$fs.pred <- fs.pred
  x
}

initIntern <- function(x, ...) {
  eps.t <- .005
  eps.x <- .01

  if (T || is.null(x$X)) { # allways rerun ICA

    x$B1 <- .5*(1+sin(2*pi*1.53*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
    x$B2 <- .5*(1+cos(2*pi*1.92*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
    if (x$numOrig > 2) {
      x$B3 <- .3*(1+cos((pi/5+2*pi*2.34)*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
    }
    if (x$numOrig > 3) {
      x$B4 <- .4*(1+sin((pi/4+2*pi*2.97)*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
    }
    if (x$numOrig > 4) {
      x$B5 <- .5*(1+sin((pi/5+2*pi*3.45)*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
    }
    if (x$numOrig > 5) {
      x$B6 <- .6*(1+cos((pi/4+2*pi*3.97)*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
    }
    x$b <- tibble(.rows = length(x$B1))
    for (name in names(x)[grep('B\\d',names(x))]) x$b <- bind_cols(x$b, x[name])

    kx = .2
    x$X <- tibble(.rows = length(x$B1))
    for (i in 1:x$numOrig) {
      idx <- 1+(i-1)%%x$numOrig
      idx2 <- 1+(i)%%x$numOrig
      x[[paste0('S',i)]] <- x[[paste0('B',idx)]] + kx * x[[paste0('B',idx2)]]
      x$X <- bind_cols(x$X, x[paste0('S',i)])
    }

    # Add filler columns to support NumSig signal evaluation
    ki <- seq(.1,2,.1)
    for (i in (x$numOrig+1):(x$numBase)) {
      idx <- 1+(i%%x$numOrig)
      idx2 <- 1+((i+1)%%x$numOrig)
      xs <- -.01 * (ki[i]*(-1)^idx * x[[paste0('B',idx)]] + ki[i+1]*(-1)^idx2 * x[[paste0('B',idx2)]])
      xs <- xs - min(xs)
      x$X <- mutate(x$X, !!paste0('S',i) := xs)
    }

    x$X <- scale(x$X, center = T, scale = F)
    x$X <- x$X * .8*x$baseX / (max(x$X)-min(x$X)) # v5:/ ncol(x$X) # X just a scaling factor?
    x$X.orig <- as_tibble(x$X)
    x$X.orig <- x$X.orig %>% mutate_all(list(~ . + .001*runif(length(.))))
    print(paste('Signals generated:',object.size(x)))

    if (x$ica) {

      x$ica.Xlog = tibble(i=integer(),ratio=numeric(), ratio.b=numeric(),ratio.xo=numeric(),ratio.x=numeric())
      ica.retries = 0
      ica.bestX <- NULL
      ica.bestRatio <- 1000000
      ica.worstRatio <- 0
      print('Runnig ICA on X.orig:')
      print(x$X.orig)

      repeat {

        x$xpca <- prcomp(x$X.orig)
        x$X.pca <- as.matrix(x$X.orig) %*% x$xpca$rotation[,1:min(x$numSigs,dim(x$xpca$rotation)[2])]
        tryCatch({
          x$xica <- fastICA::fastICA(x$X.orig, n.comp = x$numSigs, maxit=100, tol = .01, verbose = T)
        }, error = function(e) { print(paste('cannot execute ICA:',e)); print(x$X.orig); stop('cannot execute ICA:',e) })
        print('ICA calculated..')
        # Rescaling to positive values:
        #   - move restricting to positive gradient -> would move to negative values
        x$X <- as_tibble(x$xica$S, .name_repair = ~paste0('S',1:ncol(x$xica$S))) %>%
          transmute_all(list(~./max(abs(c(max(.),min(.)))) * x$baseX*.95)) # negative States

        freqRes <- eval_freq2(x, high = ica.high)

        ## No useful Freq-Analysis:
        if (is.na(freqRes$ratio)) {
          ica.retries <- ica.retries+1
          if (ica.retries > ica.maxRetries) break
          print(paste('Ratio is NA:',toString(freqRes)))
          next
        }

        x$ica.Xlog <- add_row(x$ica.Xlog, i = ica.retries, ratio = freqRes$ratio,
                              ratio.b=freqRes$ratio.b, ratio.xo=freqRes$ratio.xo, ratio.x=freqRes$ratio.x)

        if (freqRes$ratio > ica.worstRatio) {
          ica.worstRatio <- freqRes$ratio
          x$ica.worstX <- x$X
        }
        if (freqRes$ratio <= ica.ratio) {
          break
        }
        if (freqRes$ratio < ica.bestRatio) {
          ica.bestRatio <- freqRes$ratio
          ica.bestX <- x$X
        }
        if (ica.retries > ica.maxRetries) {
          x$X <- ica.bestX
          print(ica.bestX)
          print(paste('ICA BestRatio: ', ica.bestRatio))
          break
        } else {
          print(paste('retrying ICA: ratio',freqRes$ratio,'best',ica.bestRatio,'retry',ica.retries))
        }
        ica.retries <- ica.retries+1
      }
    }
  }

  #ggplot(as_tibble(x$X) %>% rowid_to_column('t') %>%  gather(key,value,-t))+geom_line(aes(t,value, color=key))
  #print(x$X)

  x$X$SC <- x$baseX*.95
  x$xColNames <- names(x$X)[str_detect(names(x$X),'^S(\\d|C)$')]
  x$xColNum <- length(x$xColNames)
  x$xSigColNames <- x$xColNames[str_detect(x$xColNames,'^S\\d$')]
  x$xSigColNum <- length(x$xSigColNames)
  x$vColNames <- paste0(x$xColNames,'i')
  x$allColNames <- c(x$xColNames,'S')
  x$itv <- x$sSteps+1

  ####### Interval definitions
  x$valSteps <- x$baseX * c(0,1:(x$sSteps))/x$sSteps
  x$SIntv <- tibble(sid=1:x$itv)
  for (n in x$xColNames) x$SIntv <- bind_cols(x$SIntv, !!n := x$valSteps)
  x$SIntv <- x$SIntv %>% select(-sid)

  ####### X preparation
  x$X$Si <- rowSums(x$X[x$xColNames])
  for (n in x$xColNames) x$X[paste0(n,'.i')] <- findInterval(pull(x$X,n), pull(x$SIntv))
  x$X$S <- findInterval(pull(x$X,'Si'), pull(x$SIntv))
  x$X <- rowid_to_column(x$X, 'xid')
  x$xCols <- which(names(x$X) %in% x$xColNames)
  x$xSigCols <- which(names(x$X) %in% x$xSigColNames)

  ###### Signal preparation
  x$Si <- tibble(sid = 1)
  for (n in x$xSigColNames) x$Si <- mutate(x$Si, !!n := 0)
  x$Si <- x$Si %>% select(-sid)
  names(x$Si)<-paste0(names(x$Si),'i')
  x$Si <- x$Si %>% mutate(SCi = 0) %>% rowid_to_column('sid')

  x$S <- tibble(sid = 0)
  for (n in x$xSigColNames) x$S <- mutate(x$S, !!n := 1)
  x$S <-  x$S %>% mutate(SC = 0) #%>% rowid_to_column('sid')
  x$S <- bind_cols(x$S, x$Si %>% select(-sid)) %>% mutate(n=0)
  x$S <- x$S %>% select(-sid) %>% rowid_to_column('sid')
  x$sCols <- which(names(x$S) %in% x$xColNames)
  x$vCols <- which(names(x$S) %in% paste0(x$xColNames,'i'))
  x$vSigCols <- which(names(x$S) %in% paste0(x$xSigColNames,'i'))
  x$maxSid <- x$itv^x$xSigColNum
  print(paste('States generated:',object.size(x)))

  x$FS <- mapState(x) %>% mutate(p = Si/x$baseX, dp = c(0,diff(Si))) %>% mutate(sig.p = sign(dp))

  x$XRev <- x$X %>% arrange(Si)

  sv <- x$X[1,x$xCols]
  sv$SC <- 0
  sv$SC <- x$baseX/x$numOrig
  # correlation of signals:
  xb <- bind_cols(x$X %>% select(matches('^S\\d$')), FS=x$FS$Si)
  x$W.init <- c((corrr::correlate(xb))$FS[1:x$xSigColNum],1/x$numSigs)
  x$W.init <- x$W.init/max(abs(x$W.init))/x$numSigs
  x$W <- x$W.init
  x$W.prev <- x$W
  x$SValues.init <- abs(sv)
  x$SValues <- x$SValues.init
  x$SValues.prev <- x$SValues

  acts <- c(-1, 0, 1)
  x$Acts <- tibble(S1=acts)
  for (n in x$xSigColNames[-1]) x$Acts <- crossing(x$Acts, !!n := acts)
  x$Acts$SC <- 0

  x$A <- crossing(x$Acts, sid = x$S$sid) %>% rowid_to_column('aid') %>% mutate(t=0,n=0)
  print(paste('Actions generated:',object.size(x)))

  x$sRows <- nrow(x$S)
  x$aRows <- nrow(x$A)

  x$Z <- rep(0, x$xColNum)
  x$Qold <- 0

  x$brokeR <- -10000
  x$offLimR <- -1
  x$A$dV <- 0
  x$SumR <- .0
  x$numSigs <- numSigs
  x$R <- 0
  if (!has_name(x, 't')) {
    x$t <- 1
    x$i <- 0
  }
  x$currEpsilon <- getEpsilon(x)
  if (!has_name(x, 'trace')) {
    x$trace <- c()
  }
  x$n.t <- microbenchmark::get_nanotime()
  if (!has_name(x, 'log')) {
    x$logSVals.cols <- paste0('SVal.',x$xColNames)
    x$logR.orig.cols <- paste0('R.orig.',x$xColNames)
    x$logR.cols <- paste0('R.',x$xColNames)
    x$logSumR.cols <- paste0('SumR.',x$xColNames)
    x$logQpo.cols <- paste0('Qpo.',x$xColNames)
    x$logdQ.cols <- paste0('dQ.',x$xColNames)
    x$logdQn.cols <- paste0('dQn.',x$xColNames)
    x$logQa.cols <- paste0('Qa.',x$xColNames)
    x$logQan.cols <- paste0('Qan.',x$xColNames)
    x$logDelta.cols <- paste0('delta.',x$xColNames)
    x$logW.cols <- paste0('W.',x$xColNames)
    x$logZ.cols <- paste0('Z.',x$xColNames)
    x$log <- tibble(run = integer(), lid=integer(), time=numeric(), i=integer(), t=integer(), sid=integer(), aid=integer(), Qpo=numeric(), epsilon=numeric())
    for (n in c(x$logSVals.cols,x$logR.cols,x$logR.orig.cols,x$logSumR.cols,x$logQpo.cols,x$logdQ.cols,x$logdQn.cols,
                x$logQa.cols,x$logQan.cols,x$logDelta.cols,x$logW.cols,x$logZ.cols)) x$log[n] <- numeric()
  } else {
    x$log <- head(x$log, 0)
  }
  print(paste('Environment generated:',object.size(x)))
  x
}

new_finDemo <- function(..., class = character()) {
  new_world(..., class = c(class, 'finDemo')
  )
}

#' initialize environment (e.g. reset Q to zero)
init.finDemo <- function(x, ...) {
  x <- initIntern(x, ...)
  x$W <- x$W.init
  x$W.prev <- x$W
  x <- getStart(x)
  x
}

#' get the starting state of this environment
getStart.finDemo <- function(x, ...) {
  x$t = 1
  x$currEpsilon <- getEpsilon(x)
  x$SValues.init <- updateInitialState(x)
  x$SValues <- x$SValues.init
  x$SValues.prev <- x$SValues
  x$SumR <- .0
  x$Z <- x$Z * 0
  x$Qold <- 0
  x$trace <- c()
  x$n.t <- microbenchmark::get_nanotime()
  print(paste0('updateInitial: SValues=',toString(x$SValues),', init=',toString(x$SValues.init)))
  x
}

updateInitialState <- function(x) {
  if (deblog) print(paste0('updateInitial: SValues.init=',toString(x$SValues.init)))
  if (nrow(x$log) > 1000) {
    S <- getState(x, x$SValues.init)
    if (nrow(S) == 0) {
      stop(paste0('updateInitial: bad InitialState: init=',toString(x$SValues.init)))
    }

    sv.mean <- x$log %>% filter(t==2) %>% select(matches('SVal')) %>% tail(100) %>% summarize_all(mean)
    sv.diff <- x$alpha * (sv.mean - x$SValues.init)
    svinit <- x$SValues.init + sv.diff

    S <- getState(x, svinit)
    if (nrow(S)>0 && !any(is.na(S$sid))) {
      print(paste0('updateInitial: diff=',toString(sv.diff),', W=',toString(x$W)))
      return(svinit)
    } else  {
      print(paste0('updateInitial: no State: SValues=',toString(svinit),', init=',toString(x$SValues.init),', m=',toString(m$SValues),', topA=',toString(topA)))
    }
  }
  x$SValues.init
}

#' find the real state for the possibly updated location values (Sx,Sy,..) of state S
#' @param S the possbly updated input state
#' @return the real state element matching the example state
#' @example sum(x$S$sid - (getState(x, x$S[,x$vCols]))$sid) == 0
getState.finDemo <- function(x, SValues = x$SValues, ...) {
  SvSigs <- as.matrix(SValues[,1:x$xSigColNum])
  s <- matrix(findInterval(SvSigs, pull(x$SIntv)), nrow(SvSigs), byrow = F)
  sid = toSid(x, s)
  if (length(sid) == 1 && (sid<0 || sid>x$maxSid)) { #nrow(x$S)))
    print(paste0('getState: Invalid sid=',sid,', SValues=',toString(SValues)))
    return(NA)
  }
  S <- tibble()
  for (i in 1:length(sid)) S <- bind_rows(S, x$S)
  S$sid <- sid
  S[,x$xSigCols] <- s
  S[,x$vSigCols] <- s*x$baseX/x$sSteps
  S
}

mapState.finDemo <- function(x, X = NULL, ...) {
  if (!is.null(X) && is.numeric(X)) {
    si <- x$XRev[findInterval(X,pull(x$XRev,Si)),x$vCols]
    names(si)<-x$xColNames
    si
  } else {
    b <- x$B1+x$B2
    for (i in 3:8) {
      if (has_name(x, paste0('B',i))) b <- b + x[[paste0('B',i)]]
    }
    b <- scale(b)
    b <- b - min(b) # offset to >=0
    # # negative States
    tibble(Si = as.numeric(b/max(abs(c(max(b),min(b)))) * .95*x$baseX)) %>%
      mutate(S = as.numeric(findInterval(Si, x$valSteps)))
  }
}

toQid <- function(sid, aid, aRows) {
  aRows*(sid-1)+aid
}

#' calculate S-index by Si and sig.p
#' @example sum(x$S$sid-toSid(x, x$S[,x$xCols])) == 0
toSid <- function(x, s, itv = x$itv) {
  m <- as.matrix(s)-1
  m <- m[,1:x$xSigColNum] # remove S-Constant column
  as.numeric(m %*% itv^((x$xSigColNum-1):0) + 1)
}

move <- function(x, SValues, A, t = x$t, k = 1) {
  if (is.numeric(A)) A <- x$A[A,]
  if (is.vector(SValues)) SValues <- matrix(SValues, nrow=1)

  SvSigs <- SValues[,1:x$xSigColNum]
  m.v <- A[,x$xSigCols] * k * x$baseX/x$sSteps
  svals <- t(apply(m.v, 1, `+`, as.numeric(SvSigs)))

  if (deblog) { print('move: m.v:'); print(m.v) }
  if (deblog) print(paste0('move: SValues=',toString(SvSigs)))
  if (deblog) print(svals)
  if (deblog) { print('move: x[t:t+1]:'); print(x$X[t:(t+1),]) }
  if (deblog) { print('move: diff:'); print(x$X[t:(t+1),x$xSigCols]-as.numeric(SvSigs)) }

  if (any(svals<0 | svals>x$baseX)) {
    m.v <- A[,x$xSigCols] * 1 * x$baseX/x$sSteps
    svals <- t(apply(m.v, 1, `+`, as.numeric(SvSigs)))
  }

  svals <- replace(svals, svals<0 | svals>x$baseX, NA)
  svals <- as.data.frame(svals)
  svals$SC <- SValues$SC

  S.moves <- getState(x, svals)
  if (deblog) print(paste0('move[',x$run,'/',x$i,']: A=', toString(A)))
  if (deblog) { print(paste0('move[',x$run,'/',x$i,']: SValues:')); print(svals) }
  if (deblog) print(paste0('move[',x$run,'/',x$i,']: next=',toString(S.moves$sid)))
  list(S.moves = S.moves, SValues = svals)
}

#'Error function providing linear positive gradient and quadratic center
#'@usage
#'fe <- tibble(e=seq(-40,40,.01),y=fErr(x, seq(-40,40,.01))); ggplot(fe)+geom_line(aes(e,y))+
#'   geom_hline(aes(yintercept=.5), linetype='longdash')+
#'   geom_vline(aes(xintercept=sqrt(x$fErr.a*x$fErr.base/2)), linetype='longdash')+
#'   geom_vline(aes(xintercept=-sqrt(x$fErr.a*x$fErr.base/2)),linetype='longdash')
#'   geom_hline(yintercept = x$fErr.v0/2/x$fErr.base, linetype='dashed')
fErr <- function(x, err, v0 = x$fErr.v0, base = x$fErr.base, a = x$fErr.a) {
   ifelse(abs(err) < sqrt((a*base)/2),
          v0 * (1 - err^2/(a*base)),
          v0*a*base/4 / err^2
   ) / base
}

estimateR.all <- function (x, SValues, t = x$t+1) {
  apply(SValues, 1, function(sv) { estimateR(x, sv, t) })
}

estimateR <- function (x, SValues, t = x$t+1) {
  W <- x$W
  xn <- x$vCols

  if (!is.na(x$fs.pred) && 'S'==x$fs.pred) {
    FS.S <- x$FS$Si[max(1,t-1)] + sum((x$W-x$W.prev)*SValues + (x$SValues-x$SValues.prev)*x$W)
    FS.t <- ifelse(t>2, 2*x$FS$Si[t-1] - x$FS$Si[t-2], x$FS$Si[max(1,t-1)])
    if (x$i > 50 || is.na(FS.S) || FS.S < 0 || FS.S > x$baseX) {
      FS <- FS.S
    } else {
      FS <- (FS.t*(50-x$i) + FS.S*x$i)/50
      if (FS < 0 || FS > x$baseX) FS <- FS.t
    }
    if (FS < 0) FS <- 0
    if (FS > x$baseX) FS <- x$baseX
    print(paste0('estimateR: FS=',toString(FS),', FS.t=',x$FS$Si[t],', dW=',toString(x$W-x$W.prev),', dSV=',toString(x$SValues-x$SValues.prev)))
  } else if (T) { # predict
    FS <- ifelse(t>2, 2*x$FS$Si[t-1] - x$FS$Si[t-2], x$FS$Si[max(1,t-1)])
  } else {
    FS <- x$FS$Si[t]
  }

  ws <- SValues * W
  norm.ws <- sqrt(sum(ws^2))
  if (is.na(norm.ws) || norm.ws==0) print(paste0('WARN: W=',toString(W),', SValues=',toString(SValues)))
  if (deblog) print(paste('SValues=',toString(SValues),', FS=',toString(FS),', W=',toString(W),', ws=',toString(ws)))

  ifelse(!is.na(norm.ws) && norm.ws > 0, fErr(x, (norm.ws - FS)), -x$baseX*x$numSigs)
}

#' R: reward of the State/Action pair
#' t: time t in State.next (after action A)
calculateR <- function(x, SValues, A, SValues.prev = x$SValues.prev, t = x$t) {
  W <- x$W
  FS <- x$FS$Si[t]
  ws <- SValues * W
  norm.ws <- sqrt(sum(ws^2))
  norm.ws - FS
}

#' differenced Value-Function
#' dfs <- seq(-50,50,.1);
#' ggplot()+geom_line(aes(dfs,dFerr(dfs,x$fErr.a,x$fErr.base)))+
#'   geom_vline(xintercept = -sqrt(x$fErr.a*x$fErr.base/2), linetype='dashed')+
#'   geom_vline(xintercept = sqrt(x$fErr.a*x$fErr.base/2), linetype='dashed')+ylab("f'(err)")+xlab('err')
dFerr <- function(dFs, a, base) {
   ifelse(abs(dFs) < sqrt(a*base/2),
          -2*dFs/(a*base),
          -(a*base)/(2*dFs^3))
}

calculateDx <- function(x, SValues, A, SValues.prev = x$SValues.prev, t = x$t, v0 = x$fErr.v0, base = x$fErr.base, a = x$fErr.a, alpha = x$alpha) {
  norm.ws <- sqrt(sum((x$W*SValues)^2))
  Fs <- x$FS$Si[t]
  deltaFs <- norm.ws-Fs

  dx <- ifelse(norm.ws == 0, 1, dFerr(deltaFs, a, base)) * v0/a/base * ((x$W*SValues^2) / norm.ws)

  print(paste0('calculateDx[',x$run,'/',x$i,']: t=',x$t,' dx=',toString(dx),', SValues=',toString(SValues),', W=',toString(x$W)))

  dx
}

#' get the real State/Action for the possibly updated S and A
#' @param state S as tibble/list
#' @param action A as tibble/list or factor
#' @return a merged object of S and A
#' @seealso getState
getStateAction.finDemo <- function(x, S, A, ...) {
  if (is.null(x) || is.null(S) || is.null(A)) stop('getStateAction: missing input')
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
  if (deblog) print(paste0('possible[',x$run,'/',x$i,']: t=',x$t,', SValues=', toString(SValues)))
  if (is.null(x) || is.null(S)) stop('possibleActions: missing input')
  if (is.numeric(S)) {
    S <- x$S[S,]
  }

  SValues.in <- SValues

  xA <- x$A

  m <- move(x, SValues, xA)
  good <- !is.na(m$S.moves$sid)
  if (sum(good) == 0) {
    print(SValues.in)
    print(xA)
    print(S.moves)
    print(SValues)
    print(m$S.moves)
    print(m$SValues)
    print(paste0('no good Action in State: ',toString(S)))
    return(NULL)
  }
  S.moves <- m$S.moves
  SValues = m$SValues

  xA <- xA[good,]
  S.moves <- S.moves[good,]
  SValues <- SValues[good,]
  if (deblog) { print(paste0('possible[',x$run,'/',x$i,']: t=',x$t,', SValues.good:')); print(SValues) }

  R <- estimateR.all(x, SValues)
  xA <- bind_cols(xA, R=R)

  xA$Q <- R

  if (nlog) {
    print(paste0('possible[',x$run,'/',x$i,']: t=',x$t,', S=',toString(list(S))))
  }

  # by: Value increase by Delta-X, sell: Value decreased by Delta-X:
  po <- xA %>% arrange(desc(Q))
  print(paste0('possible[',x$run,'/',x$i,']: t=',x$t,', po=',toString(head(po,1))))

  if (deblog) print(tibble(m=paste0('possible[',x$run,'/',x$i,']:'), po), width=Inf)
  if (deblog) x$n.t <- nanoTime(paste('possible',x$numSigs), n.t)

  po
}

#' plot(2*x$epsilon*(1-e1071::sigmoid(((1:5000)-5000)/1000)))
#' t<-1:7000; d<-tibble(t,y=getEpsilon(x,t)); ggplot(d,aes(t,y))+geom_line()+ylab('Epsilon')
getEpsilon.finDemo <- function(x, i = x$i) {
  2*x$epsilon*(1-e1071::sigmoid((i-5000)/500))
}

#' Updates state S according action A (as tibble/list or as factor/character)
#' @param state S as tibble/list
#' @param action A as tibble/list
#' @return updated real (S,A,SValues) object
evalAction.finDemo <- function(x, S, A, SValues = x$SValues, ...) {
  if (deblog) n.t <- nanoTime(paste('evalAction',x$numSigs), x$n.t)
  if (is.null(x) || is.null(S) || is.null(A)) stop('evalAction: missing input')
  if (is.factor(A)||is.character(A)) { print(A); stop('evalAction: A cannot be primitive') }
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$A[A,]

  AIn <- A
  dMFee <- 0
  SValues.in <- SValues

  R <- rep(0,x$xColNum)

  if (all(A[x$xCols] == 0)) { # no action
    sv <- SValues + x$alpha * (SValues-x$SValues.prev)
    if (all(sv > 0 & sv < x$baseX)) {
      SValues <- sv
      print(paste0('ACTION[',x$run,'/',x$i,']: t=',x$t,', ZeroAction: SValues=',toString(SValues)))
    }
    R <- calculateR(x, SValues, A, SValues.in, t=x$t+1)

  } else { # try to trade
    m <- move(x, SValues, A)
    S.move <- m$S.move
    SValues <- m$SValues

    if (is.null(S.move) || nrow(S.move) == 0 || is.na(S.move$sid)) {
      R <- calculateR(x, SValues, A, SValues.in, t=x$t+1)
      A <- NULL
    } else {
      S.next <- S.move
      R <- calculateR(x, SValues, A, SValues.in, t=x$t+1)
      S <- S.next
    }
  }

  if (!is.null(A)) {
    print(paste0('ACTION[',x$run,'/',x$i,']: t=',x$t,', A=',toString(list(A))))
    print(paste0('ACTION[',x$run,'/',x$i,']: t=',x$t,', S=',toString(list(S)),', SValues=',toString(SValues),'/',toString(x$SValues)))
  } else {
    if (nlog) print(paste0('DENIED[',x$run,'/',x$i,']: t=',x$t,', S=',toString(list(S))))
    if (nlog) print(paste0('DENIED[',x$run,'/',x$i,']: t=',x$t,', A=',toString(list(AIn)),', SValues=',toString(SValues),'/',toString(x$SValues)))
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
  if (is.null(x) || is.null(S) || is.null(A)) stop('act: missing input')
  if (is.numeric(S)) S <- x$S[S,]
  if (is.numeric(A)) A <- x$A[A,]

  terminal <- F
  x$resA <- NULL
  SIn <- S
  x$trace <- bind_rows(x$trace, c(t = x$t, A = A))

  res <- evalAction(x, S, A, SValues)
  S <- res$S
  SValues <- res$SValues

  if (F&&sum(SValues) <= 1) { # BROKE!! end of episode
    print(paste0('Broke[',x$run,'/',x$i,']! t=',x$t,', S=',toString(list(S)),', SValues=',toString(SValues)));
    R <- rep(-1000, x$xColNum)
    terminal = T
  } else {
    R <- res$R
    if (!is.null(res$A)) {
      A <- res$A
      A$sid <- res$S$sid
    }
  }

  R.r <- R
  x$t <- x$t + 1
  if (x$t >= x$tSteps) { # end of episode!
    terminal = T
    R <- x$SumR # 0 # True Online SARA # x$SumR
  } else {
    x$SumR <- R + x$SumR
  }

  if (terminal) {
    print(paste0('END[',x$run,'/',x$i,']: t=',x$t,', R=',toString(R),', SValues=',toString(x$SValues),', S=',toString(list(S))))
  }

  A$in.sid <- SIn$sid
  A$sid <- S$sid
  A$terminal <- terminal
  A$t <- x$t

  x$A$t[A$aid] <- x$t
  x$SValues.prev <- x$SValues
  x$SValues <- SValues
  x$resA <- list(A=A, R=R, S = S, SValues = SValues, R.r = R.r)
  if (!terminal & nlog) {
    print(paste0('act[',x$run,'/',x$i,']: t=',x$t,', resA=',toString(list(x$resA)),', SIn=',toString(list(SIn))));
  }
  print(paste0('act[',x$run,'/',x$i,']: t=',x$t,', R=',toString(R.r),', SumR=',toString(x$SumR),', SValues=',toString(x$resA$SValues),', S=',toString(list(S))))
  if (deblog) x$n.t <- nanoTime(paste('act',x$numSigs), n.t)
  return (x)
}

#' check if acting out action A in state S will leave the world's context
#' @param state S as tibble/list
#' @param action A as tibble/list or factor/character
#' @return true if next state would be outside world's context
#' @seealso evalAction
offLimits.finDemo <- function(x, S, A, ...) {
  if (is.null(x) || is.null(S) || is.null(A)) stop('offLimits: missing input')
  S <- evalAction(x, S, A)
  is.null(S) || nrow(S) == 0
}

#' penalize koefficients if exceeding interval {-1,1}
#' @see w <- seq(-2,5,.1); ggplot(tibble(w,y=penalize(w)),aes(w,y))+geom_line()
#' @see s <- seq(-4,4,.1);ggplot()+geom_line(aes(s,5*(-.5+e1071::sigmoid(s))))+geom_line(aes(s,penalize(s,a=1,d=5)))
penalize <- function(w, k=4) {
  k*(-.5+e1071::sigmoid(w))
}

#' prohibit parameter decay
#' @see s<-seq(-.2,.5,.01);m<-matrix(c(rep(1,length(s)),s),ncol=2);
#' md<-t(apply(m,1,deDegenerate)); ggplot()+geom_line(aes(s,md[,1]))+geom_line(aes(s,md[,2]))
deDegenerate <- function(w, k=.2, w0=.02) {
  sw <- w / max(abs(w))
  w + k * (-.5+e1071::sigmoid(sw/w0))
}

#' update W for this action A
#' @param action A as tibble/list (match is done by A$aid) or as aid
#' @param return R from last action
#' @return the environment
updateQ.finDemo <- function(x, A, R, Anext = NULL, alpha = x$alpha, gamma = x$gamma, lambda = .2, n = 8, ...) {
  if (deblog) n.t <- nanoTime(paste('updateQ',x$numSigs), x$n.t)
  if (deblog) {
    print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', A=',toString(list(A)),', R=',toString(R)));
    if(!is.null(Anext)) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', Anext=',toString(Anext)));
  }
  if (is.null(x) || is.null(A) || is.null(R)) stop('updateQ: missing input')
  if (is.numeric(A)) A <- x$A[A,]
  if (is.numeric(Anext)) Anext <- x$A[Anext,]

  aid <- A$aid
  R.orig <- R
  R <- fErr(x, x$SumR)
  if (deblog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', R.orig=',toString(R.orig),', R=',toString(R)))

  Rc <- R
  failed <- !is.null(x$resA) && has_name(x$resA$A,'failed') && x$resA$A$failed
  W <- x$W
  dQ <- NA
  Qa <- NA
  Qa.next <- NA
  dQ.next <- NA
  delta <- NA

  # Update W
  dQ <- as.numeric(calculateDx(x, x$SValues, A, x$SValues.prev, t = x$t))
  Qa <- estimateR(x, x$SValues, t = x$t)
  Qa.next <- NA
  dQ.next <- NA
  delta <- NA

  if (!failed && !A$terminal) {
    if (is.null(Anext)) {
      if (F) { # True Online SARSA: no update in S==terminal
        dW <- (sum(Rc - Qa)) * dQ #/ x$baseX^2
        W <- W + alpha * dW # / max(1,abs(as.numeric(dW)))
        if (nlog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', dQ=',toString(dQ),', Qa=',toString(Qa),', dW=',toString(dW)))
      } else {
        if (nlog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', Anext=',toString(Anext)))
      }
    } else {

      m <- move(x, x$SValues, Anext)
      if (!any(is.na(m$SValues))) {

        dQ.next <- as.numeric(calculateDx(x, m$SValues, Anext, x$SValues, t = x$t+1))
        Qa.next <- estimateR(x, x$SValues)

        if (T|deblog) print(paste0('     dQ=',toString(dQ),', Qa=',toString(Qa)))
        if (T|deblog) print(paste0('dQ.next=',toString(dQ.next),', Qa.next=',toString(Qa.next)))

        delta <- Rc + gamma * Qa.next - Qa
        if (deblog) print(paste0('  delta=',toString(delta),' Z=',toString(x$Z)))

        Ze <- x$Z
        x$Z <- gamma*lambda * x$Z + penalize((1 - alpha*gamma*lambda*sum(x$Z*dQ)) * dQ *.5, k=2)
        if (max(abs(x$Z)) > 2) {
          print(paste0('recalibrated: Z=',toString(x$Z),', Ze=',toString(Ze)))
        }

        dw.dqz <- as.numeric(sum(delta + Qa - x$Qold) * x$Z)
        dw.qz <- as.numeric(sum(Qa - x$Qold) * dQ) # / x$baseX

        We <- W
        Wa <- (alpha * (dw.dqz - dw.qz)) #/ 10 # not too much volatility
        W <- (W + Wa)
        if (max(abs(W)) > .8) {
          Wp <- W
          W <- penalize(W/1.5, k=6)
          print(paste0('penalized: W=',toString(W),', Wp=',toString(Wp),', Wa=',toString(Wa)))
        }
        if (min(abs(W)) < .01) {
          Wd <- W
          W <- deDegenerate(W)
          print(paste0('deDegenerated: W=',toString(W),', Wd=',toString(Wd),', Wa=',toString(Wa)))
        }

        if (T|deblog) print(paste0('   Qold=',toString(x$Qold)))
        if (T|deblog) print(paste0('      Z=',toString(x$Z),', Ze=',toString(Ze)))
        if (T|deblog) print(paste0(' dw.dqz=',toString(dw.dqz)))
        if (T|deblog) print(paste0('  dw.qz=',toString(dw.qz)))
        if (T|deblog) print(paste0('      W=',toString(W),', Wa=',toString(Wa),', We=',toString(We)))

        if (nlog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', dQ=',toString(dQ),', Qa=',toString(Qa),', delta=',toString(delta),', dQ=',toString(dQ)))
        if (nlog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', dQ.n=',toString(dQ.next),', Qa.n=',toString(Qa.next),', Z=',toString(x$Z),', Qold=',toString(x$Qold)))
        if (nlog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', dw.dqz=',toString(dw.dqz),'/',toString(dw.dqz/max(1,abs(dw.dqz))),', dw.qz=',toString(dw.qz),'/',toString(dw.qz/max(1,abs(dw.qz)))))
        x$Qold <- Qa.next

      } else {
        if (T|nlog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', invalid move: SValues=',toString(x$SValues),', Anext=',toString(Anext)))
      }
    }

    if (nlog) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', x$W=',toString(x$W),', W=',toString(W)))
  } else {
    if (nlog) if (failed) print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,' failed: x$W=',toString(x$W),', W=',toString(W)))
  }

  po <- x$chosenAction
  lid <- nrow(x$log)+1
  x$log <- add_row(x$log, run = x$run, lid = lid, time = microbenchmark::get_nanotime())
  x$log$i[lid] <- x$i
  x$log$t[lid] <- x$t
  x$log$sid[lid] <- A$sid
  x$log$aid[lid] <- A$aid
  x$log[lid,x$logSVals.cols] <- x$SValues
  x$log[lid,x$logR.orig.cols] <- x$resA$R.r
  x$log[lid,x$logR.cols] <- R
  x$log[lid,x$logSumR.cols] <- x$SumR
  x$log[lid,x$logdQ.cols] <- dQ
  x$log[lid,x$logdQn.cols] <- dQ.next
  x$log[lid,x$logQa.cols] <- Qa
  x$log[lid,x$logQan.cols] <- Qa.next
  x$log[lid,x$logDelta.cols] <- delta
  x$log[lid,x$logW.cols] <- W
  x$log[lid,x$logZ.cols] <- x$Z
  x$log$Qpo[lid] <- po$Q
  x$log$epsilon[lid] <- x$currEpsilon

  if (nlog) {
    print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', A=',toString(list(A))))
    print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', Rc=',toString(Rc),', Qa=',toString(Qa),', Qa.next=',toString(Qa.next)))
  }

  wd <- as.numeric(W)
  if (any(is.na(wd) | is.nan(wd) | wd == Inf | wd == -Inf)) {
    print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', A=',toString(A),', x$W=',toString(x$W)));
    print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', R=',R,', Rc=',toString(Rc)));
    print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', Qa=',toString(Qa),', W=',toString(W)));
    stop(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,' invalid W: ', toString(W)))
  }
  if (any(x$SValues < 0)) {
    print(paste0('updateQ[',x$run,'/',x$i,']: t=',x$t,', aid=',aid,', SValues:',toString(x$SValues)));
    print('Invalid SValues: going SHORT!')
  }

  x$W.prev <- x$W
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


######################### The FinDemo ##########################
##############################################################
