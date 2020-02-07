##############################################################
source('worlds.R')
source('utils.R')
######################### Finance Demo #######################
sSteps = 20; tSteps = 20; baseX = 100; X = NULL; deltaX = .1;
ica=T; numSigs=2; numBase = 7; numOrig = 4
fErr.v0 = baseX; fErr.base = baseX; fErr.a = 1; fs.pred='S'

finDemo <- function(sSteps = 10, tSteps = 20, baseX = 100,
                    epsilon = .1, alpha = .1, gamma = .9, X = NULL, deltaX = .1,
                    ica = T, numSigs = 2, numBase = 7, numOrig = 3, ica.ratio = 10, ica.maxRetries = 500, ica.high = 20,
                    fErr.v0 = baseX, fErr.base = baseX, fErr.a = 1, fs.pred = 'S', ...) {

  x <- new_finDemo(sSteps = sSteps, tSteps = tSteps, baseX = baseX,
                   epsilon = epsilon, alpha = alpha, gamma = gamma,
                   X = X, deltaX = deltaX, ...)
  x$numSigs <- numSigs
  x$ica <- ica
  x$numBase <- numBase
  x$numOrig <- numOrig
  x$fErr.v0 <- fErr.v0
  x$fErr.base <- fErr.base
  x$fErr.a <- fErr.a
  x$fs.pred <- fs.pred
  eps.t <- .005
  eps.x <- .01

  x$B1 <- .5*(1+sin(2*pi*1.5*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
  x$B2 <- .5*(1+cos(2*pi*2.5*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
  if (x$numOrig > 2) {
    x$B3 <- .5*(1+cos(pi/5+2*pi*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
  }
  if (x$numOrig > 3) {
    x$B4 <- .5*(1+sin(pi/4+2*pi*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
  }
  if (x$numOrig > 4) {
    x$B5 <- .5*(1+cos(pi/5+3*pi*(1:x$tSteps)/x$tSteps + runif(x$tSteps,-eps.t,eps.t))+runif(x$tSteps,-eps.x,eps.x))
  }
  if (is.null(x$X)) {
    if (numOrig == 5) {
      x$X <- tibble(S1 = .5*(x$B1+x$B2), S2 = 1+.5*(x$B2-x$B3), S3 = .7*(x$B3-x$B4), S4 = .5*(x$B4-x$B5), S5 = .5*(x$B5-x$B1))
    } else if (numOrig == 4) {
      x$X <- tibble(S1 = .5*(x$B1+x$B2), S2 = 1+.5*(x$B2-x$B3), S3 = .7*(x$B3-x$B4), S4 = .5*(x$B4-x$B1))
    } else if (numOrig == 3) {
      x$X <- tibble(S1 = .5*(x$B1+x$B2), S2 = 1+.5*(x$B2-x$B3), S3 = .7*(x$B3-x$B1))
    } else {
      x$X <- tibble(S1 = x$B1, S2 = x$B2)
    }
    for (i in (x$numOrig+1):x$numBase) x$X <-
        bind_cols(x$X, !!paste0('S',i) :=
                    i/10 * x$B1 * (1+runif(x$tSteps,-eps.t,eps.t)) +
                    (6-i)/10 * x$B2 * (1+runif(x$tSteps,-eps.x,eps.x))
                  )
    x$X <- .95*x$baseX * scale(x$X) # v5:/ ncol(x$X) # X just a scaling factor?
    x$X.orig <- as_tibble(x$X)
  }
  print(paste('Signals generated:',object.size(x)))

  #ggplot(as_tibble(x$X) %>% rowid_to_column('t') %>%  gather(key,value,-t))+geom_line(aes(t,value, color=key))

  if (x$ica) {

    ica.retries = 0
    ica.bestX <- NULL
    ica.bestRatio <- 10000

    repeat {

      x$xpca <- prcomp(x$X.orig)
      x$X.pca <- as.matrix(x$X.orig) %*% x$xpca$rotation[,1:min(x$numSigs,dim(x$xpca$rotation)[2])]
      x$xica <- fastICA::fastICA(x$X.orig, n.comp = x$numSigs, maxit=50, tol=.01, fun = 'exp',
                                 row.norm = T,
                                 method = 'C', alg.typ = 'parallel', verbose = T)
      print('ICA calculated..')
      #x$X <- sweep(x$xica$S,2, apply(x$xica$S,2,min))
      #ica.s <- apply(x$xica$S, 2, function(s) { s <- s-min(s); s/max(s) * x$baseX*.95 })
      #ica.s <- apply(x$xica$S, 2, function(s) { s/(max(s)-min(s)) * 2*x$baseX*.95 })
      # Rescaling to positive values:
      #   - move restricting to positive gradient -> would move to negative values
      x$X <- as_tibble(x$xica$S, .name_repair = ~paste0('S',1:ncol(x$xica$S))) %>%
        #transmute_all(list(~(.-min(.))/max(abs(c(max(.),min(.))))/2 * x$baseX*.95))
        transmute_all(list(~./max(abs(c(max(.),min(.)))) * x$baseX*.95)) # negative States
        #names(x$X) <- paste0('S',1:ncol(x$X))
      #x$X <- as_tibble(x$X / (max(x$X)-min(x$X)) * .98*x$baseX) # v5: / ncol(x$X))

      freqRes <- eval_freq2(x, high = ica.high)
      if (freqRes$ratio <= ica.ratio) {
        break
      }
      if (freqRes$ratio < ica.bestRatio) {
        ica.bestRatio <- freqRes$ratio
        ica.bestX <- x$X
      }
      if (ica.retries > ica.maxRetries) {
        x$X <- ica.bestX
        print(paste('ICA BestRatio: ', ica.bestRatio))
        break
      } else {
        print(paste('retrying ICA: ratio',freqRes$ratio,'best',ica.bestRatio,'retry',ica.retries))
      }
      ica.retries <- ica.retries+1
    }
  }
  #ggplot(as_tibble(x$X) %>% rowid_to_column('t') %>%  gather(key,value,-t))+geom_line(aes(t,value, color=key))

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
  #x$valSteps <- x$baseX * -(x$sSteps):(x$sSteps)/x$sSteps # negative States
  x$SIntv <- tibble(sid=1:x$itv)
  #x$SIntv <- tibble(sid=1:(2*x$itv-1)) # negative States
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
  #for (n in x$xSigColNames) x$Si <- crossing(x$Si, !!n := x$SIntv[n])
  for (n in x$xSigColNames) x$Si <- mutate(x$Si, !!n := 0)
  x$Si <- x$Si %>% select(-sid)
  names(x$Si)<-paste0(names(x$Si),'i')
  x$Si <- x$Si %>% mutate(SCi = 0) %>% rowid_to_column('sid')

  x$S <- tibble(sid = 0)
  #for (n in x$xColNames) x$S <- crossing(x$S, !!n := 1:(2*x$itv+1))
  #for (n in x$xSigColNames) x$S <- crossing(x$S, !!n := 1:x$itv)
  for (n in x$xSigColNames) x$S <- mutate(x$S, !!n := 1)
  x$S <-  x$S %>% mutate(SC = 0) #%>% rowid_to_column('sid')
  #x$S <- crossing(x$S, sig.p = c(-1,1))
  x$S <- bind_cols(x$S, x$Si %>% select(-sid)) %>% mutate(n=0)
  #x$S['spid'] <- 1:nrow(x$S)
  x$S <- x$S %>% select(-sid) %>% rowid_to_column('sid')
  x$sCols <- which(names(x$S) %in% x$xColNames)
  x$vCols <- which(names(x$S) %in% paste0(x$xColNames,'i'))
  x$vSigCols <- which(names(x$S) %in% paste0(x$xSigColNames,'i'))
  x$maxSid <- x$itv^x$xSigColNum
  #x$maxSid <- (2*x$itv-1)^x$xSigColNum # negative States
  print(paste('States generated:',object.size(x)))
  print(x$S)

  x$FS <- mapState(x) %>% mutate(p = Si/x$baseX, dp = c(0,diff(Si))) %>% mutate(sig.p = sign(dp))

  x$XRev <- x$X %>% arrange(Si)
  #for (n in xColNames) { XRev[paste0(n,'.i')] <- findInterval(pull(XRev,n),pull(SIntv,1)) }

  sv <- x$X[1,x$xCols] # /sum(abs(x$X[1,x$xCols]))*x$FS$Si[1] # v3: x$X$S[1] + x$SIntv[2,]
  sv$SC <- 0
  sv$SC <- x$baseX/x$numOrig
  #x$W.init <- rep(1/x$xColNum, x$xColNum) * replace(sign(sv),sv==0,1)
  # correlation of signals:
  xb <- bind_cols(x$X %>% select(matches('^S\\d$')), FS=x$FS$Si)
  x$W.init <- c((corrr::correlate(xb))$FS[1:x$xSigColNum],1/x$numSigs)
  x$W.init <- x$W.init/max(abs(x$W.init))/x$numSigs
  x$W <- x$W.init
  x$W.prev <- x$W
  x$SValues.init <- abs(sv)
  #x$SValues.init <- sv # negative States
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
  #q <- crossing(sid=x$S$sid,aid=x$A$aid) %>% mutate(qid = toQid(sid,aid,x$aRows))
  #x$Q <- tibble(aid = x$A$aid)
  #for (n in x$xColNames) x$Q[n] <- 0
  #x$Q <- bind_cols(x$Q, sid=q$sid,aid=q$aid)

  x$Z <- rep(0, x$xColNum)
  x$Qold <- 0

  x$brokeR <- -10000
  x$offLimR <- -1
  x$A$dV <- 0
  x$SumR <- .0
  x$numSigs <- numSigs
  x$R <- 0
  x$t <- 1
  x$i <- 0
  x$currEpsilon <- getEpsilon(x)
  x$trace <- c()
  x$n.t <- microbenchmark::get_nanotime()
  x$logSVals.cols <- paste0('SVal.',x$xColNames)
  x$logR.cols <- paste0('R.',x$xColNames)
  x$logQpo.cols <- paste0('Qpo.',x$xColNames)
  x$logdQ.cols <- paste0('dQ.',x$xColNames)
  x$logdQn.cols <- paste0('dQn.',x$xColNames)
  x$logQa.cols <- paste0('Qa.',x$xColNames)
  x$logQan.cols <- paste0('Qan.',x$xColNames)
  x$logDelta.cols <- paste0('delta.',x$xColNames)
  x$logW.cols <- paste0('W.',x$xColNames)
  x$logZ.cols <- paste0('Z.',x$xColNames)
  x$log <- tibble(lid=integer(), time=numeric(), i=integer(), t=integer(), sid=integer(), aid=integer(), Qpo=numeric(), epsilon=numeric())
  for (n in c(x$logSVals.cols,x$logR.cols,x$logQpo.cols,x$logdQ.cols,x$logdQn.cols,
              x$logQa.cols,x$logQan.cols,x$logDelta.cols,x$logW.cols,x$logZ.cols)) x$log[n] <- numeric()
  print(paste('Environment generated:',object.size(x)))
  x
}

new_finDemo <- function(..., class = character()) {
  new_world(..., class = c(class, 'finDemo')
  )
}

#' initialize environment (e.g. reset Q to zero)
init.finDemo <- function(x, ...) {
  x$W <- x$W.init
  x$W.prev <- x$W
  x <- getStart(x)
  x
}

#' get the starting state of this environment
getStart.finDemo <- function(x, ...) {
  x$t = 1
  x$currEpsilon <- getEpsilon(x)
  #x$SValues.init <- updateInitialState(x)
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
  if (nrow(x$log) > 10) {
    S <- getState(x, x$SValues.init)
    if (nrow(S) == 0) {
      stop(paste0('updateInitial: bad InitialState: init=',toString(x$SValues.init)))
    }
    po <- possibleActions(x, S, x$SValues.init)
    if (is.null(po)) {
      print(paste0('updateInitial: no Action: init=',toString(x$SValues.init)))
      return(x$SValues.init)
    }

    if (deblog) print(S)
    if (deblog) print(po)
    topA <- po %>% arrange(desc(Q)) %>% head(1)
    if (sum(topA[x$xCols])==0) {
      #dX <- calculateDx(x, x$SValues.init, topA, x$SValues.init)/x$baseX
      #svinit <- x$SValues.init + x$alpha * dX
      #print(paste0('updateInitial: SValues=',toString(svinit),', init=',toString(x$SValues.init),', dX=',toString(dX)))
      if (deblog) print(paste0('updateInitial: init=',toString(x$SValues.init),', topA=',toString(topA)))
      return(x$SValues.init)
    }
    m <- move(x, x$SValues.init, topA, k = 1)
    svinit <- x$SValues.init + x$alpha * (m$SValues-x$SValues.init)
    S <- getState(x, svinit)
    if (nrow(S)>0 && !any(is.na(S$sid))) {
      print(paste0('updateInitial: SValues=',toString(svinit),', init=',toString(x$SValues.init),', m=',toString(m$SValues),', topA=',toString(topA)))
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
  #SValues <- replace(SValues, SValues<0, NA)
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
  #S[,x$vSigCols] <- (s-x$itv)*x$baseX/x$sSteps # negative States
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
  #as.numeric(m %*% (2*itv-1)^((x$xSigColNum-1):0) + 1) # negative States
}

move <- function(x, SValues, A, t = x$t, k = 2) {
  if (is.numeric(A)) A <- x$A[A,]
  if (is.vector(SValues)) SValues <- matrix(SValues, nrow=1)

  SvSigs <- SValues[,1:x$xSigColNum]
  m.v <- A[,x$xSigCols] * k * x$baseX/x$sSteps
  svals <- t(apply(m.v, 1, `+`, as.numeric(SvSigs)))

  if (deblog) print(m.v)
  if (deblog) print(paste0('move: SValues=',toString(SvSigs)))
  if (deblog) print(svals)
  if (deblog) print(x$X[t:(t+1),])
  if (deblog) print(x$X[t:(t+1),x$xSigCols]-as.numeric(SvSigs))

  #if (any(svals < -x$baseX | svals > x$baseX)) { # negative States
  if (any(svals<0 | svals>x$baseX)) {
    m.v <- A[,x$xSigCols] * 1 * x$baseX/x$sSteps
    svals <- t(apply(m.v, 1, `+`, as.numeric(SvSigs)))
  }

  # Action in direction of gradient of X
  #dS.sv <- t(apply(A[,x$xCols], 1, `*`, as.numeric(x$X[t+1,x$xCols]-x$X[t,x$xCols])*x$baseX/x$sSteps))

  # Action in direction of SValues torwards X[t+1]
  if (F && nrow(A) > 1) {
    dS.X <- as.numeric(x$X[t+1,x$xSigCols] * x$W[1:x$xSigColNum] - SvSigs)
    dS.sv <- t(apply(A[,x$xSigCols], 1, `*`, dS.X * x$baseX/x$sSteps))
    poss <- rowSums(dS.sv) >= 0
    svals[!poss,] <- NA
  }

  svals <- replace(svals, svals<0 | svals>x$baseX, NA)
  #svals <- replace(svals, svals < -x$baseX | svals > x$baseX, NA) # negative States
  svals <- as.data.frame(svals)
  svals$SC <- SValues$SC
  if (deblog) print(svals)

  S.moves <- getState(x, svals)
  if (deblog) { print(paste0('move[',x$i,']: A')); print(A) }
  if (deblog) { print(paste0('move[',x$i,']: SValues')); print(SValues); } # print(svals) }
  if (deblog) print(paste0('move[',x$i,']: next=',toString(S.moves$sid)))
  list(S.moves = S.moves, SValues = svals)
}

#'Error function providing linear positive gradient and quadratic center
#'@usage
#' fe <- tibble(e=seq(-3,3,.01),y=fErr(x, seq(-3,3,.01))); ggplot(fe)+geom_line(aes(e,y))
#'   +geom_hline(aes(yintercept=(1-1/(2*x$numSigs*x$baseX))), linetype='longdash')+
#'    geom_vline(aes(xintercept=1), linetype='longdash')+
#'    geom_vline(aes(xintercept=-1),linetype='longdash')
#' fe <- tibble(e=seq(-20,20,.01),y=fErr(x, seq(-20,20,.01))); ggplot(fe)+geom_line(aes(e,y))
#'   +geom_vline(xintercept = x$fErr.a*x$fErr.base/2, linetype='dashed')+
#'    geom_vline(xintercept = -x$fErr.a*x$fErr.base/2, linetype='dashed')+
#'    geom_hline(yintercept = x$fErr.v0/2/x$fErr.base, linetype='dashed')
#fErr <- function(x, err, v0 = x$fErr.v0, base = x$fErr.base, a = x$fErr.a) {
#  ifelse(abs(err) > 1,
#         ifelse(abs(err) < (a*base)/2,
#                v0 * (1 - err/(a*base)),
#                a*base*v0/4 / err
#                ),
#         v0 * (1 - (1/(2*a*base) * (1 + err^2)))
#  ) / base
#}

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
  #p <- x$FS$p[t]
  W <- x$W # / sum(1,abs(x$W))*(1+length(x$W))
  xn <- x$vCols
  #-((SValues-FS.pred)^2) / x$baseX^2

  #-((sum(W*s[xn])-FS.pred)^2) / x$baseX^2 * s[,x$xCols]/(1+sum(s[,x$xCols]))
  #1/(.01+(sum(W*s[,xn])-FS.pred)^2) * x$baseX * abs(s[,x$xCols]/(1+sum(s[,x$xCols])))
  #1 - abs(FS.pred - (W*s[,xn])) / x$baseX/x$numSigs #* abs((A[,x$xCols])/(1+sum(abs(A[,x$xCols]))))
  #fErr(x, rowSums(W*s[,xn])-FS.pred) * abs(A[,x$xCols])/(1+sum(abs(A[,x$xCols])))
  #kw <-
  #fErr(x, rowSums(W*s[,xn])-FS.pred)# %*% matrix(W, nrow=1)/sum(abs(W)) #* (s[,x$xCols]-(x$itv+1))/(1+sum(abs(s[,x$xCols]-(x$itv+1))))

  #2 * (fErr(x, rowSums(W*s[,xn])-FS.pred)-fErr(x, rowSums(W*s.prev[,xn])-FS.pred))*W*A[,x$xCols]

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
    #s <- getState(x, SValues)
    #s.prev <- getState(x, x$SValues.prev)
    FS <- ifelse(t>2, 2*x$FS$Si[t-1] - x$FS$Si[t-2], x$FS$Si[max(1,t-1)])
  } else {
    FS <- x$FS$Si[t]
  }

  #print(paste('SValues=',toString(SValues),'A=',toString(A),'W=',toString(W)))
  ws <- SValues * W
  norm.ws <- sqrt(sum(ws^2))
  if (is.na(norm.ws) || norm.ws==0) print(paste0('WARN: W=',toString(W),', SValues=',toString(SValues)))

  ifelse(!is.na(norm.ws) && norm.ws > 0, fErr(x, (norm.ws - FS)), -x$baseX*x$numSigs)

  #if (is.null(S.next)) {
    #return(rowSums(S[,x$vCols]) * x$FS$p[t])
    #return(-(abs(S[,x$vCols] - x$FS$Si[t])))
  #} else {
    #p.t <- ifelse(t > 1, x$FS$p[t-1], p)
    #dp <- p - p.t
    #dS.tp <- sweep(S.next[,xn], 2, as.numeric(S[xn])) * p.t
    #s <- S.next[,xn] * dp - dS.tp
    #s <- rowSums(s)
    #-((S.next[,xn]-FS.pred)^2) / x$baseX
  #}
}

#' R: reward of the State/Action pair
#' t: time t in State.next (after action A)
calculateR <- function(x, SValues, A, SValues.prev = x$SValues.prev, t = x$t) {
  #if (is.null(S.next)) {
  #  #return(rowSums(S[,x$vCols] * as.numeric(x$W)) * x$FS[t,]$p - x$FS$Si[t])
  #  #return(rowSums(S[,x$vCols]) * x$FS[t,]$p)
  #  #return(S[,x$vCols] * x$FS[t,]$p)
  #  #return(-abs(S[,x$vCols] - x$FS$Si[t]))
  #  return(0)
  #}

  #xn <- x$vCols
  #p.t <- x$FS$p[t-1]
  #dp <- x$FS$p[t] - p.t
  #v4: the to-be-expected target
  #target.s <- x$FS$Si[t]*dp - (x$FS$Si[t]-x$FS$Si[t-1])*x$FS$p[t-1]

  #dS.tp <- sweep(S.next[,xn], 2, as.numeric(S[xn])) * p.t
  #s <- S.next[,xn] * dp - dS.tp
  #s <- t(apply(s, 1, `*`, as.numeric(x$W)))
  #s <- rowSums(s) - target.s
  #v4: s <- s - target.s
  #SValues <- matrix(SValues, ncol=x$xColNum)
  W <- x$W #/ sum(1,abs(x$W))*(1+length(x$W))
  #-((sum(W*SValues)-x$FS$Si[t])^2) / x$baseX^2
  #1/(.01+(sum(W*SValues)-x$FS$Si[t])^2) * x$baseX * abs(SValues/(1+sum(SValues)))
  #1 - abs(x$FS$Si[t] - (W*SValues)) / x$baseX/x$numSigs
  #fErr(x, rowSums(W*SValues)-x$FS$Si[t]) #* abs((A[,x$xCols])/(1+sum(abs(A[,x$xCols]))))
  #print(paste('class=',class(SValues),'SValues=',toString(SValues),'A=',toString(A),'class=',class(W),'W=',toString(W)))
  FS <- x$FS$Si[t]
  ws <- SValues * W
  norm.ws <- sqrt(sum(ws^2))
  ifelse(norm.ws > 0, fErr(x, (norm.ws - FS)), -x$baseX*x$numSigs)
}

#' differenced Value-Function
#' dfl <- seq(-20,20,.1); ggplot()+geom_line(aes(dfl,dFerr(dfl,2,100)))
#'   +geom_vline(xintercept = 100, linetype='dashed')+geom_vline(xintercept = -100, linetype='dashed')
#' dfs <- seq(-2,2,.1); ggplot()+geom_line(aes(dfs,dFerr(dfs,2,100)))
#'   +geom_vline(xintercept = 1, linetype='dashed')+geom_vline(xintercept = -1, linetype='dashed')
#dFerr <- function(dFs, a, base) {
#  ifelse(abs(dFs) > 1,
#         ifelse(abs(dFs) < a*base/2,
#                -dFs,
#                (a*base)^2/(4*dFs^3)),
#         -(dFs^2))
#}
dFerr <- function(dFs, a, base) {
   ifelse(abs(dFs) < sqrt(a*base/2),
          -2*dFs/(a*base),
          -(a*base)/(2*dFs^3))
}

calculateDx <- function(x, SValues, A, SValues.prev = x$SValues.prev, t = x$t, v0 = x$baseX, base = x$baseX, a = x$numSigs, alpha = x$alpha) {
  norm.ws <- sqrt(sum((x$W*SValues)^2))
  Fs <- x$FS$Si[t]
  deltaFs <- norm.ws-Fs

  dx <- ifelse(norm.ws == 0, 1, dFerr(deltaFs, a, base)) * v0/a/base * ((x$W*SValues^2) / norm.ws)

  print(paste0('calculateDx[',x$i,']: t=',x$t,' dx=',toString(dx),', SValues=',toString(SValues),', W=',toString(x$W)))

  dx
}

#calculateDx <- function(x, S, A, S.moves = NULL, alpha = x$alpha, t = x$t) {
#  R.estimate <- estimateR(x, S, S.moves)
#  R.estimate[is.na(R.estimate)] <- -x$baseX
#  #dX <- S.moves[x$vCols]/x$xColNum * R.estimate
#  Q <- x$Q[A$aid,x$xCols]
#  Q[is.na(Q)] <-  0
#  #dx <- (dX/(1+rowSums(abs(dX))) * alpha / (1+rowSums(abs(Q)))) + Q
#  #dx <- dX + Q
#  dx <- R.estimate
#  if (deblog) { print(paste0('calculateDx[',x$i,']:')); print(Q); print(dx) }
#  dx
#}

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
  if (is.null(x) || is.null(S)) stop('possibleActions: missing input')
  if (is.numeric(S)) {
    S <- x$S[S,]
  }

  SValues.in <- SValues

  #xA <- x$A[x$A$sid == S$sid,] %>% filter(tx == x$t) %>% mutate(z=ifelse(n==0,1,n)/sum(n))
  xA <- x$A #[x$A$sid == S$sid,]

  #dX <- x$X[x$t+1,] - x$X[x$t,] # calculate Delta-X
  m <- move(x, SValues, xA)
  S.moves <- m$S.moves
  SValues = m$SValues

  good <- !is.na(S.moves$sid)
  if (sum(good) == 0) {
    print(SValues.in)
    print(xA)
    print(S.moves)
    print(SValues)
    print(paste0('no good Action in State: ',toString(S)))
    return(NULL)
  }
  xA <- xA[good,]
  S.moves <- S.moves[good,]
  SValues <- SValues[good,]
  #xA$spid <- S.moves$spid

  R <- estimateR.all(x, SValues)
  xA <- bind_cols(xA, R=R)

  #W <- x$W # rep(1,x$xColNum) # x$W
  # remember: apply-transpsed-rowSums!
  #Q <- sqrt(colSums(apply(dX,1,`*`,W)^2)) # the predicted cost/value of a move
  #Q <- rowSums(dX) # the predicted cost/value of a move
  xA$Q <- R

  #if (any(abs(Q)>x$baseX*10)) browser()
  if (nlog) {
    print(paste0('possible[',x$i,']: t=',x$t,', S=',toString(list(S))))
  }

  # by: Value increase by Delta-X, sell: Value decreased by Delta-X:
  po <- xA %>% arrange(desc(Q))
  print(paste0('possible[',x$i,']: t=',x$t,', po=',toString(head(po,1))))

  if (deblog) print(tibble(m=paste0('possible[',x$i,']:'), po), width=Inf)
  if (deblog) x$n.t <- nanoTime(paste('possible',x$numSigs), n.t)

  po
}

#' plot(2*x$epsilon*(1-e1071::sigmoid((-300:500)/100)))
getEpsilon.finDemo <- function(x, i = x$i) {
  #1/max(1/x$epsilon, i-(100*(1+(5-x$numSigs)^3)))
  2*x$epsilon*(1-e1071::sigmoid((i-300)/100))
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
    #dX <- x$X[x$t+1,] - x$X[x$t,] # calculate Delta-X
    #A$R <- sum((1 - x$X[,x$t]) * S[x$vCols])
    sv <- SValues + x$alpha * (SValues-x$SValues.prev)
    if (all(sv > 0 & sv < x$baseX)) {
      SValues <- sv
      print(paste0('ACTION[',x$i,']: t=',x$t,', ZeroAction: SValues=',toString(SValues)))
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

  R.r <- R
  x$t <- x$t + 1
  if (x$t >= x$tSteps) { # end of episode!
    terminal = T
    # v2: R <- calculateR(x, SIn) + x$SumR
    R <- x$SumR # 0 # True Online SARA # x$SumR
  } else {
    x$SumR <- R + x$SumR
    # absolute delayed reward
    #R <- R * 0
    if (F && sum(SValues) == 0) { # penalize no investment
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

  x$A$t[A$aid] <- x$t
  x$SValues.prev <- x$SValues
  x$SValues <- SValues
  x$resA <- list(A=A, R=R, S = S, SValues = SValues, R.r = R.r)
  if (!terminal & nlog) {
    print(paste0('act[',x$i,']: t=',x$t,', resA=',toString(list(x$resA)),', SIn=',toString(list(SIn))));
  }
  print(paste0('act[',x$i,']: t=',x$t,', R=',toString(R.r),', SumR=',toString(x$SumR),', SValues=',toString(x$resA$SValues),', S=',toString(list(S))))
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
  #w * a*(1/(1+exp(-k*(w+1+d)))*(1-1/(1+exp(-k*(w-1-d)))))
  k*(-.5+e1071::sigmoid(w))
}

#' prohibit parameter decay
#' @see s<-seq(-.2,.5,.01);m<-matrix(c(rep(1,length(s)),s),ncol=2);
#' md<-t(apply(m,1,deDegenerate)); ggplot()+geom_line(aes(s,md[,1]))+geom_line(aes(s,md[,2]))
deDegenerate <- function(w, k=.1, w0=.01) {
  sw <- w / max(abs(w))
  w + k * (-.5+e1071::sigmoid(sw/w0))
}

#' update W for this action A
#' @param action A as tibble/list (match is done by A$aid) or as aid
#' @param return R from last action
#' @return the environment
updateQ.finDemo <- function(x, A, R, Anext = NULL, alpha = x$alpha, gamma = x$gamma, lambda = .8, n = 8, ...) {
  if (deblog) n.t <- nanoTime(paste('updateQ',x$numSigs), x$n.t)
  if (deblog) {
    print(paste0('updateQ[',x$i,']: t=',x$t,', A=',toString(list(A)),', R=',toString(R)));
    if(!is.null(Anext)) print(paste0('updateQ[',x$i,']: t=',x$t,', Anext=',toString(Anext)));
  }
  if (is.null(x) || is.null(A) || is.null(R)) stop('updateQ: missing input')
  if (is.numeric(A)) A <- x$A[A,]
  if (is.numeric(Anext)) Anext <- x$A[Anext,]

  aid <- A$aid
  #S <- x$S[A$sid,]
  Rc <- R  #* replace(A[x$xCols],A[x$xCols]==0,1) / 2
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
        if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dQ=',toString(dQ),', Qa=',toString(Qa),', dW=',toString(dW)))
      } else {
        if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', Anext=',toString(Anext)))
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
        x$Z <- gamma*lambda * x$Z + penalize((1 - alpha*gamma*lambda*sum(x$Z*dQ)) * dQ) # /x$baseX)) * dQ) # / x$baseX)
        if (max(abs(x$Z)) > 2) {
          print(paste0('recalibrated: Z=',toString(x$Z),', Ze=',toString(Ze)))
        }

        dw.dqz <- as.numeric(sum(delta + Qa - x$Qold) * x$Z)
        dw.qz <- as.numeric(sum(Qa - x$Qold) * dQ) # / x$baseX

        We <- W
        Wa <- deDegenerate(alpha * (dw.dqz - dw.qz)) #/ 10 # not too much volatility
        W <- W + Wa
        if (max(abs(W)) > 1.5) {
          #W <- W / max(c(1,abs(as.numeric(W))))
          print(paste0('recalibrated: W=',toString(W),', We=',toString(We),', Wa=',toString(Wa)))
        }

        if (T|deblog) print(paste0('   Qold=',toString(x$Qold)))
        if (T|deblog) print(paste0('      Z=',toString(x$Z),', Ze=',toString(Ze)))
        if (T|deblog) print(paste0(' dw.dqz=',toString(dw.dqz)))
        if (T|deblog) print(paste0('  dw.qz=',toString(dw.qz)))
        if (T|deblog) print(paste0('      W=',toString(W),', We=',toString(We)))

        if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dQ=',toString(dQ),', Qa=',toString(Qa),', delta=',toString(delta),', dQ=',toString(dQ)))
        if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dQ.n=',toString(dQ.next),', Qa.n=',toString(Qa.next),', Z=',toString(x$Z),', Qold=',toString(x$Qold)))
        if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', dw.dqz=',toString(dw.dqz),'/',toString(dw.dqz/max(1,abs(dw.dqz))),', dw.qz=',toString(dw.qz),'/',toString(dw.qz/max(1,abs(dw.qz)))))
        x$Qold <- Qa.next

      } else {
        if (T|nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', invalid move: SValues=',toString(x$SValues),', Anext=',toString(Anext)))
      }
    }

    if (nlog) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', x$W=',toString(x$W),', W=',toString(W)))
  } else {
    if (nlog) if (failed) print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,' failed: x$W=',toString(x$W),', W=',toString(W)))
  }

  #if (any(abs(W)>x$baseX)) browser()

  po <- x$chosenAction
  lid <- nrow(x$log)+1
  x$log <- add_row(x$log, lid = lid, time = microbenchmark::get_nanotime())
  x$log$i[lid] <- x$i
  x$log$t[lid] <- x$t
  x$log$sid[lid] <- A$sid
  x$log$aid[lid] <- A$aid
  x$log[lid,x$logSVals.cols] <- x$SValues
  x$log[lid,x$logR.cols] <- x$resA$R.r
  x$log[lid,x$logdQ.cols] <- dQ
  x$log[lid,x$logdQn.cols] <- dQ.next
  x$log[lid,x$logQa.cols] <- Qa
  x$log[lid,x$logQan.cols] <- Qa.next
  x$log[lid,x$logDelta.cols] <- delta
  x$log[lid,x$logW.cols] <- W
  x$log[lid,x$logZ.cols] <- x$Z
  x$log$Qpo[lid] <- po$Q
  x$log$epsilon[lid] <- x$currEpsilon
  #x$log[lid,x$logQpo.cols] <- po[x$vColNames]

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
    print(paste0('updateQ[',x$i,']: t=',x$t,', aid=',aid,', SValues:',toString(x$SValues)));
    print('Invalid SValues: going SHORT!')
    #stop('Invalid SValues: going SHORT!')
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
