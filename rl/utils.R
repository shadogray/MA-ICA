library(tidyverse)

soq <- function (x) { sum(x^2) }
debug = T

eval_freq2 <- function(myEnv, rng = 5, high = 10, ignore = 5) {
  x <- myEnv$X %>% select(-matches('SC|S.\\.i'))
  xo <- myEnv$X.orig
  b <- tibble(.rows = length(myEnv$B1))
  for (n in names(myEnv)[grep('B\\d',names(myEnv))]) b <- bind_cols(b, myEnv[n])
  numSigs <- myEnv$numSigs
  numOrig <- myEnv$numOrig
  sigs <- paste0(numSigs,numOrig)
  if (any(is.na(x))) {
    print("x isNa:"); print(x)
  }
  if (any(is.na(xo))) {
    print("xo isNa:"); print(xo)
  }
  if (any(is.na(b))) {
    print("b isNa:"); print(b)
  }

  i <- 1
  fmbsft <- rep(0, length(myEnv$tSteps))
  fmxosft <- rep(0, length(myEnv$tSteps))
  fmxsft <- rep(0, length(myEnv$tSteps))

  for (i in 1:(myEnv$numOrig)) {
    bsft <- e1071::stft(pull(b,paste0('B',i)),win=5,inc=1)
    upper = round(ncol(bsft$values)-high-ignore):(ncol(bsft$values)-ignore)
    bsft.vals <- bsft$values[, upper]
    if (any(is.na(bsft.vals))) { print("bsft.vals isNa:"); print(bsft.vals) }
    fmbsft <- fmbsft + (Mod(fft(apply(bsft.vals,1,mean))))
    if (any(is.na(fmbsft))) { print("fmsft isNa:"); print(fmbsft) }
    #fmbsft <- fmbsft + var(scale(apply(bsft.vals,1,mean), scale=F))
    if (debug) {
      bsft$values <- bsft$values[,upper]
      plot(bsft)
      title(paste('B-Upper',i))
    }
  }
  for (ni in names(xo)) {
    xosft <- e1071::stft(pull(xo,ni),win=5,inc=1)
    xosft.vals <- xosft$values
    if (any(is.na(xosft.vals))) { print("xosft.vals isNa:"); print(xosft.vals) }
    upper = round(ncol(xosft$values)-high-ignore):(ncol(xosft$values)-ignore)
    xosft.vals <- xosft$values[, upper]
    fmxosft <- fmxosft + (Mod(fft(apply(xosft.vals,1,mean))))
    if (any(is.na(fmxosft))) { print("fmxosft isNa:"); print(fmxosft) }
    #fmxosft <- fmxosft + var(scale(apply(xosft.vals,1,mean), scale=F))
    if (debug) {
      print(paste('S.orig',ni,toString(fmxosft)))
      plot(xosft)
      title(paste('S-Orig-Upper',ni))
    }
  }
  for (i in 1:myEnv$numSigs) {
    xsft <- e1071::stft(pull(x,paste0('S',i)),win=5,inc=1)
    upper = round(ncol(xsft$values)-high-ignore):(ncol(xsft$values)-ignore)
    xsft.vals <- xsft$values[, upper]
    fmxsft <- fmxsft + Mod(fft(apply(xsft.vals,1,mean)))
    #fmxsft <- fmxsft + var(scale(apply(xsft.vals,1,mean), scale=F))
    if (debug) {
      print(paste('S',i,toString(fmxsft)))
      plot(xsft)
      title(paste('X',i))
      #Sys.sleep(2)
      xsft$values <- xsft$values[,upper]
      plot(xsft)
      title(paste('X-UPPER',i))
    }
  }

  len <- length(fmbsft)
  lower <- 1:(rng)
  upper <- (rng+1):len

  ratio.b <- sum(Mod(fmbsft[upper]))/sum(Mod(fmbsft[lower]))
  ratio.xo <- sum(Mod(fmxosft[upper]))/sum(Mod(fmxosft[lower]))
  ratio.x <- sum(Mod(fmxsft[upper]))/sum(Mod(fmxsft[lower]))
  #ratio.b <- fmbsft/myEnv$numOrig
  #ratio.xo <- fmxosft/myEnv$numBase
  #ratio.x <- fmxsft/myEnv$numSigs
  print(paste('Ratio-B:',ratio.b,'Ratio-X.orig:',ratio.xo,'Ratio-X:',ratio.x,'Ratio:',ratio.x/ratio.xo))
  return(list(ratio = ratio.x/ratio.xo, ratio.b = ratio.b, ratio.xo = ratio.xo, ratio.x = ratio.x))
}
