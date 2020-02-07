
soq <- function (x) { sum(x^2) }

eval_freq2 <- function(env, rng = 5, high = 10) {
  x <- env$X
  xo <- env$X.orig
  b <- tibble(.rows = length(env$B1))
  for (n in names(env)[grep('B\\d',names(env))]) b <- bind_cols(b, env[n])
  numSigs <- env$numSigs
  numOrig <- env$numOrig
  sigs <- paste0(numSigs,numOrig)

  i <- 1
  fmbsft <- rep(0, length(env$tSteps))
  fmxosft <- rep(0, length(env$tSteps))
  fmxsft <- rep(0, length(env$tSteps))

  for (i in 1:(env$numOrig)) {
    bsft <- e1071::stft(pull(b,paste0('B',i)),win=5,inc=1)
    upper = round(ncol(bsft$values)-high):ncol(bsft$values)
    bsft.vals <- bsft$values[, upper]
    fmbsft <- fmbsft + (Re(fft(apply(bsft.vals,1,mean))))^2
    #bsft$values <- bsft$values[,upper]
    #plot(bsft)
    #title(paste('B',i))
  }
  for (i in 1:env$numSigs) {
    xosft <- e1071::stft(pull(xo,paste0('S',i)),win=5,inc=1)
    xosft.vals <- xosft$values
    upper = round(ncol(xosft$values)-high):ncol(xosft$values)
    xosft.vals <- xosft$values[, upper]
    fmxosft <- fmxosft + (Re(fft(apply(xosft.vals,1,mean))))^2
    #print(paste('S.orig',i,toString(fmxosft)))

    xsft <- e1071::stft(pull(x,paste0('S',i)),win=5,inc=1)
    upper = round(ncol(xsft$values)-high):ncol(xsft$values)
    xsft.vals <- xsft$values[, upper]
    fmxsft <- fmxsft + Re(fft(apply(xsft.vals,1,mean)))^2
    #print(paste('S',i,toString(fmxsft)))
    #xsft$values <- xsft$values[,upper]
    #plot(xsft)
    #title(paste('X',i))
  }

  len <- length(fmbsft)
  inner <- (rng+1):(len-rng)
  outer <- c(1:rng,(len-rng):len)

  ratio.b <- soq(Re(fmbsft)[inner])/soq(Re(fmbsft)[outer])
  ratio.xo <- soq(Re(fmxosft)[inner])/soq(Re(fmxosft)[outer])
  ratio.x <- soq(Re(fmxsft)[inner])/soq(Re(fmxsft)[outer])
  print(paste('Ratio-B:',ratio.b,'Ratio-X.orig:',ratio.xo,'Ratio-X:',ratio.x,'Ratio:',ratio.x/ratio.xo))
  return(list(ratio = ratio.x/ratio.xo, ratio.b = ratio.b, ratio.xo = ratio.xo, ratio.x = ratio.x))
}
