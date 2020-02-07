library(tidyverse)
#library(signal)

source('utils.R')

load('ssteps55.rdata')
load('ssteps33.rdata')

x <- result$myEnv$X
xo <- result$myEnv$X.orig
b <- tibble(.rows = length(result$myEnv$B1))
for (n in names(result$myEnv)[grep('B\\d',names(result$myEnv))]) b <- bind_cols(b, result$myEnv[n])
sigs <- paste0(result$myEnv$numSigs,result$myEnv$numOrig)

png(paste0('base_stft_',sigs,'.png'))
par(mfrow=c(2,2))
for (i in 1:4) { plot(e1071::stft(pull(b,paste0('B',i)),win=4,inc=1)); title(main=paste0('STFT: B',i)) }
dev.off()
png(paste0('signals_stft_',sigs,'.png'))
par(mfrow=c(2,2))
for (i in 1:4) { plot(e1071::stft(pull(x,paste0('S',i)),win=4,inc=1)); title(main=paste0('STFT: S',i)) }
dev.off()

i <- 1
bsft <- e1071::stft(pull(b,paste0('B',i)),win=5,inc=1)
bsft.vals <- bsft$values*0
xosft <- e1071::stft(pull(xo,paste0('S',i)),win=5,inc=1)
xosft.vals <- xosft$values*0
xsft <- e1071::stft(pull(x,paste0('S',i)),win=5,inc=1)
xsft.vals <- xsft$values*0

for (i in 1:(result$myEnv$numOrig)) {
  bsft.vals <- bsft.vals + e1071::stft(pull(b,paste0('B',i)),win=5,inc=1)$values
}
for (i in 1:result$myEnv$numSigs) {
  xosft.vals <- xosft.vals + e1071::stft(pull(xo,paste0('S',i)),win=5,inc=1)$values
  xsft.vals <- xsft.vals + e1071::stft(pull(x,paste0('S',i)),win=5,inc=1)$values
}

mbsft <- apply(bsft.vals,1,mean)
mxosft <- apply(xosft.vals,1,mean)
mxsft <- apply(xsft.vals,1,mean)

fmbsft <- fft(apply(bsft$values,1,mean))
fmxosft <- fft(apply(xosft$values,1,mean))
fmxsft <- fft(apply(xsft$values,1,mean))
#plot(Mod(fmbsft)[4:12])
plot(Mod(fmxosft)[4:12])
plot(Mod(fmxsft)[4:12])

rng <- 4
len <- length(fmbsft)
inner <- (rng+1):(len-rng)
outer <- c(1:rng,(len-rng):len)

ratio.b <- soq(Mod(fmbsft)[inner])/soq(Mod(fmbsft)[outer])
ratio.xo <- soq(Mod(fmxosft)[inner])/soq(Mod(fmxosft)[outer])
ratio.x <- soq(Mod(fmxsft)[inner])/soq(Mod(fmxsft)[outer])
print(paste('Ratio-B:',ratio.b,'Ratio-X.orig:',ratio.xo,'Ratio-X:',ratio.x))


eval_freq <- function(res) {
  x <- res$myEnv$X
  xo <- res$myEnv$X.orig
  b <- tibble(.rows = length(res$myEnv$B1))
  for (n in names(res$myEnv)[grep('B\\d',names(res$myEnv))]) b <- bind_cols(b, res$myEnv[n])
  numSigs <- res$myEnv$numSigs
  numOrig <- res$myEnv$numOrig
  sigs <- paste0(numSigs,numOrig)

  png(paste0('base_stft_',sigs,'.png'))
  par(mfrow=c(2,2))
  for (i in 1:min(4,(numOrig))) { plot(e1071::stft(pull(b,paste0('B',i)),win=4,inc=1)); title(main=paste0('STFT: B',i)) }
  dev.off()
  png(paste0('signals_stft_',sigs,'.png'))
  par(mfrow=c(2,2))
  for (i in 1:min(4,numSigs)) { plot(e1071::stft(pull(x,paste0('S',i)),win=4,inc=1)); title(main=paste0('STFT: S',i)) }
  dev.off()

  i <- 1
  bsft <- e1071::stft(pull(b,paste0('B',i)),win=5,inc=1)
  bsft.vals <- bsft$values*0
  xosft <- e1071::stft(pull(xo,paste0('S',i)),win=5,inc=1)
  xosft.vals <- xosft$values*0
  xsft <- e1071::stft(pull(x,paste0('S',i)),win=5,inc=1)
  xsft.vals <- xsft$values*0

  for (i in 1:(res$myEnv$numOrig)) {
    bsft.vals <- bsft.vals + e1071::stft(pull(b,paste0('B',i)),win=5,inc=1)$values
  }
  for (i in 1:res$myEnv$numSigs) {
    xosft.vals <- xosft.vals + e1071::stft(pull(xo,paste0('S',i)),win=5,inc=1)$values
    xsft.vals <- xsft.vals + e1071::stft(pull(x,paste0('S',i)),win=5,inc=1)$values
  }

  mbsft <- apply(bsft.vals,1,mean)
  mxosft <- apply(xosft.vals,1,mean)
  mxsft <- apply(xsft.vals,1,mean)

  fmbsft <- fft(apply(bsft.vals,1,mean))
  fmxosft <- fft(apply(xosft.vals,1,mean))
  fmxsft <- fft(apply(xsft.vals,1,mean))
  #plot(Mod(fmbsft)[4:12])
  #plot(Mod(fmxosft)[4:12])
  #plot(Mod(fmxsft)[4:12])

  rng <- 4
  len <- length(fmbsft)
  inner <- (rng+1):(len-rng)
  outer <- c(1:rng,(len-rng):len)

  ratio.b <- soq(Mod(fmbsft)[inner])/soq(Mod(fmbsft)[outer])
  ratio.xo <- soq(Mod(fmxosft)[inner])/soq(Mod(fmxosft)[outer])
  ratio.x <- soq(Mod(fmxsft)[inner])/soq(Mod(fmxsft)[outer])
  print(paste('Ratio-B:',ratio.b,'Ratio-X.orig:',ratio.xo,'Ratio-X:',ratio.x))
}

load('ssteps55.rdata')
load('ssteps33.rdata')
eval_freq(result)


load('ssteps42.rdata')
eval_freq2(result$myEnv)
load('ssteps43.rdata')
eval_freq2(result$myEnv)
