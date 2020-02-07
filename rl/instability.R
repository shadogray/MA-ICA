library(tidyverse)

load('ssteps33.rdata')
log <- result$myEnv$log

w.max <- tibble(.rows=1)
W.diff <- tibble(.rows = length(log$W.S1)-1)
for (n in grep('W\\.S.',names(log), value = T)) {
  W.diff[n] <- diff(pull(log,n))
  print(paste('Diff: ', n, 'val:',max(abs(pull(W.diff,n))),'idx:',which.max(abs(pull(W.diff,n)))))
  w.max[n] <- which.max(abs(pull(W.diff,n)))
}

wmx <- unlist(w.max)
wmx <- c(wmx-1, wmx, wmx+1) %>% sort()
log[unlist(wmx),] %>% select(lid,time,i,t,sid,epsilon,matches('dQ\\.|W'))

ggplot(W.diff %>% gather(key,val)) + geom_histogram(aes(val))+facet_wrap(~key)
