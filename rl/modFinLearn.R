getEpsilon.finDemo <- function(x, i = x$i) {
  1/max(1/x$epsilon, i-500)
}

