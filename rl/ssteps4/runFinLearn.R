library(tidyverse)
options(keep.source = T)
args = commandArgs(trailingOnly=TRUE)

runs = 20
numSigs = 2

if (length(args) > 0) {
  for (idx in 1:length(args)) {
    arg = args[idx]
    if (str_detect(arg,'runs=')) runs = as.numeric(str_extract(arg, '\\d+'))
    if (str_detect(arg,'numSigs=')) numSigs = as.numeric(str_extract(arg, '\\d+'))
  }
}

for (run in 1:runs) {
  file <- paste0('runFinLearn',numSigs,'.R')
  print(paste('Running: run=',run,'numSigs=', numSigs, 'runs=', runs, 'file=',file))
  source(file = file)
}