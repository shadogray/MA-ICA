library(tidyverse)

source('finlearn.R')

nlog = T
deblog = F
if (deblog) {
  nlog = T
}
options(tibble.width = Inf)
options(width = 2000)

numSigs = 5
numOrig = 3

sSteps = 10
tSteps = 20
numIter = 100000
deltaSumR = 2*.001
epsilon = .2
alpha = .1
gamma = .8
ctrl.types = factor('Q','SARSA')

type = 'SARSA'
args <- commandArgs(trailingOnly = T)
if (length(args) > 0) type = args[1]

envFinDemo <- finDemo(sSteps = sSteps, tSteps = tSteps, numSigs = numSigs, numOrig = numOrig)
envFinDemo$offLimR = -1
envFinDemo <- init(envFinDemo)

stTime <- Sys.time()
envFinDemo <- learn(envFinDemo, numIter, type)
print(paste('execution: ',Sys.time()-stTime,'SumR:',envFinDemo$SumR))

envFinDemo$SumR <- 0
result <- findWay(envFinDemo)
print(result$track, n=nrow(result$track))
print(paste('Type:',type,'Iterations:',result$myEnv$Iterations,'SumR:',result$myEnv$SumR))

