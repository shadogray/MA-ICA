library(tidyverse)
#library(tryCatchLog)

source('finlearn.R')

nlog = F
deblog = F
if (deblog) {
  nlog = T
}
options(tibble.width = Inf, keep.source = T)
options(width = 2000, digits = 4)
#options(error = function(x=NA){traceback(1);print(x);stop(x)})

numSigs = 3
numOrig = 3

sSteps = 20
tSteps = 20
numIter = 120
ica.maxRetries = 5
deltaSumR = 2*.001
epsilon = .1
alpha = .1
gamma = .8
ctrl.types = factor('Q','SARSA')

type = 'SARSA'

envFinDemo <- finDemo(sSteps = sSteps, tSteps = tSteps, numSigs = numSigs, numOrig = numOrig,
                      epsilon = epsilon)
envFinDemo$offLimR = -1
envFinDemo <- init(envFinDemo)

stTime <- Sys.time()
tryCatchLog::tryCatchLog({
  envFinDemo <- learn(envFinDemo, numIter, type)
  print(paste('execution: ',Sys.time()-stTime,'SumR:',envFinDemo$SumR))
}, write.error.dump.file = T)

envFinDemo$SumR <- 0
result <- findWay(envFinDemo)
print(result$track, n=nrow(result$track))
print(paste('Type:',type,'Iterations:',result$myEnv$Iterations,'SumR:',result$myEnv$SumR))

