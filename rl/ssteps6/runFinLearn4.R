library(tidyverse)
#library(tryCatchLog)

source('finlearn.R')

nlog = F
deblog = F
if (deblog) {
  nlog = T
}
options(tibble.width = Inf)
options(width = 2000)
#options(error = function(x=NA){traceback(1);print(x);stop(x)})

numSigs = 4
numOrig = 6

sSteps = 40
tSteps = 20
numIter = 6000
deltaSumR = 2*.001
epsilon = .2
alpha = .1
gamma = .8
ctrl.types = factor('Q','SARSA')

type = 'SARSA'

envFinDemo <- finDemo(sSteps = sSteps, tSteps = tSteps, numSigs = numSigs, numOrig = numOrig, epsilon = epsilon, fs.pred='F')
envFinDemo$offLimR = -1
envFinDemo <- init(envFinDemo)

stTime <- Sys.time()
#tryCatchLog::tryCatchLog({
    envFinDemo <- learn(envFinDemo, numIter, type)
#}, write.error.dump.file = T)

envFinDemo$SumR <- 0
result <- findWay(envFinDemo)
print(result$track, n=nrow(result$track))
print(paste('Type:',type,'Iterations:',result$myEnv$Iterations,'SumR:',result$myEnv$SumR))
