library(tidyverse)

source('worlds.R')
source('learn.R')

deblog = F
numIter = 700
deltaSumR = 2*.001
epsilon = .1
alpha = .1
gamma = 1
ctrl.types = factor('Q','SARSA')

type = 'Q' # 'Q' #'SARSA'
args <- commandArgs(trailingOnly = T)
if (length(args) > 0) type = args[1]

envTheCliff <- theCliff(1:12,1:3)
envTheCliff$offLimR = -1
envTheCliff <- init(envTheCliff)

stTime <- Sys.time()
envTheCliff <- learn(envTheCliff)
print(paste('execution: ',Sys.time()-stTime,'SumR:',envTheCliff$SumR))

envTheCliff$SumR <- 0
result <- findWay(envTheCliff)
print(result$track, n=nrow(result$track))
print(paste('Type:',type,'Iterations:',result$myEnv$Iterations,'SumR:',result$myEnv$SumR))

