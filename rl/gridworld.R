library(tidyverse)
library(ReinforcementLearning)
library(MDPtoolbox)

load('etf.syms.rdata')

N = 1000
# Sampling data (1000 grid sequences)
data <- sampleGridSequence(N)
# Setting reinforcement learning parameters
control <- list(alpha = 0.1, gamma = 0.1, epsilon = 0.1)

states <- c("s1", "s2", "s3", "s4")
actions <- c("up", "down", "left", "right")

model <- NULL
actionSelection = 'random' # 'epsilon-greedy'

for (i in 1:20) {
  data <- sampleExperience(N = N,
                           states = states,
                           actions = actions,
                           env = gridworldEnvironment,
                           actionSelection = actionSelection,
                           control = control,
                           model = model)

  # Performing reinforcement learning
  model <- ReinforcementLearning(data,
                                 s = "State", a = "Action", r = "Reward",
                                 s_new = "NextState",
                                 control = control,
                                 model = model)

  # prepare for next iteration:
  actionSelection = 'epsilon-greedy'
  # Printing model
  print(model)
}

# Plotting learning curve
plot(model)
