library(tidyverse)

init <- function(x) { UseMethod('init') }
getStart <- function(x) { UseMethod('getStart') }
getState <- function(x, s, ...) { UseMethod('getState') }
mapState <- function(x, X, ...) { UseMethod('mapState') }
getStateAction <- function(x, s, a, ...) { UseMethod('getStateAction') }
possibleActions <- function(x, s, ...) { UseMethod('possibleActions') }
getEpsilon <- function(x, ...) { UseMethod('getEpsilon') }
act <- function(x, s, a, ...) { UseMethod('act') }
offLimits <- function(x, s, a, ...) { UseMethod('offLimits') }
evalAction <- function(x, s, a, ...) { UseMethod('evalAction')}
predictQ <- function(x, a, ...) { UseMethod('predictQ')}
updateQ <- function(x, a, v, ...) { UseMethod('updateQ') }
test <- function(x, ...) { UseMethod('test') }

new_world <- function(..., class = character()) {
  structure(list(...), class = c(class,'world'))
}

world <- function(stateActions, ...) {
  new_world(stateActions, ...)
}

init.world <- function(x, ...) { NULL }
getStart.world <- function(x, ...) { NULL }
getState.world <- function(x, S, ...) { NULL }
mapState.world <- function(x, X, ...) { NULL }
getStateAction.world <- function(x, S, A, ...) { NULL }
possibleActions.world <- function(x, S, ...) { NULL }
getEpsilon.world <- function(x, ...) { x$epsilon }
act.world <- function(x, S, A, ...) { NULL }
offLimits.world <- function(x, S, A, ...) { NULL }
evalAction.world <- function(x, S, A, ...) { NULL }
predictQ.world <- function(x, A, ...) { NULL }
updateQ.world <- function(x, A, q, ...) { x }

##############################################################
######################### template ##########################
new_template <- function(a,b,c, ..., class = character()) {
  new_world(class = c(class, 'template'))
}
template <- function(a,b,c, ...) { new_template(a,b,c, ...) }
init.template <- function(x, ...) { NULL }
getStart.template <- function(x, ...) { NULL }
getState.template <- function(x, S, ...) { NULL }
mapState.template <- function(x, X, ...) { NULL }
getStateAction.template <- function(x, S, A, ...) { NULL }
possibleActions.template <- function(x, S, ...) { NULL }
getEpsilon.template <- function(x, ...) { x$alpha }
evalAction.template <- function(x, S, A) {}
offLimits.template <- function(x, S, A, ...) {}
act.template <- function(x, S, A, ...) {}
updateQ.template <- function(x, A, q, ...) { x }
test.template <- function(x, ...) { print('not implemented') }
######################### template ##########################
##############################################################
