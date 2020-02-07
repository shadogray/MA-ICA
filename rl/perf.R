compiler::enableJIT(0)
##### Functions #####
 
is.compile <- function(func)
{
	# this function lets us know if a function has been byte-coded or not
	#If you have a better idea for how to do this - please let me know...
    if(class(func) != "function") stop("You need to enter a function")
    last_2_lines <- tail(capture.output(func),2)
    any(grepl("bytecode:", last_2_lines)) # returns TRUE if it finds the text "bytecode:" in any of the last two lines of the function's print
}
 
# old R version of lapply
slow_func <- function(X, FUN, ...) {
   FUN <- match.fun(FUN)
   if (!is.list(X))
    X <- as.list(X)
   rval <- vector("list", length(X))
   for(i in seq(along = X))
    rval[i] <- list(FUN(X[[i]], ...))
   names(rval) <- names(X)          # keep `names' !
   return(rval)
}

f <- function(n=100) { for (i in 1:n) slow_func(1:100, is.null) }

microbenchmark::microbenchmark(f(), times=1000)

compiler::enableJIT(3)
microbenchmark::microbenchmark(f(), times=1000)
