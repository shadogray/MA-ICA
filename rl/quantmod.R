library(quantmod)
library(reticulate)

getSymbols('GOOGL')
chartSeries(GOOGL)
addMACD()
addBBands()

getSymbols('GOGL')
chartSeries(GOGL)
addMACD()
addBBands()

getSymbols('AAPL')
chartSeries(AAPL)
addMACD()
addBBands()

getSymbols('XPT/USD',src='oanda')
chartSeries(XPTUSD)
addMACD()
addBBands()

fin <- reticulate::import('yahoo_fin')
fin$stock_info$get_live_price('aapl')

load('etf.quotes.rdata')

