#library(curl)
#library(jsonlite)

URL <- 'https://query2.finance.yahoo.com/v1/finance/screener'
types <- list(
  INDEX = c("Germany","United States","United Kingdom"),
  ETF = c("Large Growth", "Bear Market", "Commodities Agriculture",
                "Convertibles", "Corporate Bond", "Equity Energy", "Europe Stock", "Global Real Estate",
                "High Yield Bond", "Multicurrency", "Multisector Bond", "Multialternative",
                "Real Estate", "Preferred Stock", "World Bond", "World Stock")
)

getETFs <- function(type = 'ETF', category = 'Global Real Estate') {
  data <- list(size = 240, offset = 0,
               sortField = "fundnetassets", sortType = 'DESC',
               quoteType=type,
               query=list(
                 operator='or',
                 operands=tibble(operator='EQ', operands=c(list(c('categoryname',category))))
                ),
               userId = '', userIdType = 'guid'
  )
  json <- jsonlite::toJSON(data, auto_unbox = T)
  #json
  res <- httr::POST(URL, body = json)
  rdat <- jsonlite::fromJSON(httr::content(res, as='text'))
  quotes <- as_tibble(rdat$finance$result$quotes[[1]]) %>% mutate(category=category)
  #quotes$symbol
  return (list(type = type, category = category, quotes = quotes))
}

etf.quotes <- NULL
for (cat in types$ETF) {
  qres <- getETFs(category = cat)
  print(paste("retrieved quotes:",cat,'rows:',nrow(qres$quotes)))
  if (is.null(etf.quotes)) etf.quotes <- qres$quotes
  else etf.quotes <- bind_rows(etf.quotes, qres$quotes)
  save(etf.quotes, file='etf.quotes.rdata.tmp')
  Sys.sleep(5)
}
save(etf.quotes, file='etf.quotes.rdata')

etf.syms <- list()
for (sym in etf.quotes$symbol) {
  data <- quantmod::getSymbols(sym, env=NULL, auto.assign = F)
  print(paste("retrieved symbols:",sym,'nrow:',nrow(obj)))
  etf.syms[[sym]] <- data
  save(etf.syms, file = paste0('etf.syms.rdata'))
  Sys.sleep(5)
}
