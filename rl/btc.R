library(tidyverse)
library(lubridate)
library(xts)

btc <- read_csv('~/Downloads/Bitfinex_BTCUSD_1h.csv') %>%
  mutate(Date = ymd_h(Date))
btc.x <- xts(btc, order.by = btc$Date)

ggplot(btc, aes(Date,Close))+geom_line()+geom_smooth()

btc.f.c <- fft(btc$Close)
btc.f.v <- fft(btc$`Volume BTC`)
btc.f <- tibble(C.r = Re(btc.f.c), C.i = Im(btc.f.c), C = btc.f.c,
                V.r = Re(btc.f.v), V.i = Im(btc.f.v), V = btc.f.v) %>%
  rowid_to_column('f') %>% mutate(f = f/nrow(btc.f)*12)

ggplot(btc.f)+
  geom_path(aes(f,C.r), color='red')+
  geom_path(aes(f,C.i), color='green')+
  scale_y_log10()+scale_x_log10()

ggplot(btc.f)+
  geom_smooth(aes(f,C.r), color='red', span=.001)+
  geom_smooth(aes(f,C.i), color='green', span=.001)+
  scale_y_log10()+scale_x_continuous()

btc.periods <- xts::periodicity(btc.x)
