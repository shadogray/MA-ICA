library(tidyverse)
library(FactoMineR)

nSig <- 4
a <- 10
b <- .2
t <- 2*pi*seq(0,3,.01)
noise <- .2

X1 <- 1/2 + cos(1.0 * t)
Y1 <- .8 + sin(1/3 * t)

X2 <- -1/2 + cos(1.0 * t)
Y2 <- -.8 + sin(1/3 * t)
X <- c(X1,X2)
Y <- c(Y1,Y2)
DATA <- tibble(X,Y) %>% rowid_to_column('t')
ggplot(DATA %>% gather(key,value,-t), aes(t,value, color=key))+geom_line()+
  scale_color_discrete('Input\nSignals')+ylab('Amplitude')
ggsave('ICA_Demo_InputSignals.png')

x1 <- a/2 + a * cos(1.0 * t + runif(length(t), 0, .001)) *
  runif(length(t), 1-2*noise,1+noise)
y1 <- 2*b + b*sin(1/3 * t + runif(length(t), 0, .001)) *
  runif(length(t), 1-noise/2,1+noise/2)

x2 <- -a/2 + a*cos(1.0 * t + runif(length(t), 0, .001)) *
  runif(length(t), 1-2*noise,1+noise)
y2 <- -2*b + b*sin(1/3 * t + runif(length(t), 0, .001)) *
  runif(length(t), 1-noise/2,1+noise/2)

x <- c(x1,x2)
y <- c(y1,y2)

base.data <- tibble(x,y)
data <- tibble(.rows = nrow(base.data))
for (i in 1:(nSig-1)) {
  data <- data %>% bind_cols(tibble(x = i/nSig*x, y = (1-i/nSig)*y))
}

t.data <- data %>% rowid_to_column('t')
t.data %>% ggplot(aes(x,y))+geom_point()+scale_y_continuous(limits = c(-2,2))
ggsave('ICA_Demo_InputSignals_Noise_XY.png')

ggplot(data %>% rowid_to_column('t') %>%
         gather(key,value,-t), aes(t,value, color=key))+geom_line()+
  scale_color_discrete('Mixed\nSignals')+ylab('Amplitude')
ggsave('ICA_Demo_InputSignals_Noise_Mixed.png')



pca <- prcomp(data, center = T, scale. = T)
summary(pca)
#autoplot(pca)
FactoMineR::plot.PCA(FactoMineR::PCA(data, graph = F), choix = 'var')

data.pca <- as_tibble(as.matrix(data)%*%pca$rotation[,1:2]) %>%
  rowid_to_column('t')
ggplot(data.pca %>% gather(key,value,-t), aes(t,value, color=key))+
  geom_line()+
  scale_color_discrete('PCA\nSignals')+ylab('Amplitude')
ggsave('ICA_Demo_PCASignals.png')


ica <- fastICA::fastICA(data, n.comp = 2, row.norm = T)
summary(ica)
ggplot(as_tibble(ica$S) %>% rowid_to_column('t') %>%
         gather(key,value,-t), aes(t,value, color=key))+
  geom_line()+
  scale_color_discrete('ICA\nSignals')+
  ylab('Amplitude')#+ggtitle('Independent Components ')
ggsave('ICA_Demo_ICASignals.png')

S <- ica$S
for (i in 1:20) {
  ica <- fastICA::fastICA(data, n.comp = 2, row.norm = T)
  S <- rbind(S, ica$S)
}

ggplot(as_tibble(S) %>% rowid_to_column('t') %>%
         gather(key,value,-t), aes(t,value, color=key))+
  geom_smooth(span=.1)+
  scale_color_discrete('ICA\nSignals')+
  ylab('Amplitude')#+ggtitle('Independent Components ')
