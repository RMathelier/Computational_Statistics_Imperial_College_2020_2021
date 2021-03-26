###   MATH97125 - Computational Statistics, Coursework final, question 1  ###

### code of Robin Mathelier ###


setwd("C:/Users/robin/Dropbox/Applications/Overleaf/Coursework Final Computational") # Ã  personnaliser
getwd()
.libPaths("C:/Users/robin/Documents/R/win-library/3.6")
rm(list=objects())

library(xtable)
library(latex2exp)
library(expm)
library(coda)
library(ggplot2)
library(MASS)
library('pracma')

set.seed(42)

## Question 1 ##

f = function(z) exp(-z[1]^2/100 - (z[2]+3/100*z[1]^2 -3)^2)

log_f = function(z) -z[1]^2/100 - (z[2]+3/100*z[1]^2 -3)^2

f = function(x,y) exp(-x^2/100 - (y+3/100*x^2 -3)^2)
IntY <- function(x) sapply(x, function(b) integrate(function(y) f(b,y),-Inf,Inf)$value)
integrate(function(x) x*IntY(x), -Inf, Inf)$value
IntX <- function(y) sapply(y, function(b) integrate(function(x) f(x,b),-Inf,Inf)$value)
k = integrate(function(y) IntX(y), -Inf, Inf)$value
integrate(function(y) y*IntX(y), -Inf, Inf)$value/k

mvdnorm = function(vec,mu,corr_matrix) {
  x=vec[1]
  y=vec[2]
  x_mu = as.matrix(c(x,y)-mu)
  return( (1/(2*pi*sqrt(det(corr_matrix))))*
  exp((-1/2)*t(x_mu)%*%solve(corr_matrix)%*%x_mu) )
}

create_chains = function(init, corr_matrix){
  q <- function(vec_1,vec_2) mvdnorm(vec_1,mu=vec_2,corr_matrix=corr_matrix)
  MHpost <- function(n=1e4,init){
    xall <- matrix(NA,ncol=n, nrow=2)
    x <- init
    xall[,1] <- x
    for (t in seq(2,n)){
      y <- mvrnorm(1,mu=x,Sigma=corr_matrix)
      if (runif(1)<=(exp(log_f(y)-log_f(x))*q(y,x)/q(x,y))){
        x<-y
      }
      xall[,t] <- x
    }
    return(xall)
  }
  return(MHpost(init=init))
}

theta = c(0,3)

create_corr_mat = function(sx,sy,ro) rbind(c(sx^2,ro*sx*sy),
                                           c(ro*sx*sy,sy^2))

corr_matrix = create_corr_mat(sx=5,sy=1,ro=0.5)

chains = create_chains(init=theta,corr_matrix=corr_matrix)
dim(chains)

chain_x = chains[1,]
traceplot(mcmc(chain_x))

chain_y = chains[2,]
traceplot(mcmc(chain_y))

effectiveSize(chain_x)
effectiveSize(chain_y)

acf(chain_x,lag.max=100)
acf(chain_y)

list_sx = linspace(7,9,n=4)
list_sy = linspace(1,2,n=4)
list_ro = c(0.005,0.01,0.05,0.1,0.15)

set.seed(42)

ESS_x = 0
ESS_y = 0
for (sx in list_sx){
  for (sy in list_sy){
    for (ro in list_ro){
      corr_matrix = create_corr_mat(sx,sy,ro)
      if(eigen(corr_matrix)$values[1]>0 & eigen(corr_matrix)$values[2]>0){
        chains = create_chains(init=theta,corr_matrix=corr_matrix)
        chain_x = chains[1,]
        chain_y = chains[2,]
        new_ESS_x = effectiveSize(mcmc(chain_x))
        new_ESS_y = effectiveSize(mcmc(chain_y))
        if (new_ESS_x > ESS_x & new_ESS_y > ESS_y){
          ESS_x = new_ESS_x
          ESS_y = new_ESS_y
          sx_chosen = sx
          sy_chosen = sy
          ro_chosen = ro
          chain_x_chosen = chain_x
          chain_y_chosen = chain_y
        }
      }
    }
  }
}


cat('sx:',sx_chosen,'sy:',sy_chosen,'rho',ro_chosen)
chain_x = chain_x_chosen
chain_y = chain_y_chosen

acf(chain_x,lag.max=100)
acf(chain_y,lag.max=100)

effectiveSize((chain_x))
effectiveSize((chain_y))

traceplot(mcmc(chain_x))
traceplot(mcmc(chain_y))

plot(chain_x,chain_y)

chain_x_ts = ts(chain_x)
chain_y_ts = ts(chain_y)

acf(ts.union(chain_x_ts,chain_y_ts))

create_chains_def = function(init, corr_matrix,n_iter=1e4){
  q <- function(vec_1,vec_2) mvdnorm(vec_1,mu=vec_2,corr_matrix=corr_matrix)
  MHpost <- function(n=n_iter,init){
    xall <- matrix(NA,ncol=n, nrow=2)
    x <- init
    xall[,1] <- x
    for (t in seq(2,n)){
      y <- mvrnorm(1,mu=x,Sigma=corr_matrix)
      if (runif(1)<=(exp(log_f(y)-log_f(x)))){
        x<-y
      }
      xall[,t] <- x
    }
    return(xall)
  }
  return(MHpost(init=init))
}

set.seed(42)
theta=c(0,3)

corr_matrix = create_corr_mat(sx=sx_chosen,sy=sy_chosen,ro=ro_chosen)
chains = create_chains_def(init=theta,corr_matrix,n_iter=1e5)

chain_x = chains[1,]
chain_y = chains[2,]

acf(chain_x,lag.max=300)
acf(chain_y,lag.max=300)

png(filename='acf_x.png')
acf(chain_x,lag.max=300)
dev.off()

png(filename='acf_y.png')
acf(chain_y,lag.max=300)
dev.off()

effectiveSize((chain_x))
effectiveSize((chain_y))

plot(chain_x,chain_y,
     col=rgb(red=0,green=0,blue=0.4,alpha=0.1),
     xlim=c(-25,25)
     ,xlab = 'X',ylab='Y',main = 'Simulated joint distribution of (X,Y)')

png(filename='plot_x_y.png')
plot(chain_x,chain_y,
     col=rgb(red=0,green=0,blue=0.4,alpha=0.1),
     xlim=c(-25,25),ylim=c(-12,7),
     xlab = 'X',ylab='Y',main = 'Simulated joint distribution of (X,Y)')
dev.off()

traceplot(mcmc(chain_x))
traceplot(mcmc(chain_y))

list_theta = list(cbind(-20,-20),
                  cbind(-10,-10),
                  cbind(0,0),
                  cbind(10,10),
                  cbind(20,20))

set.seed(42)

chains_x=c()
chains_y=c()

for (i in 1:length(list_theta)){
  theta = as.vector(list_theta[[i]])
  chain = create_chains_def(init=theta,corr_matrix,n_iter = 5e4)
  chain_x=chain[1,]
  chain_y=chain[2,]
  chains_x=rbind(chains_x,chain_x)
  chains_y=rbind(chains_y,chain_y)
}

set.seed(42)

chains_x_mcmc <- lapply(1:length(list_theta), 
                    function(i) mcmc(chains_x[i,]))

chains_y_mcmc <- lapply(1:length(list_theta), 
                    function(i) mcmc(chains_y[i,]))

png(filename='traceplot_x_init_values.png')
traceplot(chains_x_mcmc,main='Traceplot of the Markov chain simulated for X',
          ylab = 'X_t',xlab = 'iteration t')
abline(v=2000,lwd=3,lty=2)
dev.off()

png(filename='traceplot_y_init_values.png')
traceplot(chains_y_mcmc,main='Traceplot of the Markov chain simulated for Y',
          ylab = 'Y_t',xlab = 'iteration t')
abline(v=2000,lwd=3,lty=2)
dev.off()

gelman.diag(chains_x_mcmc)
gelman.diag(chains_y_mcmc)

### Question 2

chains_x_thin = sapply(seq(5), function(i) mean(chains_x[i,seq(1000,length(chains_x[i,]),100)]))
mean(chains_x_thin)

means_mc_x <- sapply(seq(5), function(i) mean(chains_x[i,-seq(2000)]))
mean(means_mc_x)

means_mc_y <- sapply(seq(5), function(i) mean(chains_y[i,-seq(2000)]))
mean(means_mc_y) 

chains_x_bi = chains_x[,-seq(2000)]
b = 100
n = length(chains_x_bi)
mub_x = diff(cumsum(chains_x_bi)[(0:(n/b))*b])/b
acf(mub_x)
mean(mub_x)

chains_y_bi = chains_y[,-seq(2000)]
b = 100
n = length(chains_y_bi)
mub_y = diff(cumsum(chains_y_bi)[(0:(n/b))*b])/b
acf(mub_y)
mean(mub_y)

### Question 3

d=2
H=200
c=2.4/sqrt(2)


create_chains_haario = function(init,
                                H=200,U=200,
                                c=2.4/sqrt(2),
                                n_iter=1e4,
                                corr_matrix=diag(c(1,1))){
  MHpost <- function(init){
    xall <- matrix(NA,ncol=n_iter, nrow=2)
    x <- init
    xall[,1] <- x
    for (t in seq(2,H)){
      y <- mvrnorm(1,mu=x,Sigma=corr_matrix)
      if (runif(1)<=(exp(log_f(y)-log_f(x)))){
        x<-y
      }
      xall[,t] <- x
    }
    return(xall)}
    x=init
    xall = MHpost(init=init)
    K = t(xall[,1:H])
    exp_K = sapply(1:2,function(j) mean(K[,j]))
    tilde_K = sapply(1:2,function(j) K[,j]-exp_K[j])
    for(t in (H+1):n_iter){
      if(t%%U==0){
        K = t(xall[,(t-H):(t-1)])
        exp_K = sapply(1:2,function(j) mean(K[,j]))
        tilde_K = sapply(1:2,function(j) K[,j]-exp_K[j])
      }
      y <- x+(c/sqrt(H-1))*t(tilde_K)%*%cbind(rnorm(H))
      if (runif(1)<=(exp(log_f(y)-log_f(x)))){x<-y}
      xall[,t] <- x
      }
  return(xall)
    }

set.seed(44)

chains=create_chains_haario(init=c(0,3),n_iter=1e5)

x_haario = chains[1,]
y_haario = chains[2,]

png(filename='traceplot_x_haario.png')
traceplot(mcmc(x_haario),main='Traceplot of the Markov chain simulated for X',
          ylab = 'X_t',xlab = 'iteration t')
dev.off()

png(filename='traceplot_y_haario.png')
traceplot(mcmc(y_haario),main='Traceplot of the Markov chain simulated for Y',
          ylab = 'X_t',xlab = 'iteration t')
dev.off()

plot(x_haario,y_haario,
     col=rgb(red=0,green=0,blue=0.4,alpha=0.3),
     xlim=c(-25,25)
     ,xlab = 'X',ylab='Y',main = 'Simulated joint distribution of (X,Y)')

png(filename='plot_x_y_haario.png')
plot(x_haario,y_haario,
     col=rgb(red=0,green=0,blue=0.4,alpha=0.3),
     xlim=c(-25,25),ylim=c(-12,7),
     xlab = 'X',ylab='Y',main = 'Simulated joint distribution of (X,Y)')
dev.off()

acf(x_haario,lag.max=300)
acf(y_haario,lag.max=300)

png(filename='acf_x_haario.png')
acf(x_haario,lag.max=300)
dev.off()

png(filename='acf_y_haario.png')
acf(y_haario,lag.max=300)
dev.off()

x_haario_ts = ts(x_haario)
y_haario_ts = ts(y_haario)

acf(ts.union(x_haario_ts,y_haario_ts))

x_haario_thin = x_haario[seq(500,length(x_haario),10)]
mean(x_haario_thin)

y_haario_thin = y_haario[seq(500,length(x_haario),10)]
mean(y_haario_thin)

mean(x_haario[-seq(1,5000)])
mean(y_haario[-seq(1,5000)])

effectiveSize(x_haario)
effectiveSize(y_haario)

