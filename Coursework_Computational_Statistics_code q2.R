###   MATH97125 - Computational Statistics, Coursework final, question 2 ###

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
library(pracma)
library('purrr')


set.seed(42)

## part a

a = 4

g = function(theta) (4/10)*exp(-(a-theta)^2) + (6/10)*exp(-(30-theta)^2)

plot(g,from=0,to=40,main='plot of g(x)',xlab='x',ylab='g(x)')

mh_rw = function(n=1e4,sigma=1,init=0){
  Xres = rep(NA,n)
  X = init
  for (i in 1:n){
    Y = X+rnorm(1,mean=0,sd=sigma)
    if (runif(1) <= g(Y)/g(X)){X = Y}
    Xres[i] = X
  }
  Xres
}

chain = mh_rw()

inits = c(0,5,25,35) ## starting values
inits

# sigma = 1

ess_s1 = lapply(1:length(inits), 
                function(i) effectiveSize(mh_rw(init=inits[i], sigma=1)))
ess_s1

chains_1 = lapply(1:length(inits), 
                  function(i) mcmc(mh_rw(init=inits[i], sigma=1)))
chains_1 = mcmc.list(chains_1)
gel_1 = gelman.diag(chains_1)[[1]][1]
gel_1_CI = gelman.diag(chains_1)[[1]][2]

chain_1 = mh_rw(init=inits[1], sigma=1)
n_diff_1 = 1
for (i in 1:(length(chain_1)-1)){
  if(chain_1[i+1] != chain_1[i]){n_diff_1= n_diff_1+1}
}
accept_rate_1 = n_diff_1/length(chain_1)
accept_rate_1

chains = lapply(1:length(inits), function(i) mcmc(mh_rw(init=inits[i], sigma=1)))
png(filename='traceplot_s1.png')
traceplot(chains,main='Traceplot of the Markov chain simulated for theta - sigma = 1',
          ylab='theta_t',xlab='iteration t')
dev.off()

# sigma = 2

ess_s2 = lapply(1:length(inits), 
                function(i) effectiveSize(mh_rw(init=inits[i], sigma=2)))
ess_s2

chains_2 = lapply(1:length(inits), 
                  function(i) mcmc(mh_rw(init=inits[i], sigma=2)))
chains_2 = mcmc.list(chains)
gel_2 = gelman.diag(chains_2)[[1]][1]
gel_2_CI = gelman.diag(chains_2)[[1]][2]

chain_2 = mh_rw(init=inits[1], sigma=2)
n_diff_2 = 1
for (i in 1:(length(chain_2)-1)){
  if(chain_2[i+1] != chain_2[i]){n_diff_2= n_diff_2+1}
}
accept_rate_2 = n_diff_2/length(chain_2)
accept_rate_2

# sigma = 5

ess_s3 = lapply(1:length(inits), 
                function(i) effectiveSize(mh_rw(init=inits[i], sigma=5)))
ess_s3

chains_3 = lapply(1:length(inits), 
                  function(i) mcmc(mh_rw(init=inits[i], sigma=5)))
chains_3 = mcmc.list(chains_3)
gel_3 = gelman.diag(chains_3)[[1]][1]
gel_3_CI = gelman.diag(chains_3)[[1]][2]

chain_3 = mh_rw(init=inits[1], sigma=5)
traceplot(mcmc(chain_3))
n_diff_3 = 1
for (i in 1:(length(chain_3)-1)){
  if(chain_3[i+1] != chain_3[i]){n_diff_3= n_diff_3+1}
}
accept_rate_3 = n_diff_3/length(chain_3)
accept_rate_3

# sigma = 10

ess_s4 = lapply(1:length(inits), 
                function(i) effectiveSize(mh_rw(n=1e5,init=inits[i], sigma=10)))
ess_s4

chains_4 = lapply(1:length(inits), 
                  function(i) mcmc(mh_rw(init=inits[i], sigma=10)))
chains_4 = mcmc.list(chains_4)
gel_4 = gelman.diag(chains_4)[[1]][1]
gel_4_CI = gelman.diag(chains_4)[[1]][2]

chain_4 = mh_rw(init=inits[1], sigma=10)
n_diff_4 = 1
for (i in 1:(length(chain_4)-1)){
  if(chain_4[i+1] != chain_4[i]){n_diff_4= n_diff_4+1}
}
accept_rate_4 = n_diff_4/length(chain_4)
accept_rate_4

chains = lapply(1:length(inits), function(i) mcmc(mh_rw(init=inits[i], sigma=10)))
png(filename = 'traceplot_s2.png')
traceplot(chains,main='Traceplot of the Markov chain simulated for theta - sigma=10',
          ylab='theta_t',xlab='iteration t')
dev.off()

# sigma = 20

ess_s5 = lapply(1:length(inits), 
                function(i) effectiveSize(mh_rw(init=inits[i], sigma=20)))
ess_s5

chains_5 = lapply(1:length(inits), 
                  function(i) mcmc(mh_rw(init=inits[i], sigma=20)))
chains_5 = mcmc.list(chains_5)
gel_5 = gelman.diag(chains_5)[[1]][1]
gel_5_CI = gelman.diag(chains_5)[[1]][2]

chain_5 = mh_rw(init=inits[1], sigma=20)
n_diff_5 = 1
for (i in 1:(length(chain_5)-1)){
  if(chain_5[i+1] != chain_5[i]){n_diff_5= n_diff_5+1}
}
accept_rate_5 = n_diff_5/length(chain_5)
accept_rate_5


# table

xtable(data.frame(rbind(cbind(t(ess_s1),accept_rate_1,gel_1,gel_1_CI),
                        cbind(t(ess_s2),accept_rate_2,gel_2,gel_2_CI),
                        cbind(t(ess_s3),accept_rate_3,gel_3,gel_3_CI),
                        cbind(t(ess_s4),accept_rate_4,gel_4,gel_4_CI),
                        cbind(t(ess_s5),accept_rate_5,gel_5,gel_5_CI))),
       digits=1)


## part b


T = c(1,2,4,5,10,20,50)
M=length(T)
h = function(theta) -log(g(theta))
g_m = function(theta,m) exp(-h(theta)/T[m])
k = 1/integrate(g,lower=-10,upper=50)$value
g_norm = function(x) k*g(x)

abs = seq(-20,50,0.01)
ord = g_m(abs,m=1)

png(filename='unnormalized_gm.png')
plot(g,xlim=c(-20,50),ylim=c(0,1)
     ,main='Unnormalized densities of g_1, ... , g_M',
     ylab = 'g_m(x)')

for (i in 2:M){
  ord = g_m(abs,m=i)
  lines(abs,ord,col=i)
}

legend('topleft',paste('Tm = ',T),col=1:M,lwd=1,cex=0.7)
dev.off()

png(filename='normalized_gm.png')
plot(g_norm,xlim=c(-20,50),ylim=c(0,0.4)
     ,main='Normalized densities of g_1, ... , g_M',
     ylab = 'g_m(x)')

for (i in 2:M){
  k_i = 1/integrate(function(x) g_m(x,m=i),lower=-10,upper=50)$value
  ord = g_m(abs,m=i)*k_i
  lines(abs,ord,col=i)
}

legend('topleft',paste('Tm = ',T),col=1:M,lwd=1,cex=0.7)
dev.off()

q_b = function(l,m){
  if(l==1 & m==2){return(1)}
  if(l==M & m==M-1){return(1)}
  if(m==l-1 || m==l+1){return(0.5)}
  else{return(0)}
}

list_sigma = c(0.25,0.5,1,2,5,10,15)

xtable(data.frame(rbind(1:M,T,list_sigma)))

n=1e5

parallel_tempering = function(n_iter=n,init,liste_sigma=list_sigma){
  Xres = matrix(NA,ncol=n_iter,nrow=M)
  X = cbind(init)
  Xres[,1] = X
  for(t in 2:n_iter){
    X = Xres[,t-1]
    for(m in 1:M){
      Y_m = X[m]+rnorm(1,mean=0,sd=liste_sigma[m])
      if (runif(1) <= min((g_m(Y_m,m=m)/g_m(X[m],m=m)),1)) {X[m] = Y_m}
      Xres[m,t] = X[m]
    }
    for (p in 1:M){
    l = rdunif(1,a=1,b=M)
    if(l==1){m=2}
    if(l==M){m=M-1}
    if(l!= 1 & l!=M){
      u=runif(1)
      if(u<0.5){m=l+1}
      else{m=l-1}
    }
    v = runif(1)
    rapp = exp((h(Xres[l,t])-h(Xres[m,t])) * (1/T[l] - 1/T[m]))
    if (v < min(1,rapp*q_b(m,l)/q_b(l,m))){
      tilde_x = Xres[m,t]
      Xres[m,t] = Xres[l,t]
      Xres[l,t] = tilde_x
    }
    }
  }
  return(Xres)}

inits = rbind(rep(0,M),
             rep(10,M),
             rep(20,M),
             rep(30,M))

set.seed(42)

chains = sapply(1:dim(inits)[1], function(i) parallel_tempering(init=inits[i,])[1,])

chains_mcmc = lapply(1:dim(inits)[1], function(i) mcmc(chains[,i]))

chains_mcmc_list = mcmc.list(chains_mcmc)
gelman.diag(chains_mcmc_list)

init = rep(0,M)
chains = parallel_tempering(init=init,n_iter = 1e5)

chain_theta = chains[1,-seq(1,5000)]
png(filename='traceplot_theta.png')
traceplot(mcmc(chain_theta),main='Traceplot of the Markov chain simulated for theta',
          ylab='theta_t',xlab='iteration t')
abline(h=30,col='red',lty=2)
abline(h=a,col='red',lty=2)
dev.off()

effectiveSize(chain_theta)

png(filename='acf_theta.png')
acf(chain_theta,lag.max=300)
dev.off()

n_diff = 1
for (i in 1:(length(chain_theta)-1)){
  if(chain_theta[i+1] != chain_theta[i]){n_diff= n_diff+1}
}
accept_rate = n_diff/length(chain_theta)
accept_rate

png(filename = 'histogram_theta.png')
hist(chain_theta,breaks=1000,
     xlab='theta',
     main='Histogram of the Markov chain simulated for theta')
abline(v=a,lty=2,col='red')
abline(v=30,lty=2,col='red')
dev.off()

png(filename = 'density_plot.png')
plot(density(chain_theta,adjust=0.15),ylim=c(0,0.5),
     col='red',xlab='theta',ylab='g(theta)',
     main='True density vs sampled density')
ord = k*g(abs)
lines(abs,ord)
legend('topright',col=c(1,2),lwd=1,cex=0.8,legend=c('true density','sampled density'))
dev.off()

mean(chain_theta>20)

## Question 3


k = 1/integrate(g,lower=-10,upper=50)$value
k*integrate(g,lower=10,upper=60)$value




