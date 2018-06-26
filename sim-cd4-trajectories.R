#!/usr/bin/env Rscript

# description:  IPCW-MNAR data simulation
# programmer:   Dustin J. Rabideau

start <- proc.time()

args <- commandArgs(trailingOnly=T)

# user input
totalIter <- as.numeric(args[1]) # total simulations to be split across nCores
seed <- as.numeric(args[7]) ## 60000 for primary JSM 2018 poster figure

n <- as.numeric(args[2]) # number of obs ## 25000
K <- as.numeric(args[3]) # max follow up visits (base=0) ## 4
p.or <- as.numeric(args[4]) # P(outreach) ## 0.1

save.pdf <- T

# the following 2 log(OR) lead to MNAR (if either not equal log(1))
g3 <- log(as.numeric(args[5])) # log OR of cens for death ## 5
g4 <- log(as.numeric(args[6])) # log OR of cens for treatment ## 0.4
################

# Read in params from SLURM
nCores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
j <- Sys.getenv('SLURM_JOB_ID')
wd <- paste0(j, '/')
filename <- paste0(wd, 'out-', j, '.RData')

library(abind)
library(ggplot2)
library(doMC)

source('simMNAR.R') # function for running entire simulation

# run simulation
################################################################################

registerDoMC(nCores)
acomb <- function(...) abind(..., along=3)

# run
start1 <- proc.time()
results <- foreach(i=1:totalIter, .combine=acomb) %dopar% {
                     set.seed(seed+i)
                     simMNAR(n=n, K=K, p.or=p.or, gamma3=g3, gamma4=g4)
                   }
time1 <- proc.time() - start1
cat(totalIter, 'simulations took', time1['elapsed'], 'seconds\n')


# Plot cd4.post if save.pdf==T
##################################
cd4.post <- results[1:5, -1, ]

pdf.title2 <- paste0(wd, 'figure-', j, '.pdf')

# get "true" cd4 trajectories from large data generation
start2 <- proc.time()
ds.truth <- foreach(i=1:10, .combine=rbind) %dopar% {
  set.seed(seed+i)
  simMNAR(n=1e6, K=K, p.or=p.or, gamma3=g3, gamma4=g4, dataOnly=T)
}
truth <- as.vector(tapply(ds.truth$l, ds.truth$t, median))
rm(ds.truth)
time2 <- proc.time() - start2
cat('true cd4 trajectories took', time2['elapsed'], 'seconds\n')

truthMat <- matrix(rep(truth, each=totalIter),
                nrow=4, ncol=totalIter, byrow=T)  

obs.bias <- t(cd4.post['obs', , ] - truthMat)
mar.bias <- t(cd4.post['ipcw.mar', , ] - truthMat)
mnar.bias <- t(cd4.post['ipcw.mnar', , ] - truthMat)
mnar2.bias <- t(cd4.post['ipcw.mnar2', , ] - truthMat)
b <- as.data.frame(rbind(obs.bias, mar.bias, mnar.bias, mnar2.bias), row.names=F)
b$type <- rep(c('obs','mar','mnar','mnar2'), each=totalIter)
b$type.num <- rep(1:4, each=totalIter)

plot.data <- reshape(b, direction='long', varying=names(b)[1:K], v.names='bias')
means <- as.vector(tapply(plot.data$bias, plot.data[,c('type.num','time')], mean))

# tabulate simulation cd4.post
#############################

# bias
eThetaHat <- apply(cd4.post, c(1,2), mean)
bias <- eThetaHat[-1, ] - matrix(rep(truth, 4), nrow=4, byrow=T)
bias <- tapply(plot.data$bias, plot.data[,c('type.num','time')], mean)

# variance
cd4.post2 <- cd4.post
for (i in 1:dim(cd4.post2)[3]) {
  cd4.post2[, , i] <- (cd4.post2[, , i] - eThetaHat)^2
}
variance <- apply(cd4.post2, c(1,2), mean)[-1, ]

# mse
cd4.post3 <- cd4.post
for (k in 1:dim(cd4.post3)[3]) {
  for (i in 1:dim(cd4.post3)[1]) {
    cd4.post3[i, , k] <- (cd4.post3[i, , k] - truth)^2
  }
}
mse <- apply(cd4.post3, c(1,2), mean)[-1, ]
# mse - (variance + bias^2) # check

bias_plot <- sprintf("%.2f", round(as.vector(bias), 2))
variance_plot <- sprintf("%.2f", round(as.vector(variance), 2))
mse_plot <- sprintf("%.2f", round(as.vector(mse), 2))

if (save.pdf) pdf(pdf.title2, height=7, width=9)
  ylim <- range(plot.data$bias)
  boxplot(bias ~ type.num + time, data=plot.data, xlab='follow-up visit', lty=1,
          ylim=ylim, ylab='(estimate - truth)', 
          main=paste0('n = ', n, ', iterations = ', totalIter),
          col=c('firebrick','darkseagreen','steelblue','orange'), axes=F)
  axis(2, las=1)
  axis(1, at=(1:K) * 4 - 1, labels=1:max(plot.data$time))
  box()
  abline(h=0, lty=3)
  abline(v=(1:(K-1)) * 4 + 0.5) # lines between time points
  points(seq_along(means), means, pch=16, col='grey') # add mean (bias)
  legend('topleft', legend=c('unweighted', 'IPCW-MAR', 'IPCW-MNAR (pool)', 'IPCW-MNAR (sep)', 'mean'), bty='n', 
         pch=c(rep(15,4),16), pt.cex=c(rep(2,4),1),
         col=c('firebrick','darkseagreen','steelblue','orange','grey'))
  ys <- ylim[1] + (0:2) * (diff(ylim) / 30)
  text(seq_along(means)+0.3, ys[1], bias_plot, cex=0.7, adj=1)
  text(seq_along(means)+0.3, ys[2], variance_plot, cex=0.7, adj=1)
  text(seq_along(means)+0.3, ys[3], mse_plot, cex=0.7, adj=1)
  text(0.05, ys, c('bias', 'var', 'MSE'), cex=0.7, adj=0, font=2)
if (save.pdf) graphics.off()
  
# output
save(results, n, K, p.or, totalIter, g3, g4, seed, truth, file=filename)
  
time <- proc.time() - start
cat('All R code took', time['elapsed'], 'seconds')
