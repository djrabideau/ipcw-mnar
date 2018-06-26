# This function generates CD4 trajectory data and compares 4 methods in a 
# simulation study: unweighted analysis, IPCW-MAR, IPCW-MNAR (pool),
# IPCW-MNAR (sep).
#   (sep) using separate models for uncensored/outreach 
#   (pool) using the same model with a fixed effect for outreach
#
# If gamma3 = gamma4 = 0, data is MAR. Otherwise it's MNAR.

simMNAR <- function(n=2500, K=4, p.or=0.1, gamma3, gamma4, dataOnly=F) {
  library(matrixStats)
    
  # n:          number of obs 
  # K:          max follow up visits (base=0)
  # p.or:       P(outreach)
  # dataOnly:   if T, just generates and returns data
  
  # create full set of obs
  id <- 1:n
  t <- 0:K
  ds.orig <- expand.grid(id, t, KEEP.OUT.ATTRS=F)
  names(ds.orig) <- c('id', 't')
  ds.orig <- ds.orig[with(ds.orig, order(id, t)), ]
  rownames(ds.orig) <- NULL
  
  # constants defined for loop(s)
  ##############################################################################
  
  # sim initial cd4
    shape <- 2
    scale <- 130 / shape
  # assume everyone at baseline uncens, on tx, alive
    ds.orig$c[with(ds.orig, t == 0)] <- 0
    ds.orig$d[with(ds.orig, t == 0)] <- 0
    ds.orig$a[with(ds.orig, t == 0)] <- 1
    ds.orig$p.c[with(ds.orig, t == 0)] <- 0
    ds.orig$p.d[with(ds.orig, t == 0)] <- 0
    ds.orig$p.a[with(ds.orig, t == 0)] <- 1
  # sim death: D[k] | L[k-1], A[k-1]
    alpha0 <- log(0.15) # log odds of D=1 for untreated, low cd4
    alpha1 <- log(0.3) # log OR of treatment
    alpha2 <- log(0.7) # log OR of high cd4
  # sim tx: A[k] | L[k-1], A[k-1]
    lambda0 <- log(9) # log odds of A=1 for current uncens, prev untx/lowCd4
    lambda1 <- log(4) # log OR of tx
    lambda2 <- log(2) # log OR of highCd4
  # sim cens: C[k] | L[k-1], A[k-1], D[k], A[k]
    gamma0 <- log(0.2) # log odds of C=1 for untx/prevlowCd4/alive
    gamma1 <- log(0.8) # log OR of treatment
    gamma2 <- log(0.5) # log OR of high cd4
  # sim cd4: L[k] | L[k-1], A[k-1]
    beta1 <- -10 # effect of 1-unit time increase for previous untreated
    beta2 <- 150 # effect of treatment
    sigma <- 50
    
  # begin simulation
  ##############################################################################
  
  ds <- ds.orig
  ds$l[with(ds, t == 0)] <- rgamma(n, shape=shape, scale=scale) # initial cd4
  ds$l.250 <- with(ds, as.numeric(l >= 250))
  
  for (k in 1:K) {
    # data for t[k-1]
    prev.ds <- ds[with(ds, t == (k - 1)), ]
    a <- prev.ds$a
    l <- prev.ds$l
    l.250 <- with(prev.ds, as.numeric(l >= 250))
    d <- prev.ds$d
    c <- prev.ds$c
    
    # sim death: D[k] | L[k-1], A[k-1]
    logit.d <- alpha0 + alpha1 * a + alpha2 * l.250
    p.d <- exp(logit.d) / (1 + exp(logit.d)) * (1 - d) + d # prob=1 if prev dead
    d.k <- rbinom(n, 1, p.d)
    
    # sim tx: A[k] | L[k-1], A[k-1]
    logit.a <- lambda0 + lambda1 * a + lambda2 * l.250
    p.a <- exp(logit.a) / (1 + exp(logit.a)) * (1 - d.k) # prob=0 if now dead
    a.k <- rbinom(n, 1, p.a)
    
    # sim cens: C[k] | L[k-1], A[k-1], D[k], A[k]
    logit.c <- gamma0 + gamma1 * a + gamma2 * l.250 + gamma3 * d.k + gamma4 * a.k
    p.c <- exp(logit.c) / (1 + exp(logit.c))
    c.k <- rbinom(n, 1, p.c) * (1 - c) * (1 - d) + c
  
    # sim cd4: L[k] | L[k-1], A[k-1]
    l.k <- pmax(0, (l + beta1 + beta2 * (1 / k) * a + rnorm(n, 0, sigma))) *
      (1 - d.k) - d.k  # create cd4 of -1 for dead patients
    
    # add to data
    ds$d[with(ds, t == k)] <- d.k
    ds$c[with(ds, t == k)] <- c.k
    ds$a[with(ds, t == k)] <- a.k
    ds$l[with(ds, t == k)] <- l.k
    ds$l.250[with(ds, t == k)] <- as.numeric(l.k >= 250)
    
    # add probs to data
    ds$p.d[with(ds, t == k)] <- p.d
    ds$p.c[with(ds, t == k)] <- p.c
    ds$p.a[with(ds, t == k)] <- p.a
  }
  
  ds$a[as.logical(ds$d)] <- NA # tx=NA after death
  
  # # check some numbers
  # tapply(ds$l, ds$t, median)
  # tapply(ds$l[ds$c == 0], ds$t[ds$c == 0], median)
  # # visualize data
  # p <- ggplot(data=ds[ds$id %in% sample(1:n, 100), ], aes(x=t, y=l, group=id))
  # p + geom_line()
  # tapply(ds$l, ds$t, quantile, prob=c(0.25, 0.5, 0.75))
  
  # lag variables for pooled logit models
  ds$c.lag1 <- c(NA, ds$c[-nrow(ds)])
  ds$c.lag1[which(!duplicated(ds$id))] <- NA
  ds$d.lag1 <- c(NA, ds$d[-nrow(ds)])
  ds$d.lag1[which(!duplicated(ds$id))] <- NA
  ds$a.lag1 <- c(NA, ds$a[-nrow(ds)])
  ds$a.lag1[which(!duplicated(ds$id))] <- NA
  ds$l.lag1 <- c(NA, ds$l[-nrow(ds)])
  ds$l.lag1[which(!duplicated(ds$id))] <- NA
  ds$l.250.lag1 <- c(NA, ds$l.250[-nrow(ds)])
  ds$l.250.lag1[which(!duplicated(ds$id))] <- NA
  
  # outreach
  ##############
  
  obs.or <- with(ds, c == 1 & c.lag1 + d.lag1 == 0)
  ds$o[obs.or] <- rbinom(sum(obs.or), 1, p.or)
  ds$o[!obs.or] <- 0
  
  # MAR analysis
  ###################
  
  ds.post <- ds[with(ds, t != 0), ]
  
  if (!dataOnly) {
    
    ds.c <- subset(ds.post, c.lag1 + d.lag1 == 0)
    fit.c <- glm(c ~ l.250.lag1 + a.lag1, data=ds.c, family=binomial())
    
    # put probs in dataframe
    ds.post$p.c0 <- 1 - predict(fit.c, ds.post, type='response')
    ds.post$p.c0[with(ds.post, d.lag1 == 1 & c.lag1 == 0)] <- 1 # prev obs dead
    ds.post$p.c0[with(ds.post, c.lag1 == 1)] <- 0 # prev censored
    ds.post$p.c1 <- 1 - ds.post$p.c0
    
    # cumulative probs by id over time
    ds.post$p.mar.cum <- unlist(tapply(ds.post$p.c0, ds.post$id, cumprod))
    
    # weights
    ds.post$w.mar <- 1 / ds.post$p.mar.cum
    ds.post$w.mar[with(ds.post, c == 1)] <- 0 # censored obs get 0 weight
    
    # # summary of MAR weights
    # tmp <- ds.post[with(ds.post, c == 0), ]
    # boxplot(tmp$w.mar ~ tmp$t, main='MAR weights', xlab='time', ylab='weight')
    # tapply(tmp$w.mar, tmp$t, sum)
    
    # MNAR analysis
    ###################
    
    # one model for uncensored and outreach
    #####
    ds.d <- subset(ds.c, c == 0 | o == 1)
    ds.a <- subset(ds.d, d == 0)
    fit.d <- glm(d ~ l.250.lag1 + a.lag1 + c, data=ds.d, family=binomial())
    fit.a <- glm(a ~ l.250.lag1 + a.lag1 + c, data=ds.a, family=binomial())
      
    # put probs in dataframe
    ds.predict.c0 <- ds.predict.c1 <- ds.post[, c('l.250.lag1','a.lag1','t')]
    ds.predict.c0$c <- 0
    ds.predict.c1$c <- 1
    ds.post$p.d1c0 <- predict(fit.d, ds.predict.c0, type='response')
    ds.post$p.d0c0 <- 1 - ds.post$p.d1c0
    ds.post$p.d1c1 <- predict(fit.d, ds.predict.c1, type='response')
    ds.post$p.d0c1 <- 1 - ds.post$p.d1c1
    ds.post$p.a1c0 <- predict(fit.a, ds.predict.c0, type='response')
    ds.post$p.a0c0 <- 1 - ds.post$p.a1c0
    ds.post$p.a1c1 <- predict(fit.a, ds.predict.c1, type='response')
    ds.post$p.a0c1 <- 1 - ds.post$p.a1c1
    
    # type
    ds.post$type <- 1 * with(ds.post, c.lag1 == 0 & d.lag1 == 1) + # died prior
               2 * with(ds.post, c.lag1 == 0 & d.lag1 == 0 & d == 1) + # died now
               3 * with(ds.post, c.lag1 == 0 & d == 0 & a == 1) + # alive on tx
               4 * with(ds.post, c.lag1 == 0 & d == 0 & a == 0) # alive off tx
    ds.post$type[ds.post$type == 0] <- NA
  
    # put MNAR probs in dataframe, modify censored obs probs
    type1 <- with(ds.post, type == 1 & !is.na(type))
    type2 <- with(ds.post, type == 2 & !is.na(type))
    type3 <- with(ds.post, type == 3 & !is.na(type))
    type4 <- with(ds.post, type == 4 & !is.na(type))
    
    ds.post$p.mnar[type1] <- 1
    ds.post$p.mnar[type2] <- with(ds.post[type2, ], (p.c0 * p.d1c0) / 
                            ((p.d1c1 * p.c1) + (p.d1c0 * p.c0)))
    ds.post$p.mnar[type3] <- with(ds.post[type3, ], (p.c0 * p.d0c0 * p.a1c0) / 
                            ((p.a1c1 * p.d0c1 * p.c1) + (p.a1c0 * p.d0c0 * p.c0)))
    ds.post$p.mnar[type4] <- with(ds.post[type4, ], (p.c0 * p.d0c0 * p.a0c0) / 
                            ((p.a0c1 * p.d0c1 * p.c1) + (p.a0c0 * p.d0c0 * p.c0)))
    ds.post$p.mnar[with(ds.post, c.lag1 == 1 & !is.na(c.lag1))] <- 0 # prev cens
    
    # cumulative probs by id over time
    ds.post$p.mnar.cum <- unlist(tapply(ds.post$p.mnar, ds.post$id, cumprod))
    
    # weights
    ds.post$w.mnar <- 1 / ds.post$p.mnar.cum
    ds.post$w.mnar[with(ds.post, c == 1)] <- 0 # censored obs get 0 weight
    
    # # summary of MNAR weights
    # tmp <- ds.post[with(ds.post, c == 0), ]
    # boxplot(tmp$w.mnar ~ tmp$t, main='MNAR weights', xlab='time', ylab='weight')
    # tapply(tmp$w.mnar, tmp$t, sum)
    
    
    
    # separate models for uncensored and outreach
    #####
    ds.d.unc <- subset(ds.c, c == 0)
    ds.a.unc <- subset(ds.d.unc, d == 0)
    fit.d.unc <- glm(d ~ l.250.lag1 + a.lag1, data=ds.d.unc, family=binomial())
    fit.a.unc <- glm(a ~ l.250.lag1 + a.lag1, data=ds.d.unc, family=binomial())
    
    ds.d.o <- subset(ds.c, o == 1)
    ds.a.o <- subset(ds.d.o, d == 0)
    fit.d.o <- glm(d ~ l.250.lag1 + a.lag1, data=ds.d.o, family=binomial())
    fit.a.o <- glm(a ~ l.250.lag1 + a.lag1, data=ds.d.o, family=binomial())
    
    # put probs in dataframe
    ds.post$p2.d1c0 <- predict(fit.d.unc, ds.post, type='response')
    ds.post$p2.d0c0 <- 1 - ds.post$p2.d1c0
    ds.post$p2.d1c1 <- predict(fit.d.o, ds.post, type='response')
    ds.post$p2.d0c1 <- 1 - ds.post$p2.d1c1
    ds.post$p2.a1c0 <- predict(fit.a.unc, ds.post, type='response')
    ds.post$p2.a0c0 <- 1 - ds.post$p2.a1c0
    ds.post$p2.a1c1 <- predict(fit.a.o, ds.post, type='response')
    ds.post$p2.a0c1 <- 1 - ds.post$p2.a1c1
    
    ds.post$p2.mnar[type1] <- 1
    ds.post$p2.mnar[type2] <- with(ds.post[type2, ], (p.c0 * p2.d1c0) / 
                            ((p2.d1c1 * p.c1) + (p2.d1c0 * p.c0)))
    ds.post$p2.mnar[type3] <- with(ds.post[type3, ], (p.c0 * p2.d0c0 * p2.a1c0) / 
                            ((p2.a1c1 * p2.d0c1 * p.c1) + (p2.a1c0 * p2.d0c0 * p.c0)))
    ds.post$p2.mnar[type4] <- with(ds.post[type4, ], (p.c0 * p2.d0c0 * p2.a0c0) / 
                            ((p2.a0c1 * p2.d0c1 * p.c1) + (p2.a0c0 * p2.d0c0 * p.c0)))
    ds.post$p2.mnar[with(ds.post, c.lag1 == 1 & !is.na(c.lag1))] <- 0 # prev cens
    
    # cumulative probs by id over time
    ds.post$p2.mnar.cum <- unlist(tapply(ds.post$p2.mnar, ds.post$id, cumprod))
    
    # weights
    ds.post$w2.mnar <- 1 / ds.post$p2.mnar.cum
    ds.post$w2.mnar[with(ds.post, c == 1)] <- 0 # censored obs get 0 weight
  
  
    
    # results of each method
    truth <- as.vector(tapply(ds.post$l, ds.post$t, median))
    
    obs <- as.vector(tapply(ds.post$l[ds.post$c == 0], 
                            ds.post$t[ds.post$c == 0], median))
    ipcw.mar <- unlist(lapply(split(ds.post, ds.post$t), 
                              function(x) weightedMedian(x$l, x$w.mar)))
    ipcw.mnar <- unlist(lapply(split(ds.post, ds.post$t), 
                               function(x) weightedMedian(x$l, x$w.mnar)))
    ipcw.mnar2 <- unlist(lapply(split(ds.post, ds.post$t), 
                               function(x) weightedMedian(x$l, x$w2.mnar)))
    res <- rbind(truth, obs, ipcw.mar, ipcw.mnar, ipcw.mnar2)
    res <- cbind(tapply(ds$l, ds$t, median)[1], res)
    dimnames(res)[[2]][1] <- '0'
    
    # other summary stats to store
    n.c <- tapply(ds$c, ds$t, sum)
    n.d <- tapply(ds$d, ds$t, sum)
    n.a <- tapply(ds$a, ds$t, sum, na.rm=T)
    n.l.250 <- tapply(ds$l.250, ds$t, sum)
    n.o <- tapply(ds$o, ds$t, sum)
      # prop among those just recently censored
    prop.o <- c(0, tapply(ds.post$o[ds.post$c == 1 & ds.post$c.lag1 == 0], 
                     ds.post$t[ds.post$c == 1 & ds.post$c.lag1 == 0], mean))
      # weight summary stats among those included in analysis
    mean.w.mar <- c(1, tapply(ds.post$w.mar[ds.post$w.mar != 0],
                              ds.post$t[ds.post$w.mar != 0], mean))
    min.w.mar <- c(1, tapply(ds.post$w.mar[ds.post$w.mar != 0],
                              ds.post$t[ds.post$w.mar != 0], min))
    max.w.mar <- c(1, tapply(ds.post$w.mar[ds.post$w.mar != 0],
                              ds.post$t[ds.post$w.mar != 0], max))
    med.w.mar <- c(1, tapply(ds.post$w.mar[ds.post$w.mar != 0],
                              ds.post$t[ds.post$w.mar != 0], median))
    mean.w.mnar <- c(1, tapply(ds.post$w.mnar[ds.post$w.mnar != 0],
                              ds.post$t[ds.post$w.mnar != 0], mean))
    min.w.mnar <- c(1, tapply(ds.post$w.mnar[ds.post$w.mnar != 0],
                              ds.post$t[ds.post$w.mnar != 0], min))
    max.w.mnar <- c(1, tapply(ds.post$w.mnar[ds.post$w.mnar != 0],
                              ds.post$t[ds.post$w.mnar != 0], max))
    med.w.mnar <- c(1, tapply(ds.post$w.mnar[ds.post$w.mnar != 0],
                              ds.post$t[ds.post$w.mnar != 0], median))
    mean.w2.mnar <- c(1, tapply(ds.post$w2.mnar[ds.post$w2.mnar != 0],
                              ds.post$t[ds.post$w2.mnar != 0], mean))
    min.w2.mnar <- c(1, tapply(ds.post$w2.mnar[ds.post$w2.mnar != 0],
                              ds.post$t[ds.post$w2.mnar != 0], min))
    max.w2.mnar <- c(1, tapply(ds.post$w2.mnar[ds.post$w2.mnar != 0],
                              ds.post$t[ds.post$w2.mnar != 0], max))
    med.w2.mnar <- c(1, tapply(ds.post$w2.mnar[ds.post$w2.mnar != 0],
                              ds.post$t[ds.post$w2.mnar != 0], median))
    
    stats <- rbind(n.c, n.d, n.a, n.l.250, n.o, prop.o, 
                   mean.w.mar, min.w.mar, max.w.mar, med.w.mar, 
                   mean.w.mnar, min.w.mnar, max.w.mnar, med.w.mnar,
                   mean.w2.mnar, min.w2.mnar, max.w2.mnar, med.w2.mnar)
  }

  if (dataOnly) {
    results <- ds.post
  } else {
    results <- rbind(res, stats)
  }
  
  return(results)
}
