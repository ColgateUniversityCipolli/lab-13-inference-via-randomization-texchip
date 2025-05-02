
##################################
##### Task 1 #####################
##################################

##### a) #########################

library(tidyverse)
library(e1071)
set.seed(5656)
zf <- read_csv("zebrafinches.csv")

n <- length(zf$further)
skew <- skewness(zf$further)
skew.sqn <- skew/sqrt(n)
sd.zf <- sd(zf$further)
t <- mean(zf$further) / (sd.zf / sqrt(n))
(error <- (skew.sqn)*((2*t^2+1)/6)*dnorm(t))

##### b) #########################

t.vals <- seq(-10, 10, 0.1)
errors <- tibble(error = numeric(length(t.vals)))

for(i in 1:length(t.vals)){
  errors$error[i] <- (skew.sqn)*((2*t.vals[i]^2+1)/6)*dnorm(t.vals[i])
}

plot(t.vals, errors$error, type = "l", 
     main = "Error for t statistics", 
     ylab = "Error", xlab = "t")

##### c) #########################

alpha <- 0.05
t.crit <- qnorm(alpha)
(n.star <- ((skew/(6*0.10*alpha))*(2*t.crit^2+1)*dnorm(t.crit))^2)

##################################
##### Task 2 #####################
##################################

##### a) #########################

R <- 10000
resamples <- tibble(xbars = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(zf$further,
                          size = nrow(zf),
                          replace = T)
  
  resamples$xbars[i] <- mean(curr.resample)
}
r.xbar.further <- resamples$xbars

t.further <- r.xbar.further / (sd.zf / sqrt(n))
resamples.null.further <- t.further - mean(t.further)

for(i in 1:R){
  curr.resample <- sample(zf$closer,
                          size = nrow(zf),
                          replace = T)
  
  resamples$xbars[i] <- mean(curr.resample)
}
r.xbar.closer <- resamples$xbars

t.closer <- r.xbar.closer / (sd.zf / sqrt(n))
resamples.null.closer <- t.closer - mean(t.closer)

for(i in 1:R){
  curr.resample <- sample(zf$diff,
                          size = nrow(zf),
                          replace = T)
  
  resamples$xbars[i] <- mean(curr.resample)
}
r.xbar.diff <- resamples$xbars

t.diff <- r.xbar.diff / (sd.zf / sqrt(n))
resamples.null.diff <- t.diff - mean(t.diff)

##### b) #########################

t.test.further <- t.test(zf$further, mu = 0, alternative = "less")
(p.t.further <- t.test.further$p.value)
(p.boot.further <- mean(resamples.null.further <= t.test.further$statistic))

t.test.closer <- t.test(zf$closer, mu = 0, alternative = "greater")
(p.t.closer <- t.test.closer$p.value)
(p.boot.closer <- mean(resamples.null.closer >= t.test.closer$statistic))

t.test.diff <- t.test(zf$diff, mu = 0, alternative = "two.sided")
(p.t.diff <- t.test.diff$p.value)
(p.boot.diff <- mean(resamples.null.diff <= -t.test.diff$statistic) +
  mean(resamples.null.diff >= t.test.diff$statistic))

##### c) #########################

(quantile(resamples.null.further, alpha))
(quantile(resamples.null.closer, alpha))
(quantile(resamples.null.diff, alpha))
(qt(alpha,n-1))

##### d) #########################

t.test.further <- t.test(zf$further, mu = 0, alternative = "two.sided")
(t.test.further$conf.int)
(quantile(r.xbar.further, c(0.025, 0.975)))

t.test.closer <- t.test(zf$closer, mu = 0, alternative = "two.sided")
(t.test.closer$conf.int)
(quantile(r.xbar.closer, c(0.025, 0.975)))

(t.test.diff$conf.int)
(quantile(r.xbar.diff, c(0.025, 0.975)))

##################################
##### Task 3 #####################
##################################

##### a) #########################

R <- 10000
rand.further <- tibble(xbars = rep(NA, R))
rand.closer <- tibble(xbars = rep(NA, R))
rand.diff <- tibble(xbars = rep(NA, R))

further.shift <- zf$further - mean(zf$further)
closer.shift <- zf$closer - mean(zf$closer)
diff.shift <- zf$diff - mean(zf$diff)

for(i in 1:R){
  curr.rand.further <- further.shift *
    sample(x = c(-1, 1),
           size = length(further.shift),
           replace = T)
  curr.rand.closer <- closer.shift *
    sample(x = c(-1, 1),
           size = length(closer.shift),
           replace = T)
  curr.rand.diff <- diff.shift *
    sample(x = c(-1, 1),
           size = length(diff.shift),
           replace = T)

  rand.further$xbars[i] <- mean(curr.rand.further)
  rand.closer$xbars[i] <- mean(curr.rand.closer)
  rand.diff$xbars[i] <- mean(curr.rand.diff)
}

rand.further <- rand.further |>
  mutate(xbars = xbars + mean(zf$further)) # shifting back
rand.closer <- rand.closer |>
  mutate(xbars = xbars + mean(zf$closer)) # shifting back
rand.diff <- rand.diff |>
  mutate(xbars = xbars + mean(zf$diff)) # shifting back

##### b) #########################

(p.rand.further <- mean(rand.further$xbars - mean(zf$further) <= mean(zf$further)))

(p.rand.closer <- mean(rand.closer$xbars - mean(zf$closer) >= mean(zf$closer)))

(p.rand.diff <- mean(abs(rand.diff$xbars) - mean(zf$diff) >= abs(mean(zf$diff))))

##### c) #########################

R <- 1000
mu0.iterate <- 0.01
starting.point <- mean(zf$further)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zf$further - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(zf$further) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zf$further - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(zf$further) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}
# further
c(mu.lower, mu.upper)

R <- 1000
mu0.iterate <- 0.01
starting.point <- mean(zf$closer)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zf$closer - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(zf$closer) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zf$closer - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(zf$closer) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}
# closer
c(mu.lower, mu.upper)

R <- 1000
mu0.iterate <- 0.01
starting.point <- mean(zf$diff)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zf$diff - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(zf$diff) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zf$diff - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(zf$diff) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}
# diff
c(mu.lower, mu.upper)