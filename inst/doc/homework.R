## -----------------------------------------------------------------------------
library(latex2exp)
library(scales)
library(pander)

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
#x <- b/((1-u)^(1/a)) F(x) = 1-(b/x)^a, x>=b>0, a>0
x <- 2/((1-u)^(1/2))
hist(x, prob = TRUE, xlim = c(0,40), breaks = 100, main = expression(f(x)==frac(8,x^3)))  
y <- seq(2, 50, .01)
lines(y, 8/y^3) 

## -----------------------------------------------------------------------------
# acceptance-rejection method for Beta distribution
n <- 1e3; j<-k<-0; y <- numeric(n)
armethod_beta <- function(a,b){
  while (k < n) {
    u <- runif(1)
    j <- j + 1
    x <- runif(1) #random variate from g(.)
    if (x^(a-1) * (1-x)^(b-1) > u) {
      #we accept x
      k <- k + 1
      y[k] <- x
    }
  }
  return (list(j,y))
}

## -----------------------------------------------------------------------------
# a=3,b=2 
n <- 1e3; j<-k<-0; y <- numeric(n)
while (k < n) {
  u <- runif(1)
  j <- j + 1
  x <- runif(1) #random variate from g(.)
  if (x^2 * (1-x) > 16/9*u) {
    #we accept x
    k <- k + 1
    y[k] <- x
  }
}
hist(y, prob = TRUE, main = expression(f(y)==12*y^2*(1-y)))
z <- seq(0, 1, .01)
lines(z,dbeta(z,shape1=3,shape2=2))

## -----------------------------------------------------------------------------
n <- 1e3; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
y <- rexp(n, lambda) # the length of lambda = n
hist(y, prob = TRUE, xlim = c(0,15), breaks = 100, main = "Exponential-Gamma mixture")  

## -----------------------------------------------------------------------------
n <- 1e3; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
y <- rexp(n, lambda) # the length of lambda = n
hist(y, prob = TRUE, xlim = c(0,15), breaks = 100, main = expression(f(y)==frac(64,(2+y)^5)))  
z <- seq(0, 50, .01)
lines(z, 64*(2+z)^(-5)) 

## -----------------------------------------------------------------------------
# the fast sorting algorithm
fast_sorting <- function(x){
  num = length(x)
  if (num==0||num==1){
    return(x)
  }
  else{
    a = x[1]
    y = x[-1]
    lower = y[y<a]
    upper = y[y>=a]
    return(c(fast_sorting(lower),a,fast_sorting(upper)))
  }
}

a = numeric(5)
t = c(10^4, 2*10^4, 4*10^4, 6*10^4, 8*10^4)
tn = t*log(t)
for (i in 1:5){
  sum_time = 0
  for (j in 1:100){
    test = sample(t[i])
    sum_time = sum_time + system.time(fast_sorting(test))[1]
  }
  a[i] = sum_time/100
}
a

fit = lm(a~tn)
fit$coef

plot(x=tn, y=a, main = "computation time - nlog(n)")
y <- seq(0, 10^6, 1)
lines(y, fit$coef[1]+fit$coef[2]*y, col="red")

## -----------------------------------------------------------------------------
cov = exp(1)-(exp(1)-1)^2 # compute Cov
round(cov,5)
var = 2*((exp(1)^2-1)/2-(exp(1)-1)^2) + 2*cov # compute Var
round(var,5)

round(1-var/((exp(1)^2-1)/2-(exp(1)-1)^2),5) # the percent reduction

# the calculated result
MC <- function(x, R = 10000, antithetic = FALSE) {
  u <- runif(R/2)
  if (antithetic) v <- 1 - u else v <- runif(R/2)
  u <- c(u, v)
  g <- exp(u*x) 
  cdf <- mean(g)
  cdf
}
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1
for (i in 1:m) {
  MC1[i] <- MC(x, R = 1000, antithetic = FALSE)
  MC2[i] <- MC(x, R = 1000, antithetic = TRUE)
}
c(mean(MC1),mean(MC2)) # theta's estimate：MC method and antithetic variate approach
round(c(sd(MC1)^2,sd(MC2)^2,1-sd(MC2)^2/sd(MC1)^2),7) # theta's variance：MC method and antithetic variate approach

# the theoretical value
var_antithetic = var/(2*m)
var_MC = ((exp(1)^2-1)/2-(exp(1)-1)^2)/m
round(c(var_antithetic,var_MC,1-var_antithetic/var_MC),7)

## -----------------------------------------------------------------------------
g <- function(x){
  return (x^2/sqrt(2*pi)*exp(-x^2/2))
}
f1 <- function(x){
  return (exp(1-x))
}
f2 <- function(x){
  return (1/(x^2))
}

# The true value of the integral
I <- integrate(g,1,Inf)

m <- 1e6
set.seed(516)
u <- runif(m)

# the method that uses f1 as an importance function
x1 <- 1-log(1-u)
fg1 <- g(x1) / f1(x1)
I1_hat <- mean(fg1)
sd1 <- sd(fg1)

# the method that uses f2 as an importance function
x2 <- 1/(1-u)
fg2 <- g(x2) / f2(x2)
I2_hat <- mean(fg2)
sd2 <- sd(fg2)

round(c(I[[1]], I1_hat, I2_hat),5) # The true value of the integral and two estimations
round(c(sd1,sd2),5) # stand error of f1 and f2

## -----------------------------------------------------------------------------
x <- seq(1,6,0.01)
plot(x, g(x), col="black", type="l", ylim=c(0,1), main='(A)')
lines(x, f1(x), col="orange")
lines(x, f2(x), col="green")
legend("topright", col=c("black","orange","green"), lty=1, legend=c("g","f1","f2"))

## -----------------------------------------------------------------------------
M <- 10000; k <- 5 # what if k is larger?
set.seed(516)
r <- M/k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
est <- matrix(0, N, 2)
g <- function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
p <- c(0,1/5,2/5,3/5,4/5,1)
q <- numeric(5)
for (i in 1:5){
  q[i] <- (1-exp(-1)) / (exp(-p[i]) - exp(-p[i+1]))
}
for (i in 1:N) {
  est[i,1] <- mean(g(runif(M)))
  for(j in 1:k){
    u <- runif(M) #f3, inverse transform method
    x <- - log(1 - u * (1 - exp(-1)))
    fg <- g(x) / q[j] / (exp(-x) / (1 - exp(-1)))
    T2[j] <- mean(fg)
  }
  est[i,2] <- T2[1]+T2[2]+T2[3]+T2[4]+T2[5]
}
round(c(0.5257801,apply(est,2,mean)),4)
round(c(0.0970314,apply(est,2,sd)),5)

## -----------------------------------------------------------------------------
lognormal_MC <- function(n, ml, sl){
  DCL <- numeric(1000);UCL <- numeric(1000)
  for (i in 1:1000){
    y <- rnorm(n, ml, sl)
    DCL[i] <- mean(y) - sd(y)/sqrt(n) * qt(0.975, df = n-1)
    UCL[i] <- mean(y) + sd(y)/sqrt(n) * qt(0.975, df = n-1)
  }
  D <- exp(mean(DCL));U <- exp(mean(UCL))
  rm(y);rm(DCL);rm(UCL) # clear memory
  return (c(D, U))
}

## -----------------------------------------------------------------------------
lognormal_MC(20, 0, 2) # example

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  rm(X) # clear memory
  rm(Y) # clear memory
  return(as.integer(max(c(outx, outy)) > 5))
}

Ftest <- function(x,y){
  return (as.integer(var.test(x,y,alternative = "two.sided",conf.level = 0.945)[[3]]<0.055))
}

## -----------------------------------------------------------------------------
sigma1 <- 1
sigma2 <- 1.5
compare_power <- function(m,n){
  sum_c = sum_f =0
  for (i in 1:m){
    x <- rnorm(n, 0, sigma1)
    y <- rnorm(n, 0, sigma2)
    sum_c = sum_c + count5test(x, y)
    sum_f = sum_f + Ftest(x, y)
  }
  power_c <- sum_c/m
  power_f <- sum_f/m
  rm(x) # clear memory
  rm(y) # clear memory
  return (c(power_c,power_f))
}

## -----------------------------------------------------------------------------
print("the power of the Count Five test, the power of F Test")
compare_power(1000,20) # Compare the power of the Count Five test and F test for small sample size(n=20).
compare_power(1000,100) # Compare the power of the Count Five test and F test for medium sample size(n=100).
compare_power(1000,1000) # Compare the power of the Count Five test and F test for large sample size(n=1000).

## -----------------------------------------------------------------------------
# bootstrap
x <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
B <- 10000
set.seed(516)
lambdastar <- numeric(B)
lambda <- 1/mean(x)
for(b in 1:B){
  xstar <- sample(x, replace=TRUE)
  lambdastar[b] = 1 / mean(xstar)
}

## -----------------------------------------------------------------------------
c(lambda.bootest = mean(lambdastar), lambda.MLE = lambda, bias = mean(lambdastar) - lambda, se.boot = sd(lambdastar))
rm(lambdastar) # clear memory

## -----------------------------------------------------------------------------
library(boot)
set.seed(516)
R <- 2000
skewness <- function(x,i){
  x_bar <- mean(x[i])
  x_bar
}
bo <- boot(data=x, statistic=skewness, R = R)
ci <- boot.ci(bo, type=c("norm", "basic", "perc", "bca"))

## -----------------------------------------------------------------------------
print(bo)
print(ci)
rm(bo) # clear memory
rm(ci) # clear memory

## -----------------------------------------------------------------------------
m <- 10000
ci.norm <- ci.basic <- ci.percent <- matrix(NA, m, 2)
skewness <- function(x,i){
  x_bar <- mean(x[i])
  x_bar
}
n <- 10
mu <- 0;sigma <- 1 # standard normal distribution
for(i in 1:m){
  s <- rnorm(n, mu, sigma) 
  de <- boot(s, statistic=skewness, R=1000)
  ci <- boot.ci(de, type=c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.percent[i,] <- ci$percent[4:5]
}

## -----------------------------------------------------------------------------
t <- matrix(nrow=3, ncol=3)
rownames(t) <- c("normal","basic","perc")
colnames(t) <- c("coverage.prob", "missed.left", "missed.right")

t[,2] <- c(sum(mu<ci.norm[,1]),sum(mu<ci.basic[,1]),sum(mu<ci.percent[,1]))/m
t[,3] <- c(sum(mu>ci.norm[,2]),sum(mu>ci.basic[,2]),sum(mu>ci.percent[,2]))/m
t[,1] <- 1-t[,2]-t[,3]
t
rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
n <- nrow(scor)

sigma_hat <- cov(scor)
theta.hat <- eigen(sigma_hat)$values[1] / sum(eigen(sigma_hat)$values)

# compute the jackknife replicates, leave-one-out estimates
theta.jack <- numeric(n)
for (i in 1:n){
  s <- cov(scor[-i,])
  theta.jack[i] <- eigen(s)$values[1] / sum(eigen(s)$values)
}

## -----------------------------------------------------------------------------
sigma_hat
theta.hat

bias <- (n-1) * (mean(theta.jack) - theta.hat)
print(bias) #jackknife estimate of bias

se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2))
print(se)

#The jackknife estimate of standard error is 0.04955231. From the previous result for the bias, we have the estimated coefficient of variation
bias/se

rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
library(DAAG,quietly=TRUE); attach(ironslag)
n <- length(magnetic)   #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)

# fit models on leave-two-out samples
i = 1
for (k in 1:(n-1)) {
  for (l in (k+1):n) {
    y <- magnetic[-c(k,l)]
    x <- chemical[-c(k,l)]

    J1 <- lm(y ~ x)
    yhat1k <- J1$coef[1] + J1$coef[2] * chemical[k]
    yhat1l <- J1$coef[1] + J1$coef[2] * chemical[l]
    e1[i] <- ((magnetic[k] - yhat1k)^2 + (magnetic[l] - yhat1l)^2)/2

    J2 <- lm(y ~ x + I(x^2))
    yhat2k <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
    yhat2l <- J2$coef[1] + J2$coef[2] * chemical[l] + J2$coef[3] * chemical[l]^2
    e2[i] <- ((magnetic[k] - yhat2k)^2 + (magnetic[l] - yhat2l)^2)/2

    J3 <- lm(log(y) ~ x)
    logyhat3k <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3k <- exp(logyhat3k)
    logyhat3l <- J3$coef[1] + J3$coef[2] * chemical[l]
    yhat3l <- exp(logyhat3l)
    e3[i] <- ((magnetic[k] - yhat3k)^2 + (magnetic[l] - yhat3l)^2)/2

    J4 <- lm(log(y) ~ log(x))
    logyhat4k <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    yhat4k <- exp(logyhat4k)
    logyhat4l <- J4$coef[1] + J4$coef[2] * log(chemical[l])
    yhat4l <- exp(logyhat4l)
    e4[i] <- ((magnetic[k] - yhat4k)^2 + (magnetic[l] - yhat4l)^2)/2
        
    i = i+1
  }
}

## -----------------------------------------------------------------------------
c(mean(e1), mean(e2), mean(e3), mean(e4))
    
rm(e1) # clear memory
rm(e2) # clear memory
rm(e3) # clear memory
rm(e4) # clear memory

## -----------------------------------------------------------------------------
set.seed(210)
x <- rnorm(20,0,1)
y <- rnorm(20,0,1)
z <- c(x, y) 
R <- 999
reps <- numeric(R) 
cor0 <- cor(x, y, method = "spearman")

for (i in 1:999) {
  k <- sample(1:40, size = 20, replace = FALSE)
  x1 <- z[k];y1 <- z[-k]
  reps[i] <- cor(x1, y1, method = "spearman")
}
p <- mean(c(cor0, reps) >= cor0)

## -----------------------------------------------------------------------------
round(c(p,cor.test(x,y)$p.value),5)
rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
set.seed(516)
rw.Metropolis <- function(n, sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
      if (u[i] <= exp(abs(x[i-1])-abs(y)))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
n <- 4 # degrees of freedom for target Student t dist.
N <- 2000
sigma <- c(.05, .5, 2, 8)

x0 <- 25
rw1 <- rw.Metropolis(n, sigma[1], x0, N)
rw2 <- rw.Metropolis(n, sigma[2], x0, N)
rw3 <- rw.Metropolis(n, sigma[3], x0, N)
rw4 <- rw.Metropolis(n, sigma[4], x0, N)
# number of candidate points rejected
# print(c(rw1$k, rw2$k, rw3$k, rw4$k))
print(round(c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N),3)) # the acceptance rates of each chain

#par(mfrow=c(2,2))
#plot(rw1$x, type="l", xlab="sigma=0.05", ylab="X", ylim=range(rw1$x))
#plot(rw2$x, type="l", xlab="sigma=0.5", ylab="X", ylim=range(rw2$x))
#plot(rw3$x, type="l", xlab="sigma=2", ylab="X", ylim=range(rw3$x))
#plot(rw4$x, type="l", xlab="sigma=8", ylab="X", ylim=range(rw4$x))

## -----------------------------------------------------------------------------
set.seed(516)

Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

normal.chain <- function(sigma, N, X1) {
  # generates a Metropolis chain for Normal(0,1)
  # with Normal(X[t], sigma) proposal distribution
  # and starting value X1
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)

  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
      if (u[i] <= exp(abs(x[i-1])-abs(y)))
      x[i] <- y else {
        x[i] <- x[i-1]
      }
  }
  return(x)
}

# sigma <- 2 #parameter of proposal distribution
k <- 4 # number of chains to generate
n <- 15000 # length of chains
b <- 1000 # burn-in length

# choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)

sigma.chain <- function(sigma) {
# generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- normal.chain(sigma, n, x0[i])

# trace plots
plot(1:n,X[1,],type="l")
lines(1:n,X[2,],type="l",col=2)
lines(1:n,X[3,],type="l",col=3)
lines(1:n,X[4,],type="l",col=4)

# compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

# plot psi for the four chains
par(mfrow=c(2,2))
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) # restore default

# plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)
}

## -----------------------------------------------------------------------------
#sigma.chain(.05)
#sigma.chain(.5)
#sigma.chain(2)
#sigma.chain(8)

rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
set.seed(516)
# initialize constants and parameters
N <- 10000 # length of chain
burn <- 1000 # burn-in length
X <- matrix(0, N, 2) # the chain, a bivariate sample

rho <- 0.9 # correlation
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

###### generate the chain #####

X[1, ] <- c(mu1, mu2) # initialize
for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}

Y <- t(X)
b <- burn + 1
x <- X[b:N, ]

Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
# compare sample statistics to parameters
cor(x)
cat('Means: ',round(colMeans(x),3),'\n')
cat('Standard errors: ',round(apply(x,2,sd),3),'\n')
cat('Correlation coefficients: ', round(cor(x[,1],x[,2]),3),'\n')

## -----------------------------------------------------------------------------
plot(x, main="", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(x[,2]))

plot(x[,1], type='l', col=1, lwd=2, xlab='Index', ylab='Random numbers')
lines(x[,2], col=2, lwd=2)
legend('bottomright', c(expression(X[1]),expression(X[2])), col=1:2, lwd=2)

L <- lm(x[,2] ~ x[,1])
L
summary(L)

## -----------------------------------------------------------------------------
#n <- 10000
#b <- 1000

# compute diagnostic statistics
#psi <- t(apply(Y, 1, cumsum))
#for (i in 1:nrow(psi))
#  psi[i,] <- psi[i,] / (1:ncol(psi))
#print(Gelman.Rubin(psi))

# plot psi for the four chains
#par(mfrow=c(2,2))
#for (i in 1:2)
#  plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
#par(mfrow=c(1,1)) # restore default

# plot the sequence of R-hat statistics
#rhat <- rep(0, n)
#for (j in (b+1):n)
#  rhat[j] <- Gelman.Rubin(psi[,1:j])
#plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#abline(h=1.2, lty=2)

rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
#Model <- function(alpha, beta, gamma, N){
#  x <- rnorm(N, 0, 1)
#  e_M <- rnorm(N, 0, 1)
#  e_Y <- rnorm(N, 0, 1)
#  m <- alpha*x + e_M
#  y <- beta*m + gamma*x + e_Y
#  return(data.frame(x, m, y))
#}
#
#t <- function(x, m, y){
#  data <- data.frame(x, m, y)
#  Model_M <- lm(m~x, data=data)
#  Model_Y <- lm(y~m+x, data=data)
#  t <- Model_M$coef[2]*Model_Y$coef[2]/sqrt(Model_M$coef[2]^2*(summary(Model_M)$coefficients[2,2])^2+Model_Y$coef[2]^2*(summary(Model_Y)$coefficients[2,2])^2)
#  return(as.numeric(t))
#}

#Permutation_Test_1 <- function(data){
#  N <- 10; B <- 499
#  reps <- numeric(B)
#  for (i in 1:B) {
#    k <- sample(1:N, size = N, replace = FALSE)
#    x <- data[k, 1]
#    reps[i] <- t(x, data$m, data$y)
#  }
#  p_hat <- (sum(abs(reps)>abs(t(data$x, data$m, data$y)))+1)/(1+B)
#  return(p_hat<0.05)
#}
#
#Permutation_Test_2 <- function(data){
#  N <- 10; B <- 499
#  reps <- numeric(B)
#  for (i in 1:B) {
#    k <- sample(1:N, size = N, replace = FALSE)
#    y <- data[k, 3]
#    reps[i] <- t(data$x, data$m, y)
#  }
#  p_hat <- (sum(abs(reps)>abs(t(data$x, data$m, data$y)))+1)/(1+B)
#  return(p_hat<0.05)
#}
#
#Permutation_Test_3 <- function(data){
#  N <- 10; B <- 499
#  reps <- numeric(B)
#  for (i in 1:B) {
#    k <- sample(1:N, size = N, replace = FALSE)
#    m <- data[k, 2]
#    reps[i] <- t(data$x, m, data$y)
#  }
#  p_hat <- (sum(abs(reps)>abs(t(data$x, data$m, data$y)))+1)/(1+B)
#  return(p_hat<0.05)
#}
#
#simulation <- function(alpha, beta, gamma=1, N=10){
#  Type_1_error <- matrix(0, nrow = 3, ncol = 200)
#  for (i in 1:200) {
#    data <- Model(alpha, beta, gamma, N)
#    Type_1_error[1,i] <- Permutation_Test_1(data)
#    Type_1_error[2,i] <- Permutation_Test_2(data)
#    Type_1_error[3,i] <- Permutation_Test_3(data)
#  }
#  return(c(mean(Type_1_error[1,]),mean(Type_1_error[2,]),mean(Type_1_error[3,])))
#}

## -----------------------------------------------------------------------------
#set.seed(516)
#simulation1 <- simulation(alpha=0, beta=0)
#simulation2 <- simulation(alpha=0, beta=1)
#simulation3 <- simulation(alpha=1, beta=0)
#Permutation_Test <- data.frame(Model1=simulation1, Model2=simulation2, Model3=simulation3)
#row.names(Permutation_Test) <- c("Permutation_Test_1", "Permutation_Test_2", "Permutation_Test_3")
#Permutation_Test
#
#rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
gfo <- function(N, b1, b2, b3, f0, x1, x2, x3){
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
}

alpha_estimation <- function(N, b1, b2, b3, f0){
  x1 <- rpois(N,1); x2 <- rexp(N); x3 <- sample(0:1, N, replace=TRUE)
  
  solution <- uniroot(gfo(N, b1, b2, b3, f0, x1, x2, x3),c(-15,5))
  round(unlist(solution),5)[1:3]

  alpha <- solution$root

  tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
  p <- 1/(1+tmp)
  d <- rbinom(N,1,p)
  mean(d) # estimation of alpha
}

## -----------------------------------------------------------------------------
set.seed(516)
N <- 1e6; b1 <- 0; b2 <- 1; b3 <- -1
alpha_estimation(N, b1, b2, b3, 0.1)
alpha_estimation(N, b1, b2, b3, 0.01)
alpha_estimation(N, b1, b2, b3, 0.001)
alpha_estimation(N, b1, b2, b3, 0.0001)

x = -log(c(0.1, 0.01, 0.001, 0.0001), 10)
y = -log(c(alpha_estimation(N, b1, b2, b3, 0.1), alpha_estimation(N, b1, b2, b3, 0.01), 
           alpha_estimation(N, b1, b2, b3, 0.001), alpha_estimation(N, b1, b2, b3, 0.0001)), 10)	
plot(x, y, type="p", main = "Scatter Diagram", xlab = "-log(f0)", ylab = "-log(alpha)", xlim = c(0,5), ylim = c(0,5))

rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
# maximizing likelihood function
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)

max_likelihood <- function(lam){
  return(sum((v*exp(-lam*v)-u*exp(-lam*u))/(exp(-lam*u)-exp(-lam*v))))
}

round(uniroot(max_likelihood, c(0, 10))$root, 5)

## -----------------------------------------------------------------------------
# EM algorithm
lam <- 2
i <- 0
eps <- 1e-3
lam1 <- 1
while(1){
  lam1 <- lam1 - max_likelihood(lam)/length(u) 
  if (abs(1/lam1 - lam) < eps) break
  i <- i+1
  lam <- 1/lam1
}
round(c(iter = i, lambda = lam), 5)

rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
unlist(list(list(1, 2), c(3, 4)))
as.vector(list(list(1, 2), c(3, 4)))

## -----------------------------------------------------------------------------
1 == "1"

## -----------------------------------------------------------------------------
-1 < FALSE
0 == FALSE

## -----------------------------------------------------------------------------
"one" < 2
"one" > 2

## -----------------------------------------------------------------------------
dim(c(1,2,3,4,5,6))

## -----------------------------------------------------------------------------
x <- matrix(c(1:12), nrow = 4, byrow = TRUE)

is.matrix(x)
is.array(x)

## -----------------------------------------------------------------------------
d <- data.frame(name=c("Zhang San", "Li Si", "Wang Wu"), age=c(30, 35, 28),
     height=c(182, 164, 175))
d
as.matrix(d)

## -----------------------------------------------------------------------------
# a data frame with 0 rows
# method 1
test1=data.frame(col1=character(0), col2=numeric(0), col3=logical(0))
test1
str(test1)

# method 2
test2=matrix(data=NA, nrow=0, ncol=3)
test2=as.data.frame(test2)
colnames(test2)=c("col1", "col2", "col3")
test2
str(test2)

# a data frame with 0 cols
# method 1
test3=matrix(data=NA, nrow=3, ncol=0)
test3=as.data.frame(test3)
rownames(test3)=c("row1", "row2", "row3")
test3
str(test3)

# a data frame with 0 rows and 0 cols
# method 1
data.frame()

rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
#library(tidyverse)
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
d <- data.frame(col1=c(2, 1, 0), col2=c(5, 1, 6), col3=c(4, 2, 3))
scale01_apply <- apply(d, 2, scale01)
scale01_apply

d2 <- data.frame(col1=c(2, 1, 0), col2=c(5, 1, 6), col3=c("a","b","c"))
#lapply(select_if(d2,is.numeric), scale01)

## -----------------------------------------------------------------------------
rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
d_a <- data.frame(col1=c(1, 2, 3), col2=c(1, 3, 5), col3=c(1, 4, 7))
d_a
sd_vapply <- vapply(d_a, sd, numeric(1))
sd_vapply

## -----------------------------------------------------------------------------
d_b <- data.frame(name=c("Zhang San", "Li Si", "Wang Wu"), age=c(30, 35, 28),
     height=c(182, 164, 175))
d_b
vapply(d_b[,which(vapply(
  d_b,is.numeric,FUN.VALUE = logical(1)))]
  , sd, FUN.VALUE = numeric(1))

## -----------------------------------------------------------------------------
rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
set.seed(516)
gibbsR <- function(mu1, mu2, sigma1, sigma2, rho, N, burn){
  X <- matrix(0, N, 2)
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  X[1, ] <- c(mu1, mu2) # initialize
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] <- rnorm(1, m1, s1)
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <- rnorm(1, m2, s2)
  }
  Y <- t(X)
  b <- burn + 1
  x <- X[b:N, ]
  return(x)
}

## -----------------------------------------------------------------------------
library(Rcpp)
# Can create source file in Rstudio
sourceCpp('../src/gibbsC.cpp')
library(microbenchmark)
ts <- microbenchmark(gibbsR=gibbsR(0, 0, 1, 1, 0.9, 10000, 1000),
                     gibbscpp=gibbsC(0, 0, 1, 1, 0.9, 10000, 1000))
summary(ts)[,c(1,3,5,6)]

## -----------------------------------------------------------------------------
qqplot(gibbsR(0, 0, 1, 1, 0.9, 10000, 1000), gibbsC(0, 0, 1, 1, 0.9, 10000, 1000), main="qqplot", xlab='R', ylab='Rcpp')

## -----------------------------------------------------------------------------
rm(list=ls()) # clear memory

