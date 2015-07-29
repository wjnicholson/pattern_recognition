library(MASS)

"one" <-
function()
{
sigma <- c(0.5, 1, 1.5)
jj <- 400
jj1 <- seq(0, 4, length=jj)
res <- rep(jj1,length(sigma))
m <- 1:jj

for(i in 1:length(sigma)){
  res[(m + (i-1)  * jj)] <- one_b(jj1, 4, sigma[i])
}
plot(jj1, res[1:jj], xlab="Threshold", ylab="Error Rate", type="l", col="red")
lines(jj1, res[(jj + 1):(2*jj)], type="l", col="blue")
lines(jj1, res[(2*jj + 1):(3*jj)], type="l", col="green")

}

"one_b" <-
function(l, s1, s2)
{
r <- numeric(length(l))
for(i in 1 : length(l)){
  r[i] <- 0.5 * (1 - exp((-l[i]^2) / (2 * s1^2)) + exp((-l[i]^2) / (2 * s2^2)))
}
r
}

"one_d" <-
function()
{
x <- sqrt((-32/15)*log(1/4))
print(x)
r <- one_b(x, 4, 1)
print(r)

  "class_one" <-
  function(x, sig)
  {
  (x / sig^2) * exp(- (x^2)/(sig^2))
  }

jj <- 1000
jj1 <- seq(0, 10, length=jj)
c1 <- 0.5 * class_one(jj1, 4)
c2 <- 0.5 * class_one(jj1, 1)
plot(jj1,c2,type="l", col="red")
lines(jj1, c1, type="l", col="blue")
abline(v=x)

ev <- c1 + c2
post1 <- c1 / ev
post2 <- c2 / ev

plot(jj1,post1,type="l", col="red", ylim=c(0,1))
lines(jj1, post2, type="l", col="blue")
abline(v=x)

}


"gi" <-
function(x, cov, mu, p)
{
inv = solve(cov)
w <- t(inv %*% mu)
b <- -0.5 * t(mu) %*% inv %*% mu - 0.5 * log(det(cov)) + log(p)
res <- - 0.5 * t(x) %*% inv %*% x + w %*% x + b

res
}

"disc" <-
function(x)
{
c1 <- matrix(c(1,0.2,0.2,1), nrow=2)
c2 <- matrix(c(1,-0.7,-0.7,1), nrow=2)
mu1 <- matrix(c(-1,-1), nrow=2)
mu2 <- matrix(c(1,1), nrow=2)

g1 <-  gi(x,c1,mu1,0.5)
g2 <-  gi(x,c2,mu2,0.5)

g1 - g2
}

"two" <-
function()
{
n <- 120
xx <- seq(-6, 6, length=n)

c1 <- matrix(c(1,0.2,0.2,1), nrow=2)
c2 <- matrix(c(1,-0.7,-0.7,1), nrow=2)
mu1 <- matrix(c(-1,-1), nrow=2)
mu2 <- matrix(c(1,1), nrow=2)

x_star <- matrix(runif(2,-6,6), nrow=2)
res = 1:2
res[1] <- gi(x_star,c1,mu1,0.5)
res[2] <- gi(x_star,c2,mu2,0.5)

gr <- as.matrix(expand.grid(xx,xx))
score.diff <- apply(gr,1,disc)
par(pty="s")
contour(xx,xx,matrix(score.diff,n))
text(x_star[1], x_star[2], "X")
}

"four" <-
function()
{
golden(test, -4, 0, 0.000001)
}

"test" <-
function(x)
{
exp(-x) * sin(x) + sinh(x)
}

"golden" <-
function(func, a, b, e)
{
golden <- 2 - (1 + sqrt(5)) / 2
lower <- a
upper <- b
c <- lower + golden * (upper - lower)
d <- upper - golden * (upper - lower)
f1 <- func(c)
f2 <- func(d)
while(abs(upper - lower) > e) {
  if(f1 < f2) { 
    upper <- d
    d <- c
    c <- lower + golden * (upper - lower)
    f2 <- f1
    f1 <- func(c)
  } else {
    lower <- c
    c <- d
    d <- upper - golden * (upper - lower)
    f1 <- f2
    f2 <- func(d)
  }
}
min(func(c),func(d))
}

"three" <-
function(end=1000, pend=15)
{
p <- sample(10:pend, 1)
pc1 <- 0.25
pc2 <- 0.75
E1 <- diag(p)
E2 <- diag(p) * 0.5
E2[1,3] <- -0.3
E2[3,1] <- -0.3
E2[2,3] <- 0.2
E2[3,2] <- 0.2
E2[4,8] <- 0.1
E2[8,4] <- 0.1
mu1 <- numeric(p)
mu2 <- runif(p, 0, 1)
ns <- seq(100, end, by=100)
qres = numeric(length(ns))
lres = numeric(length(ns))

for(i in 1:length(ns)){
  tr.sample <- matrix(0, nrow=ns[i], ncol=p)
  split = (ns[i]*pc1)
  tr.sample1 <- t(t(mvrnorm(split, mu1, E1)))
  tr.sample2 <- t(t(mvrnorm(ns[i] - split, mu2, E2)))
  tr.sample <- rbind(tr.sample1, tr.sample2)
  tr.cl = numeric(ns[i])
  tr.cl[1:split] <- 1
  tr.cl[(split+1):ns[i]] <- 2

  test.sample1 <- mvrnorm(1000, mu1, E1)
  test.sample2 <- mvrnorm(3000, mu2, E2)
  test.sample <- rbind(test.sample1, test.sample2)
  test.cl <- c(numeric(1000) + 1, numeric(3000) + 2)
  lda <- lda(tr.sample, tr.cl)
  qda <- qda(tr.sample, tr.cl)

  test.l <- predict(lda, test.sample)
  test.q <- predict(qda, test.sample)

  qres[i] <- 1 - sum(diag(table(test.q$class, test.cl)))/4000
  lres[i] <- 1 - sum(diag(table(test.l$class, test.cl)))/4000
}

plot(ns, lres, xlab="Training data size", ylab="Error rate", type="l", col="red", ylim=c(0,0.2))
lines(ns, qres, type="l", col="blue")

print(min(qres, lres))
print(qres)
print(lres)
}
