library(MASS)


kde <- function(x.star, data, bw) {
  num.points <- nrow(data)
  dim <- length(x.star)
  score <- 0

  normal.kernel <- function(z) {
    exp(-0.5 * z^2) / sqrt(2 * pi)
  }

  for(point.i in seq(1, num.points)) {
    temp <- 1
    for(feature.i in seq(1, dim)) {
      temp <- temp * normal.kernel((x.star - data[point.i,])[feature.i] / bw[feature.i])
    }
    score <- score + temp
  }

  #Normalize by bandwith
  band <- prod(bw)

  score / (num.points * band)
}

#Two classes assumed for classification
kde.classify <- function(train.features, train.class, test.features, test.class) {
  #Correct observations count
  correct <- 0

  classify <- function(train.f, train.c, test.f, test.c) {

    train.f.c1 <- train.f[train.c[train.f] == 0]
    train.f.c2 <- train.f[train.c[train.f] == 1]

    for (observation.i in seq(1, nrow(test.f))) {
      x <- test.f[observation.i,]
      x.class <- test.c[observation.i]

      score.c1 <- kde(x, train.f.c1, bw)
      score.c2 <- kde(x, train.f.c2, bw)
      correct <- correct + as.numeric(x.class == as.numeric(score.c2 > score.c1))
    }

    return(1 - correct / nrow(test.f))
  }

  #Loop parameters
  h.length <- 20
  hs <- seq(0.05, 1, length=h.length)

  #Checking for minimum
  min.error <- Inf

  #Matrix of average errors
  errors <- matrix(0, nrow=h.length, ncol=h.length)
  sample1 <- sample(1:

  for (h.1 in seq(1, h.length)) {
    for (h.2 in seq(1, h.length)) {
      bw <- c(hs[h.1], hs[h.2])

      error.1 <- classify(d1.class, d1.features, d2.class, d2.features, bw)
      error.2 <- classify(d2.class, d2.features, d1.class, d1.features, bw)

      #Average the error rate
      errors[h.1, h.2] <- (error.1 + error.2) / 2

      if (errors[h.1, h.2] < min.error) {
        min.bw <- bw
        min.error <- errors[h.1,h.2]
      }
    }
  }

}

em.norm <- function(x, mean, inv, det) {
  dim <- length(x)

  #Split calculation into two parts 
  firstPart <- 1 / ((2 * pi)^(dim / 2) * sqrt(det))
  secondPart <- -0.5 * t(x - mean) %*% inv %*% (x - mean) 

  return(firstPart * exp(secondPart))
}

point.likelihood <- function(point, means, mix, invs, norms) {
  k <- length(mix)
  #Calculate probability of being generated by each mixture
  probs <- mapply(em.norm, rep(list(point), k), split(means, row(means)), invs, norms) 
  #Each probability multiplied by the mixing probs
  return(sum(mix * probs))
}

em.likelihood <- function(data, means, covs, mix, invs, norms) {
  k <- length(mix)
  likelihood <- 0

  probs <- apply(data, 1, point.likelihood, means=means, mix=mix, invs=invs, norms=norms)
  #Sum the log likelihoods of each point
  likelihood <- sum(sapply(probs, log))

  return(likelihood)
}

inv <- function(cov) {
  tryCatch(
  {
    ginv(cov)
  },
  error=function(cond) {
    ginv(0.8 * cov + 0.2 * diag(dimension))
  }, finally={
  })
}

em <- function(data, means, covs, mix) {
  eps <- 0.01
  max.iterations <- 100

  #Specification of model
  num.points <- nrow(data)
  dimension <- ncol(data)
  num.dist <- length(mix)

  #Termination parameters
  num.iteration <- 0
  change.likelihood <- 5 + eps
  likelihood <- -Inf

  #Initialise responsibility matrix
  resps <- matrix(0, nrow=num.points, ncol=num.dist)

  invs <- lapply(covs, inv)
  norms <- lapply(covs, norm)

  while (change.likelihood > eps && num.iteration < max.iterations) {
    #E-Step: Evaluate Responsibilites

    resps <- t(apply(data, 1, function(x)(mix * mapply(em.norm, rep(list(x), num.dist), split(means, row(means)), invs, norms))))
    resps <- t(apply(resps, 1, function(x)(x / sum(x))))

    #M-Step: Compute parameters using new responsibilites
    n <- apply(resps, 2, sum)

    means <- t(apply(resps, 2, function(x)(apply((x * data), 2, sum))))
    means <- means / n

    covs <- lapply(1:num.dist, function(dist.i) {
          mean.centered.data <- t(t(data) - means[dist.i,])
          return (t(mean.centered.data)%*%(mean.centered.data)/(num.points - 1))
        })

    invs <- lapply(covs, inv)        
    norms <- lapply(covs, norm)

    mix <- n / num.points

    likelihood.prev <- likelihood
    likelihood <- em.likelihood(data, means, covs, mix, invs, norms)
    change.likelihood <- abs(likelihood - likelihood.prev)

    num.iteration <- num.iteration + 1

    print(paste(
    "K ",
    toString(num.dist),
    ", iter ",
    toString(num.iteration),
    ", ll ",
    toString(likelihood), sep=""))
    if(is.infinite(likelihood)) {
      return(list(resps=resps, likelihood=likelihood.prev, means=means, mix=mix, covs=covs))
    }
  }

  print(paste(
  "Finished with K ",
  toString(num.dist),
  ", on iter ",
  toString(num.iteration),
  ", ll ",
  toString(likelihood), sep=""))

  return(list(resps=resps, likelihood=likelihood, means=means, mix=mix, covs=covs))
}


em.classify <- function(train.features, train.class, test.features, test.class) {
  c1.features <- train.features[train.class[] == 0, ]
  c2.features <- train.features[train.class[] == 1, ]
  num.c1 <- nrow(c1.features)
  num.c2 <- nrow(c2.features)

  means <- rbind(apply(c1.features, 2, mean), apply(c2.features, 2, mean))

  #Covariances
  c1.features.centered <- t(t(c1.features) - means[1,])
  c2.features.centered <- t(t(c2.features) - means[2,])
  c1.cov <- (t(c1.features.centered) %*% (c1.features.centered) / (num.c1 - 1))
  c2.cov <- (t(c2.features.centered) %*% (c2.features.centered) / (num.c2 - 1))

  c1.inv <- inv(c1.cov)
  c2.inv <- inv(c2.cov)
  c1.norm <- norm(c1.cov)
  c2.norm <- norm(c2.cov)
  prior1 <- num.c1 / (num.c1 + num.c2)
  prior2 <- 1 - prior1

  class <- function(x) {
    e1 <- prior1 * em.norm(x, means[1,], c1.inv, c1.norm)
    e2 <- prior2 * em.norm(x, means[2,], c2.inv, c2.norm)
    return(as.numeric(e2 > e1))
  }

  predicted.class <- apply(test.features, 1, class)
  error.rate <- 1 - (sum(diag(table(test.class, predicted.class))) / length(test.class))

  return(list(error.rate=error.rate, predicted.class=predicted.class))

}

old.em.classify <- function(x) {
  #Covariances
  c1.features.centered <- c1.features - means[1,]
  c2.features.centered <- c2.features - means[2,]
  c1.cov <- (t(c1.features.centered) %*% (c1.features.centered) / (num.c1 - 1))
  c2.cov <- (t(c2.features.centered) %*% (c2.features.centered) / (num.c2 - 1))

  em.run <- em(train.features, means, list(c1.cov, c2.cov), rep(0.5, 2))
  c1.inv <- inv(em.run[["covs"]][[1]])
  c2.inv <- inv(em.run[["covs"]][[2]])
  c1.norm <- norm(em.run[["covs"]][[1]])
  c2.norm <- norm(em.run[["covs"]][[2]])
  means <- em.run[["means"]]

  class <- function(x) {
    e1 <- em.norm(x, means[1,], c1.inv, c1.norm)
    e2 <- em.norm(x, means[2,], c2.inv, c2.norm)
    return(as.numeric(e2 > e1))
  }
}

data <- read.table("data15.dat")
data <- as.matrix(data[complete.cases(data),])


l.fit <- lda(data[,2:ncol(data)], as.numeric(data[,1]))
l.predict <- predict(l.fit, data[,2:ncol(data)])
l.result <- table(l.predict$class, data[,1])
print( 1 - sum(diag(l.result)) / sum(l.result))

test.em(data, 2)

#princomp(data)