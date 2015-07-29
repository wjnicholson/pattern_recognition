#
# QUESTION 1
#

library(MASS)

#P dimensional normal density calculation, as given in spec
q1.norm.den <- function(x, mean, inv.covariance, determinant) {
  dim <- length(x)
  #Split calculation into two parts 
  firstPart <- 1 / ((2 * pi)^(dim / 2) * sqrt(determinant))
  secondPart <- -0.5 * t(x - mean) %*% inv.covariance %*% (x - mean) 
  firstPart * exp(secondPart) 
} 

#Calculate the likelihood of data x w.r.t given parameters 
q1.likelihood <- function(x, means, covariances, mix.prop, inverses) { 
  likelihood <- 0 
  #Number of parameters we are estimating 
  num.param <- length(mix.prop) 
  
  for (data.point in x) { 
    point.likelihood <- 0 
    #Sum the likelihoods for each parameter 
    for (index in seq(1, num.param)) {
      point.likelihood <- point.likelihood + mix.prop[index] * q1.norm.den(data.point, means[index,], inverses[[index]], norm(covariances[[index]]))
    }
    likelihood <- likelihood + log(point.likelihood)
  }
  likelihood
}

#Density of a finite mixture distribution
fin.mix.density <- function(x, means, covariances, mix.prop, inverses) {
  density <- 0
  num.param <- length(mix.prop) 
  for (dist.i in seq(1, num.param)) {
    density <- density + mix.prop[dist.i] * q1.norm.den(x, means[dist.i,], inverses[[dist.i]], norm(covariances[[dist.i]]))
  }
  density
}

#EM Algorithm Implementation
em.norm <- function(x, means, covariances, mix.prop) {
  #Control the algorithm. Terminates within max.iterations or
  #when the difference between two likelihoods is less than eps
  eps <- 0.01
  max.iterations <- 100

  #Specification of model
  num.points <- nrow(x)
  dimension <- ncol(x)
  num.dist <- length(mix.prop)

  #Termination parameters
  num.iteration <- 0
  change.likelihood <- 5 + eps
  likelihood <- -Inf

  #Initialise responsibility matrix
  resps <- matrix(0, nrow=num.points, ncol=num.dist)

  inverses <- covariances
  for(dist.i in seq(1, num.dist)) {
    inverses[[dist.i]] <- ginv(covariances[[dist.i]])
  }

  while (change.likelihood > eps && num.iteration < max.iterations) {
    
    #E-Step: Evaluate Responsibilites
    for (point.i in seq(1, num.points)) {
      for (dist.i in seq(1, num.dist)) {
        resps[point.i, dist.i] <- mix.prop[dist.i] * 
          q1.norm.den(x[point.i,], means[dist.i,], inverses[[dist.i]], norm(covariances[[dist.i]]))
      }
      #Normalisation, divide by sum over points
      resps[point.i,] <- resps[point.i,] / sum(resps[point.i,])
    }

    #M-Step: Compute parameters using new responsibilites
    for (dist.i in seq(1, num.dist)) {
      nk <- sum(resps[,dist.i])

      #Mean
      means[dist.i,] = rep.int(0, dimension)
      for (point.i in seq(1, num.points)) {
        means[dist.i,] <- means[dist.i,] + 
          resps[point.i, dist.i] * x[point.i,]
      }
      means[dist.i,] <- means[dist.i,] / nk

      #Covariance
      covariances[[dist.i]] <- matrix(0, nrow=dimension, ncol=dimension)
      for (point.i in seq(1, num.points)) {
        covariances[[dist.i]] <- covariances[[dist.i]] +
          resps[point.i, dist.i] * (x[point.i,] - means[dist.i,]) %*%
          (t(x[point.i,] - means[dist.i,]))
      }
      covariances[[dist.i]] = covariances[[dist.i]] / nk
      inverses[[dist.i]] <- ginv(covariances[[dist.i]])

      #Mixing Proportions
      mix.prop[dist.i] <- nk / num.points
    }

    #Evaluate likelihood
    likelihood.prev <- likelihood
    likelihood <- q1.likelihood(x, means, covariances, mix.prop, inverses)
    #Change is what we check for convergence
    change.likelihood <- abs(likelihood - likelihood.prev)

    num.iteration <- num.iteration + 1

    print(paste(
      "K ",
      toString(num.dist),
      ", iteration ",
      toString(num.iteration),
      ", likelihood ",
      toString(likelihood), sep=""))
  }

  #AIC is likelihood minus number of parameters
  AIC <- likelihood - 6 * num.dist + 1

  list(responsibilities=resps, likelihood=likelihood, AIC=AIC, means=means, mix.prop=mix.prop, cov=covariances, inverses=inverses)

}

question.one <- function() {
  x <- data.matrix(synth.te[,-3])

  dim <- 2
  num.points <- nrow(x)
  k.start <- 2
  k.end <- 6

  likelihood.max <- -Inf

  for (k in seq(k.start, k.end)) {
    #Means, computed randomly
    means <- cbind(runif(k, min=min(x[,1]), max=max(x[,1])), 
      runif(k, min=min(x[,2]), max=max(x[,2])))

    #Covariances
    covariances <- list(k)
    for (dist.i in seq(1, k)) {
      covariances[[dist.i]] <- matrix(0, nrow=dim, ncol=dim)
      for (point in x) {
        covariances[[dist.i]] <- covariances[[dist.i]] +
          (point - means[dist.i,]) %*% t(point - means[dist.i,])
      }

      #All distributions have the same responsibilities
      covariances[[dist.i]] <- covariances[[dist.i]] / (num.points * k)
      covariances[[dist.i]] <- covariances[[dist.i]] + diag(dim)
    }
    
    #Mixing proportions, initially each set to 1 / k
    mix.prop <- rep(1, k) / k

    em.run <- em.norm(x, means, covariances, mix.prop)

    print(paste(
      "K value is ",
      toString(k),
      " with likelihood ",
      toString(em.run[["likelihood"]]),
      " and AIC ",
      toString(em.run[["AIC"]]), sep=""))

    if (likelihood.max < em.run[["likelihood"]]) {
      likelihood.max <- em.run[["likelihood"]]
      em.run.max <- em.run
      print(paste(
        "Current best is ",
        "K ",
        toString(k),
        " likelihood ",
        toString(em.run[["likelihood"]]),
        " AIC ",
        toString(em.run[["AIC"]]), sep=""))
    }
  }

  xx <- seq(min(x[,1]), max(x[,1]), length=300)
  yy <- seq(min(x[,2]), max(x[,2]), length=300)
  grid <- as.matrix(expand.grid(xx, yy))
  means <- em.run.max[["means"]]
  cov <- em.run.max[["cov"]]
  mix.prop <- em.run.max[["mix.prop"]]
  inverses <- em.run.max[["inverses"]]
  grid <- apply(grid, 1, fin.mix.density, means=means, covariances=cov, mix.prop=mix.prop, inverses=inverses)
  contour(xx, yy, matrix(grid, 300))
  
}

#
# QUESTION 2
#

rrayleigh <- function (n,theta) {
  u <- runif(n,0,1)
  sqrt(-2*log(u))/sqrt(2*theta)
}

#Generate a sample of size n from Class C1 (theta1) and C2 (theta2)
generate.sample <- function(n, theta1, theta2) {
  
  sample.classes <- cbind(numeric(n), numeric(n) + 1)

  sample.data <- rbind(cbind(rrayleigh(n, theta1), rrayleigh(n, theta1)),
    cbind(rrayleigh(n, theta2), rrayleigh(n, theta2)))

  list(data=sample.data, class=sample.classes)

}

joint.density <- function(x, y, theta) {
  4 * x * y * theta * (exp(-theta * (x^2 + y^2)))
}

#Decision boundary for minimum error rate
boundary.min.error <- function(x, y, theta1, theta2) {
  (theta2 - theta1)* x^2 + (theta2 - theta1) * y^2 + log(theta1 / theta2)
}


#Returns estimated class - 1 (if in C1 returns 0)
discriminant <- function(x, y, theta1, theta2) {
  value <- joint.density(x, y, theta1) - joint.density(x, y, theta2)

  #If value > 0 returns true we assign to Class 1 = 0
  1 - as.numeric(value > 0)
}

estimate.best.theta <- function(sample) {
  #Most likely estimator
  2 * nrow(sample) / sum(apply(sample, 1, function(row) {row[1]^2 + row[2]^2}))
}

question.two <- function() {
  theta1 <- 4
  theta2 <- 2
  p1 <- 0.5
  p2 <- 0.5

  #Open plot
  pdf("pattern_q2_e.pdf")

  xs <- ys <- seq(0, 5, length=1000)

  densities <- outer(xs, ys, function(x, y) {p1 * joint.density(x, y, theta1) + p2 * joint.density(x, y, theta2)})

  image(x=xs, y=ys, z=densities)

  observations <- generate.sample(100, theta1, theta2)
  observations.data <- observations[["data"]]
  observations1 <- observations.data[1:100,]
  observations2 <- observations.data[101:200,]
  
  lines(observations1, type="p", pch="x")
  lines(observations2, type="p", pch="o")

  for (x in xs) {
    for (y in ys) {
      res <- boundary.min.error(x, y, theta1, theta2)
      if (abs(res) < 0.01) {
        lines(x=x, y=y, type="p", pch=".")
      }
    }
  }

  dev.off()

  best.theta1 <- estimate.best.theta(observations1)
  best.theta2 <- estimate.best.theta(observations2)

  predict.class <- apply(observations.data, 1, function(x) {
    discriminant(x[1], x[2], best.theta1, best.theta2)
  })

  correct <- 0
  for (i in seq(1, length(predict.class))) {
    if (predict.class[i] == observations[["class"]][i]) {
      #Predicted class correctly
      correct <- correct + 1
    }
  }

  error.rate <- 1 - correct / length(predict.class)

  print (paste(
    "Q2 g estimated parameters: theta1: ",
    toString(best.theta1),
    " theta2: ",
    toString(best.theta2),
    " with error: ",
    toString(error.rate), sep=""))


  train.n <- 200
  test.n <- 10000
  train.sample <- generate.sample(train.n, theta1, theta2)
  test.sample <- generate.sample(test.n, theta1, theta2)

  best.theta1 <- estimate.best.theta(train.sample[["data"]][1:train.n,])
  best.theta2 <- estimate.best.theta(train.sample[["data"]][(train.n + 1):(2 * train.n),])
  predict.class <- apply(test.sample[["data"]], 1, function(x) {
    discriminant(x[1], x[2], best.theta1, best.theta2)
  })

  correct <- 0
  for (i in seq(1, length(test.sample[["class"]]))) {
    if (predict.class[i] == test.sample[["class"]][i]) {
      correct <- correct + 1
    }
  }

  estimated.rate <- 1 - correct / length(test.sample[["class"]])

  #Linear discriminant analysis
  l.fit <- lda(train.sample[["data"]], train.sample[["class"]])
  l.prediction <- predict(l.fit, test.sample[["data"]])
  l.result <- table(l.prediction$class, test.sample[["class"]])

  l.error <- 1 - sum(diag(l.result)) / sum(l.result)

  #Quadratic discriminant analysis
  q.fit <- qda(train.sample[["data"]], train.sample[["class"]])
  q.prediction <- predict(q.fit, test.sample[["data"]])
  q.result <- table(q.prediction$class, test.sample[["class"]])

  q.error <- 1 - sum(diag(q.result)) / sum(q.result)

  print(paste(
    "Q2 h estimated theta1: ",
    toString(best.theta1),
    " theta2: ",
    toString(best.theta2),
    " with error: ",
    toString(estimated.rate), sep=""))

  print(paste(
    "Q2 h linear discriminant error: ",
    toString(l.error), sep=""))
  print(paste(
    "Q2 h quadratic discriminant error: ",
    toString(q.error), sep=""))

  list(train=train.sample, test=test.sample)
}

#
# QUESTION 3
#

normal.kernel <- function(z) {
  exp(-0.5 * z^2) / sqrt(2 * pi)
}

kde <- function(x.star, data, bw) {
  num.points <- nrow(data)
  dim <- length(x.star)
  score <- 0

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
q3.classify <- function(class.training, features.training, class.test, features.test, bw) {

  num.c2 <- sum(class.training)
  num.c1 <- length(class.training) - num.c2

  #Extract elements from training set belonging to each class
  features.c1 <- features.training[1 : num.c1,]
  features.c2 <- features.training[(num.c1 + 1) : nrow(features.training),] 

  num.points <- nrow(features.test)

  #Correct observations count
  correct <- 0

  for (point.i in seq(1, num.points)) {
    x <- features.test[point.i,]
    x.class <- class.test[point.i]

    score.c1 <- kde(x, features.c1, bw)
    score.c2 <- kde(x, features.c2, bw)

    if ((x.class == 0) && (score.c1 > score.c2)) {
      #Point is correctly classified as class 1
      correct <- correct + 1
    } else if ((x.class == 1) && (score.c2 > score.c1)) {
      #Point is correctly classified as class 2
      correct <- correct + 1
    }
  }

  #Final error rate
  1 - correct / nrow(features.test)
}

question.three <- function(class.training, features.training, class.test, features.test) {

  #Split training data randomly 
  split.at <- round(runif(1, 1, 100))
  d1.class <- c(rep(0, split.at), rep(1, 100 - split.at))
  d2.class <- c(rep(0, 100 - split.at), rep(1, split.at))

  d1.features <- rbind(features.training[1 : split.at,], features.training[(101 + split.at) : 200,])
  d2.features <- rbind(features.training[(split.at + 1) : 100,], features.training[101 : (100 + split.at),])

  #Loop parameters
  h.length <- 20
  hs <- seq(0.05, 1, length=h.length)

  #Checking for minimum
  min.error <- Inf

  #Matrix of average errors
  errors <- matrix(0, nrow=h.length, ncol=h.length)

  for (h.1 in seq(1, h.length)) {
    for (h.2 in seq(1, h.length)) {
      bw <- c(hs[h.1], hs[h.2])

      error.1 <- q3.classify(d1.class, d1.features, d2.class, d2.features, bw)
      error.2 <- q3.classify(d2.class, d2.features, d1.class, d1.features, bw)

      #Average the error rate
      errors[h.1, h.2] <- (error.1 + error.2) / 2

      if (errors[h.1, h.2] < min.error) {
        min.bw <- bw
        min.error <- errors[h.1,h.2]
      }
    }
  }

  print(paste(
   "Minimum error rate: ",
   toString(min.error),
   ", with bandwidth parameters: ",
    toString(min.bw), sep=""))

  #Use these best bandwidth parameters on training and test sets from Q2
  error.rate <- q3.classify(class.training, features.training, class.test, features.test, min.bw)
  print(paste(
    "With best bandwidth parameters, error rate on large test sample is: ",
    toString(error.rate), sep=""))
}

#
# QUESTION 4
#

#For question 4.a
knn.dist <- function(train, test, class, k) {

  #As given in coursework spec, weighting function
  tricube <- function(x) {
    a <- abs(x)
    if (a <= 1) {
      (1 - a^3)^3
    } else {
      0
    }
  }

  #Euclidean distance between two points
  e.dist <- function(x, y) {
    sqrt(sum((x - y)^2))
  }

  #Loop parameters
  num.test <- nrow(test)
  num.train <- nrow(train)

  class.weights <- matrix(0, nrow=num.test, ncol=ncol(test))

  for (point.i in seq(1, num.test)) {

    dists <- rep(1, num.train)
    for (train.i in seq(1, num.train)) {
      #Calculate distance from testing point to all training points
      dists[train.i] <- e.dist(train[train.i,], test[point.i,])
    }

    closest <- sort.int(dists, index.return=TRUE)
    k.closest <- closest[["ix"]][1 : k]
    k.closest.dists <- closest[["x"]][1 : k]

    #Divide by maximum element, to normalise so max is 1

    for (n.i in k.closest) {
      #Add 1 as our class is stored from 0, but R indexes from 1
      n.class <- class[n.i] + 1
      #Add the weight to our class, calculated using tricube function
      class.weights[point.i, n.class] <- class.weights[point.i, n.class] + tricube(dists[n.i] / k.closest.dists[k])
    }
  }

  class.weights

}

q4.error.rate <- function(weights, class) {
  num.points <- nrow(weights)
  correct <- 0

  for (point.i in seq(1, num.points)) {
    if ((class[point.i] == 0) && (weights[point.i, 1] > weights[point.i, 2])) {
      #Correctly classified class 1
      correct <- correct+1
    } else if ((class[point.i] == 1) && (weights[point.i, 2] >= weights[point.i, 1])) {
      #Correctly classified class 2
      correct <- correct+1
    }
  }

  #Return error rate
  1 - (correct / num.points)
}

question.four <- function(class.training, features.training, class.test, features.test) {
  #Split training data randomly, similar to Q3
  split.at <- 50
  d1.class <- c(rep(0, split.at), rep(1, 100 - split.at))
  d2.class <- c(rep(0, 100 - split.at), rep(1, split.at))

  d1.features <- rbind(features.training[1 : split.at,], features.training[(101 + split.at) : 200,])
  d2.features <- rbind(features.training[(split.at + 1) : 100,], features.training[101 : (100 + split.at),]
  )

  min.error <- Inf
  min.k <- 3

  for (k in seq(3, 39, by=4)) {
    d1.weights <- knn.dist(d2.features, d1.features, d2.class, k)
    d2.weights <- knn.dist(d1.features, d2.features, d1.class, k)

    d1.error <- q4.error.rate(d1.weights, d1.class)
    d2.error <- q4.error.rate(d2.weights, d2.class)

    #Average error rate between classes
    error.rate <- (d2.error + d1.error) / 2

    if (error.rate < min.error) {
      #Found a new minimum error
      min.k <- k
      min.error <- error.rate
    }
    print(k)
    print(error.rate)
  }

  print(paste(
    "K found to give minimum error: ",
    toString(min.k),
    ", error rate when applied to test set: ",
    toString(min.error), sep=""))

  #Using our best k, calculate error rate on sample from Q2
  weights <- knn.dist(features.training, features.test, class.training, min.k)
  error.rate <- q4.error.rate(weights, class.test)

  print(paste(
    "K found to give minimum error: ",
    toString(min.k),
    ", error rate when applied to large test set: ",
    toString(error.rate), sep=""))
}

#Call all the functions

#question.one()

data.set <- question.two()

class.training <- data.set[["train"]][["class"]]
class.test <- data.set[["test"]][["class"]]
data.training <- data.set[["train"]][["data"]]
data.test <- data.set[["test"]][["data"]]

question.three(class.training, data.training, class.test, data.test)
question.four(class.training, data.training, class.test, data.test)
