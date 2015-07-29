library(MASS)

fda.classify <- function(train.features, train.class, test.features, test.class) {
  #Split into different classes
  c1.features <- train.features[train.class[] == 0, ]
  c2.features <- train.features[train.class[] == 1, ]

  num.c1 <- nrow(c1.features)
  num.c2 <- nrow(c2.features)

  #Means and prior probabilities
  mu1 <- apply(c1.features, 2, mean, na.rm=TRUE)
  mu2 <- apply(c2.features, 2, mean, na.rm=TRUE)
  prior1 <- num.c1 / (num.c1 + num.c2)
  prior2 <- 1 - prior1

  num.feat <- ncol(train.features)
  S1 <- matrix(0, nrow=num.feat, ncol=num.feat)
  S2 <- matrix(0, nrow=num.feat, ncol=num.feat)

  #Calculate covariance matrices for both classes
  c1.centered <- c1.features
  c1.centered[is.na(c1.centered[,])] <- 0
  c1.centered <- t(t(c1.centered) - mu1)
  S1 <- t(c1.centered) %*% c1.centered / (num.c1 - 1)
  c2.centered <- c2.features
  c2.centered[is.na(c2.centered[,])] <- 0
  c2.centered <- t(t(c2.centered) - mu2)
  S2 <- t(c2.centered) %*% c2.centered / (num.c2 - 1)


  #Scatter within matrix
  SW <- S1 + S2

  #Optimal vector for projection
  w.star <- ginv(SW) %*% (mu1 - mu2)

  #Threshold we will classify points with, hyperplane between the two means
  threshold <- 0.5 * t(mu1 + mu2) %*% w.star - log(prior2 / prior1)

  test.features.na <- test.features
  test.features.na[is.na(test.features.na[,])] <- 0
  predicted.class <- as.numeric(apply(test.features.na, 1, function(x)(as.numeric((t(w.star) %*% x) < threshold))))
  error.rate <- 1 - (sum(diag(table(predicted.class, test.class))) / length(test.class))
  return(list(error.rate=error.rate, predicted.class=predicted.class))
}

qda.classify <- function(train.features, train.class, test.features, test.class) {
  #Implementation of QDA

  #Each class features
  c1.features <- train.features[train.class[] == 0, ]
  c2.features <- train.features[train.class[] == 1, ]

  #Means
  mu1 <- apply(c1.features, 2, mean, na.rm=TRUE)
  mu2 <- apply(c2.features, 2, mean, na.rm=TRUE)

  #Number in each class, and general; used to calculate prior probabilities
  num.c1 <- nrow(c1.features)
  num.c2 <- nrow(c2.features)
  prior1 <- num.c1 / (num.c1 + num.c2)
  prior2 <- 1 - prior1

  #Calculate covariance matrices for both classes
  c1.centered <- c1.features
  c1.centered[is.na(c1.centered[,])] <- 0
  c1.centered <- t(t(c1.centered) - mu1)
  S1 <- t(c1.centered) %*% c1.centered / (num.c1 - 1)

  c2.centered <- c2.features
  c2.centered[is.na(c2.centered[,])] <- 0
  c2.centered <- t(t(c2.centered) - mu2)
  S2 <- t(c2.centered) %*% c2.centered / (num.c2 - 1)

  S1.inv <- ginv(S1)
  S2.inv <- ginv(S2)
  S1.norm <- norm(S1)
  S2.norm <- norm(S2)

  #Discriminant function, returns predicted class of data point x
  disc <- function(x) {
    f <- function(x, norm, inv, mean, prior) {
      m.cent <- x - mean
      m.cent[is.na(m.cent)] <- 0
      p1 <- -0.5 * log(norm) + log(prior)
      p2 <- -0.5 * t(m.cent) %*% inv %*% (m.cent)
      return(p1 + p2)
    }
    return(as.numeric(f(x, S1.norm, S1.inv, mu1, prior1) < f(x, S2.norm, S2.inv, mu2, prior2)))
  }

  predicted.class <- as.numeric(apply(test.features, 1, disc))
  error.rate <- 1 - (sum(diag(table(predicted.class, test.class))) / nrow(test.features))

  return(list(error.rate=error.rate, predicted.class=predicted.class))
}

knn.classify <- function(train.features, train.class, test.features, test.class, k) {

  num.test <- nrow(test.features)
  num.train <- nrow(train.features)

  scores.c1 <- rep(0, num.test)
  scores.c2 <- rep(0, num.test)
  predicted.class <- rep(0, num.test)

  #As given in coursework spec, weighting function
  tricube <- function(x) {
    a <- abs(x)
    if (a <= 1) {
      return((1 - a^3)^3)
    } else {
      return(0)
    }
  }

  for (observation.i in 1:num.test) {
    #Compute the distances to all the training points
    dists <- rep(1, num.train)
    for (train.i in 1:num.train) {
      dists[train.i] <- dist(rbind(test.features[observation.i,], train.features[train.i,]))
    }

    closest <- (sort.int(dists, index.return=TRUE))

    #Indices and distances of k closest elements
    k.indices <- closest[["ix"]][1:k]
    k.dists <- closest[["x"]][1:k]

    #Rescale distances so they are constrained to [0,1]
    k.dists <- k.dists/k.dists[k]

    for (k.i in 1:length(k.indices)) {
      #Extract class and distance of neighbor
      neighbor.i <- k.indices[k.i]
      neighbor.class <- train.class[neighbor.i]
      neighbor.dist <- k.dists[k.i]

      #Add tricube of neighbor.distance to class scores
      if (neighbor.class == 0) {
       scores.c1[observation.i] <- scores.c1[observation.i] + tricube(neighbor.dist)
      } else {
       scores.c2[observation.i] <- scores.c2[observation.i] + tricube(neighbor.dist)
      }
    }

    predicted.class[observation.i] <- as.numeric(scores.c1[observation.i] < scores.c2[observation.i])
  }

  error.rate <- 1 - (sum(diag(table(predicted.class, test.class))) / num.test)

  return(list(predicted.class=predicted.class, error.rate=error.rate))
}

perceptron <- function(train.features, train.class, weights, steps) {

  w <- rep(0, ncol(train.features))
  b <- 0
  if(missing(steps)) {
    steps <- 40
  }

  #Rate decreases as we aim to narrow down on correct point
  rate <- seq(steps, 1) / steps
  num.points <- length(train.class)
  R.const <- max(apply(train.features, 1, function(x)(sum(x * x))))
  classes <- train.class * 2 - 1

  #Perceptron algorithm, decreasing rate value should get us close to optimal value
  for(r in rate) {
    y <- apply(train.features, 1, function(x){
        ifelse(b + sum(x * w) < 0, -1, 1)
      })
    for(i in 1:length(y)) {
      if(y[i] != classes[i]) {
         w <- w - r * y[i] * train.features[i,] * weights[i]
	 b <- b - r * y[i] * R.const * weights[i]
      }
    }
  }
  s <- sqrt(sum(w*w))
  w <- w/s
  b <- b/s

  y <- apply(train.features, 1, function(x){
      ifelse(b + sum(x * w) < 0, -1, 1)
    })
  weight <- as.numeric(y * classes)

  wrong <- sum(as.numeric(y != classes) * weights)

  #Weak classifier
  h <- function(x) {
    return(ifelse(b + sum(x * w) < 0, -1, 1))
  }

  return(list(func=h, wrong=wrong, weight=weight))
}

percept.bag <- function(train.features, train.class, test.features, test.class) {
  num.points <- nrow(train.features)
  weights <- rep(1, num.points) / num.points
  p <- perceptron(train.features, train.class, weights)
  predicted.class <- (apply(test.features, 1, p$func) + 1) / 2
  error.rate <- 1 - (sum(diag(table(predicted.class, test.class))) / nrow(test.features))
  return(list(predicted.class=predicted.class, error.rate=error.rate))
}


decision.stumps <- function(train.features, train.class, weights) {

  num.feat <- ncol(train.features)
  best.error <- 1
  best.feat <- 1
  best.s.i <- 1

  stump <- function(features) {
    results <- sapply(features, function(x){
      y <- features <= x
      wrong.below <- sum(as.numeric(y == train.class) * weights)
      wrong.above <- sum(as.numeric(y != train.class) * weights)
      return(min(wrong.below, wrong.above))})
    i <- which.min(results)
    return(list(err=results[i], index=i))
  }

  stumps <- apply(train.features, 2, stump)
  f.i <- as.numeric(which.min(sapply(stumps, function(x)(x$err))))
  if(stumps[[f.i]]$err < best.error) {
    best.error <- stumps[[f.i]]$err
    best.s.i <- stumps[[f.i]]$index
    best.feat <- f.i
  }

  threshold <- train.features[best.s.i, best.feat]
  y <- train.features[,best.feat] <= threshold

  wrong.below <- sum(as.numeric(y == train.class) * weights)
  wrong.above <- sum(as.numeric(y != train.class) * weights)
  wrong <- min(wrong.above, wrong.below)

  gt <- wrong.above > wrong.below

  if(gt) {
    weight <- (train.class * 2 - 1) * (2 * as.numeric(y) - 1)
    h <- function(x) {
      return((2 * as.numeric(x[best.feat] <= threshold)) - 1)
    }
  } else {
    weight <- (train.class * 2 - 1) * (-2 * as.numeric(y) + 1)
    h <- function(x) {
      return((-2 * as.numeric(x[best.feat] <= threshold)) + 1)
    }
  }

  return(list(func=h, wrong=wrong, weight=weight))
}

adaboost <- function(train.features, train.class, T, b.percept, b.stump) {
  num.points <- nrow(train.features)
  weights <- rep(1/num.points, num.points)
  alphas <- vector()
  funcs <- list()

  if(missing(b.percept)) b.percept <- TRUE
  if(missing(b.stump)) b.stump <- TRUE

  for(t in 1:T) {
    #Check if we are cancelling using one of the types of classifier
    if(b.percept) {
      percept <- perceptron(train.features, train.class, weights)
      h <- percept
    }
    if(b.stump) {
      stump <- decision.stumps(train.features, train.class, weights)
      h <- stump
    }
    if(b.percept && b.stump) {
      if (stump$wrong < percept$wrong) {
        h <- stump
      } else {
        h <- percept
      }
    }

    #Decision function
    funcs[[t]] <- h$func
    
    wrong <- h$wrong

    #Check our best error less than 0.5
    if(wrong >= 0.5) {
      T <- t - 1
      break;
    }

    part <- h$weight

    #Calculate alpha and update weights
    a <- sum(part * weights)
    alphas[t] <- 0.5 * log((1 + a) / (1 - a))

    weights <- weights * sapply(part, function(x)(exp(-alphas[t] * x)))
    weights <- weights / sum(weights)
  }

  #Final strong classifer formed from weak classifiers
  strong.classifier <- function(x) {
    result <- sum(sapply(1:T, function(y)(alphas[y] * funcs[[y]](x))))
    return(as.numeric(result > 0))
  }

  return(strong.classifier)

}

adaboost.classify <- function(train.features, train.class, test.features, test.class, layers) {

  if(missing(layers)) layers <- 10

  H <- adaboost(train.features, train.class, layers, TRUE, TRUE)

  predicted.class <- apply(test.features, 1, H)
  error.rate <- 1 - (sum(diag(table(predicted.class, test.class))) / nrow(test.features))
  return(list(predicted.class=predicted.class, error.rate=error.rate))
}

bagging.classify <- function(train.features, train.class, test.features, test.class, n, classifier) {

  if(missing(n)) n <- 10
  predicted <- rep(0, length(test.class))
  for(i in 1:n) {
    sample.feature <- sample.int(nrow(train.features), nrow(train.features), replace=TRUE)
    predicted <- predicted + classifier(train.features[sample.feature,], 
      train.class[sample.feature], test.features, test.class)$predicted.class
  }
  predicted <- round(predicted / n) 
  error.rate <- 1 - (sum(diag(table(predicted, test.class))) / nrow(test.features))
  return(list(predicted.class=predicted, error.rate=error.rate))
}

main <- function() {
  #Final main script
  data <- read.table("data15.dat")
  data <- as.matrix(data[complete.cases(data),])

  #Shuffle data into a permutation
  data <- data[sample(nrow(data)),]
  
  class <- as.numeric(data[,1])
  features <- as.matrix(data[,2:ncol(data)])
  num.data <- nrow(data)

  #Size of the test set
  test.size <- 100
  
  test.features.indices <- sample(seq(1:num.data), test.size, replace=FALSE)

  final.test.features <- features[test.features.indices,]
  final.test.class <- class[test.features.indices]
  train.features <- features[-test.features.indices,]
  train.class <- class[-test.features.indices]


  num.subsets <- 10

  size <- floor(nrow(train.features) / num.subsets)
  acr <- matrix(1, nrow=5, ncol=num.subsets + 1)
  timing <- matrix(0, nrow=5, ncol=num.subsets + 1)
  predicted.classes <- array(0, dim=c(5, num.subsets+1, size))
  tested.classes <- matrix(0, nrow=num.subsets+1, ncol=size)

  k.start <- 11
  k.end <- 39
  k.seq <- seq(k.start, k.end, by=4)
  knn.results <- matrix(1, nrow=num.subsets+1, ncol=(length(k.seq)))

  a.start <- 5
  a.end <- 30
  a.seq <- seq(a.start, a.end, by=5)
  adaboost.results <- matrix(1, nrow=num.subsets+1, ncol=(length(a.seq)))

  for(sub.i in (0:num.subsets)) {
    #Calculate size of subset
    if (sub.i < num.subsets) {
      subset <- (sub.i * size + 1) : (sub.i * size + size)
      test.features <- train.features[subset,]
      test.class <- train.class[subset]
      rest.features <- train.features[-subset,]
      rest.class <- train.class[-subset]
      tested.classes[sub.i,] <- test.class
    } else {
      test.features <- final.test.features
      test.class <- final.test.class
      rest.features <- train.features
      rest.class <- train.class
    }


    timing[1,sub.i + 1] <- proc.time()[3]
    fda.result <- fda.classify(rest.features, rest.class, test.features, test.class)
    timing[1,sub.i + 1] <- proc.time()[3] - timing[1,sub.i + 1]

    timing[2,sub.i + 1] <- proc.time()[3]
    qda.result <- qda.classify(rest.features, rest.class, test.features, test.class)
    timing[2,sub.i + 1] <- proc.time()[3] - timing[2,sub.i + 1]

    #Remove NA values for KNN implementation, QDA and LDA library functions
    rest.complete <- complete.cases(rest.features)
    rest.features <- as.matrix(rest.features[rest.complete,])
    rest.class <- rest.class[rest.complete]
    test.complete <- complete.cases(test.features)
    test.features <- as.matrix(test.features[test.complete,])
    test.class <- test.class[test.complete]

    q.fit <- qda(rest.features, rest.class)
    q.predict <- predict(q.fit, test.features)
    q.result <- table(q.predict$class, test.class)
    q.error <- 1 - sum(diag(q.result)) / test.size

    l.fit <- lda(rest.features, rest.class)
    l.predict <- predict(l.fit, test.features)
    l.result <- table(l.predict$class, test.class)
    l.error <- 1 - sum(diag(l.result)) / test.size

    k.result <- sapply(k.seq, function(k)(list(knn.classify(rest.features, rest.class, test.features, test.class, k))))
    knn.results[sub.i + 1,] <- sapply(k.result, function(x)(x[['error.rate']]))

    timing[4,sub.i + 1] <- system.time(knn.classify(rest.features, rest.class, test.features, test.class, k.start))[3]
    
    timing[3,sub.i + 1] <- proc.time()[3]
    ada.result <- adaboost.classify(rest.features, rest.class, test.features, test.class, 20)
    timing[3,sub.i + 1] <- proc.time()[3] - timing[3,sub.i + 1]

    a.results <- sapply(a.seq, function(a)(list(adaboost.classify(rest.features, rest.class, test.features, test.class, a))))
    adaboost.results[sub.i + 1,] <- sapply(a.results, function(x)(x[['error.rate']]))

    acr[1, sub.i + 1] <- fda.result[["error.rate"]]
    acr[2, sub.i + 1] <- l.error
    acr[3, sub.i + 1] <- qda.result[["error.rate"]]
    acr[4, sub.i + 1] <- q.error
    acr[5, sub.i + 1] <- ada.result[["error.rate"]]

    if (sub.i < num.subsets) {
      predicted.classes[1, sub.i + 1,] <- fda.result[["predicted.class"]]
      predicted.classes[2, sub.i + 1,] <- l.predict$class
      predicted.classes[3, sub.i + 1,] <- qda.result[["predicted.class"]]
      predicted.classes[4, sub.i + 1,] <- q.predict$class
      predicted.classes[5, sub.i + 1,] <- ada.result[["predicted.class"]]
    }
  }


  #McNemar
  mcnemar <- table(qda.result$predicted.class, fda.result$predicted.class)
  mcnemar.result <- mcnemar.test(mcnemar)[3]$p.value
  significance <- 0.01
  if(mcnemar.result < significance) {
    print(paste("We reject the null hypothesis with significance ",
      toString(significance),
      " with probability ",
      toString(mcnemar.result),
      sep=""))
  } else {
    print(paste("We accept the null hypothesis with significance ",
      toString(significance),
      " due to probability ",
      toString(mcnemar.result),
      sep=""))
  }
  
  


  print.stuff <- function() {

    #Graph showing error rate over test rate
    jpeg("errors.jpg")


    plot(rep(1,10), acr[1,1:10], xlab="", ylab="", xaxt="n", pch=4, xlim=c(1,4), ylim=c(0,0.3))
    z <- c(1, 0, 2, 0, 3)
    for(i in c(1, 3, 5)) {
      points(rep(z[i], 10), acr[i,1:10], xlab="", ylab="", xaxt="n", pch=4)
      points(z[i], mean(acr[i,1:10]), xlab="", col="red", ylab="", xaxt="n", pch=16)
    }
    points(rep(4, 10), knn.results[1:10, 8], xlab="", ylab="", xaxt="n", pch=4)
    points(4, mean(knn.results[1:10, 8]), xlab="", col="red", ylab="", xaxt="n", pch=16)
    axis(1, at=1:4, labels=c("FDA", "QDA", "Adaboost", "KNN"), las=0)
    for(i in c(1, 3, 5)) {
      text(z[i] + 0.2, mean(acr[i,1:10]), labels=signif(mean(acr[i,1:10]),4), cex=0.9)
    }
    text(3.8, mean(knn.results[1:10, 8]), labels=signif(mean(knn.results[1:10, 8]),4), cex=0.9)



    dev.off()


    #Graph showing KNN stuff
    jpeg("knn.jpg")


    plot(rep(1,10), knn.results[1:10,1], xlab="", ylab="", xaxt="n", pch=4, xlim=c(1,8), ylim=c(0,0.4))
    for(i in 1:8) {
      points(rep(i,10), knn.results[1:10,i], xlab="", ylab="", xaxt="n", pch=4)
      points(i, mean(knn.results[1:10,i]), xlab="", col="red", ylab="", xaxt="n", pch=16)
    }
    axis(1, at=1:8, labels=k.seq, las=0)

    dev.off()



    #Graph showing Adaboost stuff
    jpeg("adaboost.jpg")


    plot(rep(1,10), adaboost.results[1:10,1], xlab="", ylab="", xaxt="n", pch=4, xlim=c(1,6), ylim=c(0,0.3))
    for(i in 1:length(a.seq)) {
      points(rep(i,10), adaboost.results[1:10,i], xlab="", ylab="", xaxt="n", pch=4)
      points(i, mean(adaboost.results[1:10,i]), xlab="", col="red", ylab="", xaxt="n", pch=16)
    }
    axis(1, at=1:length(a.seq), labels=a.seq, las=0)
    for(i in 1:5) {
      text(i + 0.3, mean(adaboost.results[1:10, i]), labels=signif(mean(adaboost.results[1:10, i]),4), cex=0.9)
    }
    text(5.7, mean(adaboost.results[1:10, 6]), labels=signif(mean(adaboost.results[1:10, 6]),4), cex=0.9)

    dev.off()


    #Graph showing Bagging
    jpeg("bagging.jpg")


    plot(1, bagging.classify(train.features, train.class, test.features, test.class, 1, percept.bag)$error.rate, xlab="", ylab="", xaxt="n", pch=4, xlim=c(1,10), ylim=c(0,0.4))
    for(i in 2:10) {
      points(i, bagging.classify(train.features, train.class, test.features, test.class, i, percept.bag)$error.rate, xlab="", ylab="", xaxt="n", pch=4)
    }
    axis(1, at=1:10, labels=1:10, las=0)
    title(main="Graph of Bagging implementation performance", xlab="Number of resampled subsets", ylab="Error Rate")


    dev.off()

    jpeg("bagging2.jpg")


    plot(1, bagging.classify(train.features, train.class, test.features, test.class, 1, qda.classify)$error.rate, xlab="", ylab="", xaxt="n", pch=4, xlim=c(1,10), ylim=c(0,0.4))
    for(i in 2:10) {
      points(i, bagging.classify(train.features, train.class, test.features, test.class, i, qda.classify)$error.rate, xlab="", ylab="", xaxt="n", pch=4)
    }
    axis(1, at=1:10, labels=1:10, las=0)
    title(main="Graph of Bagging implementation performance", xlab="Number of resampled subsets", ylab="Error Rate")


    dev.off()

    

    #Write table of timings

    z <- ""
    for(j in 1:10) {
      for(i in 1:4) {
        z <- paste(z, signif(timing[i,j], 4), sep=" & ")
      }
      z <- paste(z, " \\\\\n ", sep="")
    }
    write(z, file="timing.txt")



    z <- ""
    for(j in 1:10) {
      for(i in c(1,3,5)) {
        z <- paste(z, signif(acr[i,j], 4), sep=" & ")
      }
      z <- paste(z, signif(knn.results[j,8], 4), sep=" & ")
      z <- paste(z, " \\\\\n ", sep="")
    }
    for(i in c(1,3,5)) {
      z <- paste(z, signif(mean(acr[i,1:10]), 4), sep=" & ")
    }
    z <- paste(z, signif(mean(knn.results[1:10,8]), 4), sep=" & ")
    z <- paste(z, " \\\\\n ", sep="")

    write(z, file="errors.txt")


    z <- ""
    for(j in 1:10) {
      z <- paste(z, signif(knn.results[j,8], 4), sep=" & ")
      z <- paste(z, " \\\\\n ", sep="")
    }
    for(i in c(1,3,5)) {
      z <- paste(z, signif(mean(acr[i,1:10]), 4), sep=" & ")
    }
    z <- paste(z, signif(mean(knn.results[1:10,8]), 4), sep=" & ")
    z <- paste(z, " \\\\\n ", sep="")

    write(z, file="errors.txt")





    z <- ""
    h <- list()
    for(i in a.seq) {
	time <- proc.time()
	h[[i]] <- adaboost(train.features, train.class, i)
	time <- proc.time() - time
        z <- paste(z, signif(time[3], 4), sep=" & ")
    }
    z <- paste(z, " \\\\\n ", sep="")
    test <- function(x) {
      predicted.class <- apply(test.features, 1, x)
      error.rate <- 1 - (sum(diag(table(predicted.class, test.class))) / nrow(test.features))
      return(list(predicted.class=predicted.class, error.rate=error.rate))
    }
    for(i in a.seq) {
        z <- paste(z, signif(system.time(test(h[[i]]))[3], 4), sep=" & ")
    }
    z <- paste(z, " \\\\\n ", sep="")
    write(z, file="adaboost_timing.txt")


    #Boring console printing
    print(paste(
      "Implementation of FDA error rate: ",
      toString(acr[1,]),
      " with average: ",
      toString(mean(acr[1,])),
      sep=""))

    print(paste(
      "Library LDA error rate: ",
      toString(acr[2,]),
      " with average: ",
      toString(mean(acr[2,])),
      sep=""))

    print(paste(
      "Implementation of QDA error rate: ",
      toString(acr[3,]),
      " with average: ",
      toString(mean(acr[3,])),
      sep=""))

    print(paste(
      "Library QDA error rate: ",
      toString(acr[4,]),
      " with average: ",
      toString(mean(acr[4,])),
      sep=""))

    print(paste(
      "Adaboost error rate: ",
      toString(acr[5,]),
      " with average: ",
      toString(mean(acr[5,])),
      sep=""))

    print(paste(
      "KNN error rate: ",
      toString(knn.results),
      sep=""))
  }

  #UNCOMMENT
  print.stuff()
}

main()
