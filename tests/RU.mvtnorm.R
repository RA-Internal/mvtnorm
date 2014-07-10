
#testing pmvnorm
test.pmvnorm <- function() {
  sigma <- diag(3)
  sigma[1,2] <- .0
  sigma[1,3] <- .0
  sigma[2,3] <- .0
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  
  x <- pmvnorm(mean = MEAN, corr = sigma, upper = upper, lower = lower)[[1]]
  print(x)
  
  checkEquals(x,
    0.125
    )
}

# testing very large dimensions 
test.dmvnorm <- function() {
  x <- dmvnorm(rep(0, 50), mean = rep(0, 50), sigma = diag(50), log = FALSE)
  print(x)
  
  checkEquals(x,
  .3989423^50,
  tolerance = .01
)
}

# qmvnorm 
test.qmvnorm <- function() {
  x <- qmvnorm(p = .975, sigma = 1, mean = 0, tail = c("upper.tail") )[[1]]
  print(x)
  
  checkEquals(x,
    -1.96,
    tolerance = .01
  )
}

# testing pmvnorm and qmvnorm inversion equivalence    

test.qmvnorm.pmvnorm <- function() {
  x1 <- qmvnorm(p = .975, mean = 0, sigma = 1, tail = c("upper.tail")  )[[1]] 
  x2 <- qmvnorm(p = 1 - pmvnorm(mean = 0, sigma = 1, upper = Inf, lower = 1.96)[[1]], 
          mean = 0, sigma = 1, tail = c("upper.tail")  )[[1]]
  print(c(x1, x2))
  checkEquals (x1, 
               x2,
               tolerance = .01
  )  
}

# testing error messages    
test.pmvnorm.error <- function(){
  checkException(
    pmvnorm(mean = c(1,"x",4), corr = sigma, upper = upper, lower = lower)[[1]]
  ) 
}
#testing pmvt
test.pmvt <- function() {
  sigma <- diag(3)
  sigma[1,2] <- .0
  sigma[1,3] <- .0
  sigma[2,3] <- .0
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  
  x <- pmvt(delta = MEAN, corr = sigma, upper = upper, lower = lower,
       df= 100000)[[1]]
  print(x)
  checkEquals(x,
    0.125
    )
}
# testing error messages
test.pmvt.error <- function(){ 
  checkException(
  pmvt(delta = c(1,"x",4), corr = sigma, upper = upper, lower = lower, 
       df = 15)[[1]]
) 
}

# Testing Algorithm equivalence
test.algos.pmvnorm <- function() {
  sigma <- diag(3)
  sigma[1,2] <- .1
  sigma[1,3] <- .2
  sigma[2,3] <- .3
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  
  checkEquals(
    pmvnorm(mean = MEAN, corr = sigma, upper = upper, lower = lower,
            algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0) )[[1]],
    pmvnorm(mean = MEAN, corr = sigma, upper = upper, lower = lower)[[1]],
            algorithm = Miwa(steps = 128),
    tolerance = .01
  ) 
}

test.algos2.pmvnorm <- function() {
  sigma <- diag(3)
  sigma[1,2] <- .2
  sigma[1,3] <- .2
  sigma[2,3] <- .1
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  
  checkEquals(
    pmvnorm(mean = MEAN, corr = sigma, upper = upper, lower = lower,
            algorithm = TVPACK(1e-6) )[[1]],
    pmvnorm(mean = MEAN, corr = sigma, upper = upper, lower = lower)[[1]],
    algorithm = Miwa(steps = 128),
    tolerance = .01
  ) 
}

test.algos3.pmvnorm <- function() {
  sigma <- diag(3)
  sigma[1,2] <- .0
  sigma[1,3] <- .0
  sigma[2,3] <- .0
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  
  checkEquals(
    pmvnorm(mean = MEAN, corr = sigma, upper = upper, lower = lower,
            algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0) )[[1]],
    pmvnorm(mean = MEAN, corr = sigma, upper = upper, lower = lower)[[1]],
    algorithm = TVPACK(abseps = 1e-6)
  ) 
}

# testing consistency of random number generation
test.random.numb.consistency <- function() { 
  sigma <- diag(3)
  sigma[1,2] <- .3
  sigma[1,3] <- .2
  sigma[2,3] <- .3
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  
  x1 <- rmvnorm(set.seed(1000), n = 5000, mean = MEAN, sigma = sigma)[,1]
  x2 <- rmvnorm(set.seed(1000), n = 5000, mean = MEAN, sigma = sigma)[,1] 
  print(head(x1))
  
  checkEquals(x1,
              x2)
}

# testing matrix decomposition, eigen, svd
test.random.numb.consistency2 <- function() { 
  sigma <- diag(3)
  sigma[1,2] <- .3
  sigma[1,3] <- .2
  sigma[2,3] <- .3
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
    set.seed(1000)
    x <-  rmvnorm(n = 5000, mean = MEAN, sigma = sigma, method = c("svd") )
    set.seed(1000)
    x2 <- rmvnorm(n = 5000, mean = MEAN, sigma = sigma, method = c("eigen") )
  
  checkEquals(
    cor(x[,1], x2[,1]),
    1,
    tolerance = .01
    )
}

# testing matrix decomposition, chol, eigen

test.random.numb.consistency3 <- function() { 
  sigma <- diag(3)
  sigma[1,2] <- .3
  sigma[1,3] <- .2
  sigma[2,3] <- .3
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  set.seed(1000)
  x <-  rmvnorm(n = 5000, mean = MEAN, sigma = sigma, method = c("svd") )
  set.seed(1000)
  x2 <- rmvnorm(n = 5000, mean = MEAN, sigma = sigma, method = c("chol") )
  
  checkEquals(
    cor(x[,1], x2[,1]),
    1,
    tolerance = .01
  )
}

test.random.numb.consistency4 <- function() { 
  sigma <- diag(3)
  sigma[1,2] <- .3
  sigma[1,3] <- .2
  sigma[2,3] <- .3
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  set.seed(1000)
  x <-  rmvnorm(n = 5000, mean = MEAN, sigma = sigma, method = c("eigen") )
  set.seed(1000)
  x2 <- rmvnorm(n = 5000, mean = MEAN, sigma = sigma, method = c("chol") )
  
  checkEquals(
    cor(x[,1], x2[,1]),
    1,
    tolerance = .01
  )
}




































