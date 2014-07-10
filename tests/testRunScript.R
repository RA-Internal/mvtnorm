# testing mtvnorm and generating output files
# to test under 3.0.1, change "mvtnorm302.Rout" and "mtvnormTestsFor3.0.2" to end in "3.0.1"
sink("mvtnorm302.Rout", type = c("output"), split=TRUE )

options(scipen=50)
install.packages("RUnit")
library("RUnit")
install.packages("mtvnorm")
library("mvtnorm")

testsuite.mvtnorm <- defineTestSuite("mvtnorm",
                     dirs = file.path(getwd(), "tests"),
                            testFileRegexp = "^RU.+\\.R", 
                            testFuncRegexp = "^test.+"
                            )
mvtnormTest <- runTestSuite(testsuite.mvtnorm)
printTextProtocol(mvtnormTest)
mvtnormTest

sink()

printHTMLProtocol(mvtnormTest, fileName = "mtvnormTestsFor3.0.2",
                  separateFailureList = TRUE)

# testing reverse dependency of mtvnorm with package tmvtnorm (truncated normal and t 
# distributions.  
# Note:  tmvtnorm's other dependencies require R 3.0.3.

install.packages("tmvtnorm")
library("tmvtnorm")

  sigma <- diag(3)
  sigma[1,2] <- .0
  sigma[1,3] <- .0
  sigma[2,3] <- .0
  sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)]
  lower <- c(0, 0, 0)
  upper <- c( Inf, Inf, Inf)
  MEAN <- c(0, 0, 0)
  
mtmvnorm(mean = MEAN, sigma = sigma, upper = upper, lower = lower)















