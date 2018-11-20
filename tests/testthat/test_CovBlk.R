library(SteppedPower)
library(Matrix)

context("Check CovBlk()")

expect_equal(construct_CovBlk(3,2,3),
             matrix(c(13,9,9,9,13,9,9,9,13),nrow=3))

expect_equal(
  construct_CovBlk(3,2,c(1,2,3))
  )
