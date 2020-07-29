

## construct_CovBlk #####

expect_equal(
  construct_CovBlk(3,2,3),
  matrix(c(13,9,9,9,13,9,9,9,13),nrow=3)
)
expect_equal(
  construct_CovBlk(3,2,c(1,2,3)),
  matrix(c(5,2,3,2,8,6,3,6,13),nrow=3)
)
expect_equal(
  construct_CovBlk(3,c(1,1,2),c(1,2,3)),
  matrix(c(2,2,3,2,5,6,3,6,13),nrow=3)
)

## construct_CovMat #####



# expect_equal_to_reference()
