library(testthat)
library(RIE)

#- test GenMatrix

test_that("GenMatrix",{
  
  m <- matrix(0, 3, 3, dimnames = list(c("x1", "x2", "x3"),
                                       c("x1", "x2", "x3")))
  diag(m) <- 1
  m[1, 2] <- m[2, 1] <- 0.5
  m[1, 3] <- m[3, 1] <- 0.3
  m[2, 3] <- m[3, 2] <- 0.7
  expect_identical(m, GenMatrix(p = 3, rlist = list(c(1, 2, 0.5),
                                                    c(1, 3, 0.3),
                                                    c(2, 3, 0.7))))
  
  expect_identical(m, GenMatrix(p = 3, rlist = list(c("x1", "x2", 0.5),
                                                    c("x1", "x3", 0.3),
                                                    c("x2", 'x3', 0.7))))
  
  expect_identical(m, GenMatrix(p = 3, rlist = list(c(1, 2, 0.5),
                                                    c(1, 3, 0.3),
                                                    c(2, 3, 0.7)), 
                                varnames = c("x1", "x2", "x3")))

})



#- test GenMatrix

test_that("GenData",{
  
  m <- GenMatrix(p = 3, rlist = list(c(1, 2, 0.5),  c(1, 3, 0.3), c(2, 3, 0.7)))
  d <- GenData(100, m, empirical = T)
  expm <- rep(0, 3)
  colm <- round(as.numeric(colMeans(d)))
  
  expsd <- rep(1, 3)
  colsd <- round(as.numeric(apply(d, 2, sd)))

  expect_identical(expm, colm)
  expect_identical(expsd, colsd)
  
})

test_that("mu and sd in GenData",{
  
  m <- GenMatrix(p = 3, rlist = list(c(1, 2, 0.5),  c(1, 3, 0.3), c(2, 3, 0.7)))
  d <- GenData(100, m, mu = 1, sd = 2, empirical = T)
  expm <- rep(1, 3)
  colm <- round(as.numeric(colMeans(d)))
  
  expsd <- rep(2, 3)
  colsd <- round(as.numeric(apply(d, 2, sd)))

  expect_identical(expm, colm)
  expect_identical(expsd, colsd)
  
})


#- test rie.c 

test_that("rie.c = error rate under equal correlations", {
  
  m <- GenMatrix(p = 3)
  rie <- rie.c(m,  rvars = c("x1", "x2"), svars = c("x1", "x3"))
  expect_equal(2/3, rie$rie)
  expect_equal(3, rie$max)
  expect_equal(1, rie$fp)
  expect_equal(1, rie$fn)
  
})


test_that("rie.c if no false selection", {
  
  m <- GenMatrix(p = 3)
  rie <- rie.c(m,  rvars = c("x1", "x2"), svars = c("x1", "x2"))
  expect_equal(0, rie$rie)
  expect_equal(3, rie$max)
  expect_equal(0, rie$fp)
  expect_equal(0, rie$fn)
  
})

test_that("rie.c if rvars and svars are c()", {
  
  m <- GenMatrix(p = 3)
  rie <- rie.c(m,  rvars = c(), svars = c())
  expect_equal(0, rie$rie)
  expect_equal(3, rie$max)
  expect_equal(0, rie$fp)
  expect_equal(0, rie$fn)
  
})

test_that("general checks", {
  
  m <- GenMatrix(p = 3)
  expect_error(rie.c(as.data.frame(m),  rvars = c(), svars = c()), 
               "cor should be matrix")
  expect_error(rie.c(m[1:2,],  rvars = c(), svars = c()), 
               "cor should have identical rows and collumns")
  m[1, 1] <- 3
  expect_error(rie.c(m,  rvars = c(), svars = c()), 
               "cor should be the correlation matrix")
})

#- test rie.d

test_that("rie less than error rate if the false selections have more correlations", 
          {
            m <- GenMatrix(p = 3, rlist = list(c(2, 3, 0.5)))
            d <- GenData(1000, m, empirical = T)
            rie <- rie.d(d,  rvars = c("x1", "x2"), svars = c("x1", "x3"))
            expect_gt(2/3, rie$rie)
          })


test_that("rie great than error rate if the false selections have less correlations", 
          {
            m <- GenMatrix(p = 4, rlist = list(c(1, 4, 0.7)))
            d <- GenData(1000, m, empirical = T)
            rie <- rie.d(d,  rvars = c("x1", "x2"), svars = c("x1", "x3"))
            expect_lt(2/4, rie$rie)
          })


test_that("rie.d if no false selection", 
          {
            m <- GenMatrix(p = 4, rlist = list(c(1, 4, 0.7)))
            d <- GenData(1000, m, empirical = T)
            rie <- rie.d(d,  rvars = c("x1", "x2"), svars = c("x1", "x2"))
            expect_equal(0, rie$rie)
            expect_equal(0, rie$fp)
            expect_equal(0, rie$fn)
          })


test_that("rie.d if rvars and svars are c()", 
          {
            m <- GenMatrix(p = 4, rlist = list(c(1, 4, 0.7)))
            d <- GenData(1000, m, empirical = T)
            rie <- rie.d(d,  rvars = c(), svars = c())
            expect_equal(0, rie$rie)
            expect_equal(0, rie$fp)
            expect_equal(0, rie$fn)
          })


test_that("If data contains binary and continuous variables", 
          {
            m <- GenMatrix(p = 4, rlist = list(c(1, 4, 0.7)))
            d <- GenData(1000, m, empirical = T)
            d <- cbind(d, x5 = rbinom(1000, 1, 0.6))
            rie <- rie.d(d,  rvars = c("x1", "x2"), svars = c("x1", "x2", "x3"),
                         contivars = c("x1", "x2", "x3", "x4"))
            expect_lt(0, rie$rie)
            expect_lt(0, rie$fp)
            expect_equal(0, rie$fn)
          })

test_that("general checks", 
          {
            m <- GenMatrix(p = 4, rlist = list(c(1, 4, 0.7)))
            d <- GenData(1000, m, empirical = T)
            expect_error(rie.d(d,  rvars = c("x1", "x2"), svars = c("x5")),
                         "rvars and svars should be contained in the colnames of X")
            
            expect_error(rie.d(as.list(d),  rvars = c("x1", "x2"), svars = c("x5")),
                         "X should be matrix or data.frame")
            
          })








