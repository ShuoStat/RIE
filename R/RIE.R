
#' Relative Information Error for continuous variables
#' 
#' @description 
#' 
#' Information error and relative information error for continuous variables following 
#' the multinomial distribution using correlation matrix.
#' 
#' @param cor covariance matrix
#' @param rvars relevant variables
#' @param svars selected variables
#' @param std logical. whether information should be standardized.
#' @author Shuo Wang
#' 
#' @examples 
#' m <- GenMatrix(3, rlist = list(c(1, 2, 0.5), c(1, 3, 0.3), c(2, 3, 0.7)))
#' rvars <- c("x1", "x2")
#' svars <- c("x2", "x3")
#' rie.c(m, rvars, svars)
#' 
#' @import Brobdingnag
#' @export rie.c


rie.c <- function(cor, rvars, svars, std = TRUE){
  
  if (!is.matrix(cor))
    stop("cor should be matrix")
  
  if (dim(cor)[1] != dim(cor)[2])
    stop("cor should have identical rows and collumns")
  
  if (any(diag(cor) != 1))
    stop("cor should be the correlation matrix")
  
  #- entropy 
  entropy.c <- function(s){
    ent <- 1 / 2 * log(2 * pi * exp(1) * s)
    return(ent)
  }
  
  #- joint entropy 
  #- Use Brobdingnag in case large values
  joint.c <- function(cor){
    if (any(dim(cor) == 0)){
      joint_ent = 0
    } else {
      m <- eigen(2 * pi * exp(1) * cor)$values
      joint_ent <- 1 / 2 * log(abs(Brobdingnag::prod(Brobdingnag::as.brob(m))))
    }
    return(joint_ent)
  }
  
  #- false positive and false negative variables
  fp <- setdiff(svars, rvars)
  tp <- intersect(rvars, svars)
  
  #- all related variables
  a <- c(fp, rvars)
  if (length(a) == 0)
    stop("No ")
  
  #- information error of the false positive
  if (length(rvars) == 1){
    i.fp <- joint.c(cor[a, a]) - entropy.c(cor[rvars, rvars, drop = FALSE])
  }else{
    i.fp <- joint.c(cor[a, a]) - joint.c(cor[rvars, rvars, drop = FALSE])
  }
  
  #- information error of the false negative
  if (length(c(tp, fp)) == 1){
    i.fn <- joint.c(cor[a, a]) - entropy.c(cor[c(tp, fp), c(tp, fp), drop = FALSE])
  } else {
    i.fn <- joint.c(cor[a, a]) - joint.c(cor[c(tp, fp), c(tp, fp), drop = FALSE])
  }
  
  #- maximum information and information error
  if (std) {
    max.inf <- as.vector(joint.c(cor) / entropy.c(1))
    ie <- as.vector((i.fp + i.fn) / entropy.c(1))
    i.fn <- i.fn / entropy.c(1)
    i.fp <- i.fp / entropy.c(1)
    
  } else {
    max.inf <- as.vector(joint.c(cor)) 
    ie <- as.vector(i.fp + i.fn)
  }
  
  #- relative information error
  rie <- ie / max.inf
  rie <- as.vector(rie)
  
  #- error rate
  error.rate <- length(c(setdiff(rvars, svars), 
                        setdiff(svars, rvars))) / ncol(cor)
  
  return(list(fp = i.fp,
              fn = i.fn,
              ie = ie,
              max = max.inf,
              rie = rie,
              error.rate = error.rate))
}

#-------------------------------------------------------------------------------

#' Relative Information Error using data sets as input
#' 
#' @description 
#' 
#' Calculate the information error and relative information error using the
#' data as the input. The data will be transformed to discretized data, from 
#' where the probability is calculated. 
#' 
#' @param X matrix or data.frame containing all variables
#' @param rvars relevant variables
#' @param svars selected variables
#' @param contivars which variables are continuous? By default, only continuous 
#' variables will be discretized. If NULL, all variables are continuous. 
#' @param method the name of the entropy estimator. We use the entropy 
#' computation of the infotheo package, where four methods are includes: 
#' "emp", "mm", "shrink", "sg". 
#' For more details, see entropy()  in "infotheo" package. 
#' @param disc the strategies to discretize data. Three options are provided: 
#' "equalfreq", "equalwidth", and "globalequalwidth". 
#' See discretize() in "infotheo" package.
#' @param nbins iteger specifying the number of bins to be used for the 
#' discretization
#' @author Shuo Wang
#' 
#' @examples 
#' 
#' m <- GenMatrix(3, rlist = list(c(1, 2, 0.5), c(1, 3, 0.3), c(2, 3, 0.7)))
#' D <- GenData(100, rmatrix = m)
#' rvars = c("x1", "x2")
#' svars = c("x2", "x3")
#' rie.d(D, rvars, svars)
#' 
#' @import infotheo
#' @export rie.d

rie.d <- function(X, rvars, svars, 
                  contivars = NULL,
                  method = "emp", 
                  disc = "equalfreq",
                  nbins = nrow(X)^(1/3)){
  
  if (!(is.matrix(X) | is.data.frame(X)))
    stop("X should be matrix or data.frame")
  
  if (!all(rvars %in% colnames(X)) | !all(svars %in% colnames(X)))
    stop("rvars and svars should be contained in the colnames of X")
  
  if (is.null(contivars))
    contivars = colnames(X)
  
  #- discretize data 
  D <- as.data.frame(X)
  D[,contivars] <- infotheo::discretize(X[,contivars], disc = disc, nbins = nbins)
  
  #- false positive and false negative variables
  fp <- setdiff(svars, rvars)
  tp <- intersect(rvars, svars)
  
  #- all relevant variables
  a <- c(fp, rvars)
  
  #- information of false positive and false negative
  i.fp <- infotheo::entropy(D[, a, drop = FALSE], method = method) - 
    infotheo::entropy(D[, rvars, drop = FALSE], method = method)
  i.fn <- infotheo::entropy(D[, a, drop = FALSE], method = method) - 
    infotheo::entropy(D[, c(tp, fp), drop = FALSE], method = method)
  
  #- maximum information 
  max.inf <- infotheo::entropy(D, method = method)
  
  #- information error
  ie <- i.fp + i.fn 
  
  #- relative information error
  rie <- ie / max.inf
  rie <- as.vector(rie)
  
  #- error rate
  error.rate <- length(c(setdiff(rvars, svars), setdiff(svars, rvars))) / ncol(D)
  
  return(list(fp = i.fp,
              fn = i.fn,
              ie = ie,
              max = max.inf,
              rie = rie,
              error.rate = error.rate))
}

#' Generate data
#' 
#' @description 
#' Simulated data using the correlation matrix
#' 
#' @param n number of the observations.
#' @param rmatrix correlation matrix
#' @param mu vector, theoretical mean of the generated data.
#' @param sd vector, the standard deviation 
#' @param empirical logical. Whether the mean and covariance matrix refers 
#' to the empirical or population. See mvrnorm() in MASS package.
#' @seealso GenMatrix
#' @author Shuo Wang
#' 
#' @examples 
#' m <- GenMatrix(3, rlist = list(c(1, 2, 0.5), c(1, 3, 0.3), c(2, 3, 0.7)))
#' D <- GenData(100, rmatrix = m)
#' head(D)
#' 
#' @import MASS
#' @export GenData

GenData <- function (n, rmatrix, mu = 0, sd = 1, empirical = FALSE) {
  
  if (!is.matrix(rmatrix))
    stop("rmatrix should be matrix")
  
  if (length(mu) == 1) {
    mu <- rep(mu, nrow(rmatrix))
  }
  
  if (length(mu) != nrow(rmatrix)) 
    stop("mu has improper length")
  
  if (length(sd) == 1) {
    sd <- rep(sd, nrow(rmatrix))
  }
  
  if (length(sd) != nrow(rmatrix)) 
    stop("se has improper length")
  rmatrix <- sd %*% t(sd) * rmatrix
  d <- MASS::mvrnorm(n, mu = mu, Sigma = rmatrix, empirical = empirical)
  colnames(d) <- paste0("x", seq_len(ncol(rmatrix)))
  return(d)
  
}


#' Generate correlation matrix
#' 
#' @description 
#' A shortcut to generate correlation matrix. The generated matrix is the 
#' input of the GenData.
#' 
#' @param p number of variables in the matrix
#' @param rlist a list to add values to matrix. The vectors in the list contain
#' three elements, with the first two being the position of names 
#' of variables, and the last one being the correlations between them.
#' @param varnames variable names. If NULL, the varnames are defaulted as 
#' x1, ..., xp.
#' @author Shuo Wang
#' 
#' @examples 
#' 
#' m <- GenMatrix(3, rlist = list(c(1, 2, 0.5), c(1, 3, 0.3), c(2, 3, 0.7)))
#' m
#' 
#' m <- GenMatrix(3, rlist = list(c("x1", "x2", 0.5), c("x1", "x3", 0.3), 
#'                                c("x2", "x3", 0.7)))
#' m
#' 
#' @export GenMatrix

GenMatrix <- function (p, rlist = list(), varnames = NULL){
  
  if(!is.list(rlist))
    warning("rlist should be a list")
  
  ma <- matrix(0, p, p)
  diag(ma) <- 1
  
  if (is.null(varnames))
    varnames <- paste0("x", 1:p)
  
  if (p != length(varnames))
    stop("p should be equal to the length of varnames")
  
  colnames(ma) <- rownames(ma) <- varnames
  
  for (i in rlist) {
    ma[i[1], i[2]] <- ma[i[2], i[1]] <- as.numeric(i[3])
  }
  
  return(ma)
}



