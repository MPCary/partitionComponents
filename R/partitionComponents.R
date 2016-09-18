#' Partition a component matrix
#' 
#' Given a numeric matrix, returns a matrix with the same dimensions, row and 
#' column names but with only 1, 0, and -1 as values.  These values indicate the
#' partitioning of each column of the matrix into positive, neither, and 
#' negative subsets.
#' 
#' @usage
#' partition(x, FUN, ...)
#'   
#' partition.functions
#' # c("fixed", "ann")
#' @param x A numeric matrix
#' @param FUN Function to use to determine threshold values for partitioning x
#' @param ... Optional arguments to FUN
#' @details Built-in partition functions are listed in the 
#'   \code{partition.functions} vector.  The ("\code{fixed}") function uses a
#'   fixed value (specified by an additional parameter, \code{t}) to create
#'   partitions, while the ("\code{ann}") function uses an artificial neural
#'   network model to create custom partition thresholds for each component. 
#'   Additional user-defined functions are allowed.  These should return a
#'   matrix with two rows and as many columns as in \code{x}; the first row of
#'   this matrix specifies the upper partition threshold for each column of
#'   \code{x}, while the second row specifies the lower partition threshold.
#' @return A matrix comprising values of 1, 0, and -1 to indicate column-wise 
#'   partitioning of the input matrix into positive, neither, and negative sets,
#'   respectively.
#' @examples
#' x = matrix(rnorm(9), 3, 3)
#' partition(x, FUN = ".fixed", t = 0.5)
partition <- function(x, FUN, ...) {
  x = .checkMatrix(x) # Check x for well-formedness

  thresholds = get(FUN)(x, ...) # Run FUN to obtain thresholds for each column
  d = dim(thresholds)
  if((d[1] != 2) | (d[2] != ncol(x))) {
    stop("Partition threshold function failed")
  } # Check thresholds
  
  x.p = x
  x.p[] = 0
  upper = t(apply(x, 1, FUN = function(x, y){x >= y}, y = thresholds[1, ]))
  lower = t(apply(x, 1, FUN = function(x, y){x <= y}, y = thresholds[2, ]))
  x.p[upper] = 1
  x.p[lower] = -1
  
  return(x.p)
}



################################
# Wrapper function definitions #

#' Wrapper function for fixed threshold partitioning
#' 
#' Wrapper function for \code{\link{partition}} for fixed threshold partitioning
#' 
#' @param x A numeric matrix
#' @param t Threshold value for use in partitioning.  Values in x >= t will
#'   become 1 in the return matrix, values in x <= -t will become -1, and values
#'   in x between -t and t will become 0
#' @return A matrix comprising values of 1, 0, and -1 to indicate column-wise 
#'   partitioning of the input matrix into positive, neither, and negative sets,
#'   respectively.
#' @examples
#' x = matrix(rnorm(9), 3, 3)
#' fixed.partition(x)
fixed.partition <- function(x, t = 3) {
  # Wrapper function for partition() using the fixed threshold method
  partition(x, FUN = ".fixed", t = t)
}


#' Wrapper function for ANN partitioning
#' 
#' Wrapper function for \code{\link{partition}} for partitioning using an
#' artificial neural network
#' 
#' @param x A numeric matrix
#' @param model ANN model to use for partitioning
#' @return A matrix comprising values of 1, 0, and -1 to indicate column-wise 
#'   partitioning of the input matrix into positive, neither, and negative sets,
#'   respectively.
#' @examples
#' x = matrix(rnorm(9), 3, 3)
#' ann.partition(x)
ann.partition <- function(x, model = ann.partition.model.1 ) {
  # Wrapper function for partition() using the fixed threshold method
  partition(x, FUN = ".ann", model = model)
}


##########################
# Non-exported functions #

.fixed <- function(x, t = 3) {
  if(t < 0) stop("t must be a positive number")
  t.upper = rep(t, ncol(x))
  t.lower = rep(t * -1, ncol(x))
  return(rbind(t.upper, t.lower))
}

.ann <- function(x, model = ann.partition.model.1) {
  result = apply(x, 2, FUN = function(y, m) {
    skew = moments::skewness(y)
    kurt = moments::kurtosis(y)
    pred = neuralnet::compute(m, matrix(c(skew, kurt), 1, 2))
    t.lower = pred$net.result[1]
    t.upper = pred$net.result[2]
    return(rbind(t.upper, t.lower))
  }, m = model)
  return(result)
}

.checkMatrix <- function(x) {
  # This function checks to make sure x is a properly formed numeric matrix or, if not,
  # coerces it into one if possible.
  # It is not exported.
  
  d = dim(x)
  if(is.null(d) | (length(d) != 2)) {
    stop("x must be a 2D matrix")
  }
  check.for.factors = sapply(x, is.factor)
  if(any(check.for.factors)){
    stop("x must not contain factor data")
  }
  
  # Preserve names
  rn = rownames(x)
  cn = colnames(x)
  
  # Convert to numeric (in case of character data)
  x = apply(x, 2, as.numeric)
  dim(x) = d
  
  if(sum(is.na(x)) > 0) {
    stop("x must not contain missing values or values coerced to NA")
  }
  
  rownames(x) = rn
  colnames(x) = cn
  return(x) 
}