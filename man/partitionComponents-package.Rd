\name{partitionComponents-package}
\alias{partitionComponents-package}
\alias{partitionComponents}
\docType{package}
\title{
Partition component matrices
}
\description{
This package partitions component matrices, e.g. the S matrix
generated in independent component analysis.  The \code{\link{partition}}
function takes a numeric matrix as input and returns a copy of it with
new values that indicate the partitioning.  Extremely positive values in
the input matrix are replaced with 1, extremely negative values are
replaced with -1, and all other values are replaced with 0.
}
\details{
\tabular{ll}{
Package: \tab partitionComponents\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2016-09-20\cr
License: \tab GPL-2\cr
Depends: \tab neuralnet, moments
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Michael Cary

Maintainer: Michael Cary <michael.cary.dev@gmail.com>
}
\references{
Cary MP, Kenyon C. Inferring module-weighted gene annotations from expression data. \emph{Nature Communications.} 2016 (in revision).
}
\examples{
# Create a random matrix
x = matrix(rnorm(25), 5, 5)

# Partition using a fixed threshold
x.par = fixed.partition(x)

# Partition using a neural network model
x.par = ann.partition(x)
}
