\name{TrioSGL}
\alias{TrioSGL}

\title{Trio Model with a Combination of Lasso and Group Lasso Regularization}
\description{Fit a trio model via penalized maximum likelihood.
The model is fit for a path of values of the penalty parameter.
This package is based on Noah Simon, et al. (2011) <doi:10.1080/10618600.2012.681250>.}
\usage{
TrioSGL(X, index = NULL, maxit = 10000, thresh = 0.0001, min.frac = 0.01, 
    nlam = 25, lambdas = NULL, alpha = 0.95, gamma = 0.8, step = 1, reset = 20, 
		standardize = FALSE, verbose = FALSE)
}

\arguments{
  \item{X}{$X$ should be an input matrix of dimension n-obs by p-vars. The number of rows must be a multiple of 4 (case followed by 3 pseudo-controls).}
  \item{index}{A p-vector indicating group membership of each covariate}
  \item{maxit}{Maximum number of iterations to convergence}
  \item{thresh}{Convergence threshold for change in beta}
  \item{min.frac}{The minimum value of the penalty parameter, as a fraction of the maximum value}
  \item{nlam}{Number of lambda to use in the regularization path}
  \item{lambdas}{A user specified sequence of lambda values for fitting. We recommend leaving this NULL and letting TrioSGL self-select values}
  \item{alpha}{The mixing parameter. \code{alpha} = 1 is the lasso penalty. \code{alpha} = 0 is the group lasso penalty.}
  \item{gamma}{Fitting parameter used for tuning backtracking (between 0 and 1)}
  \item{step}{Fitting parameter used for inital backtracking step size (between 0 and 1)}
  \item{reset}{Fitting parameter used for taking advantage of local strong convexity in nesterov momentum (number of iterations before momentum term is reset)}
  \item{standardize}{Logical flag for variable standardization prior to fitting the model.}
  \item{verbose}{Logical flag for whether or not step number will be output}
}
\details{
The sequence of models along the regularization path is fit by accelerated generalized gradient descent.
}
\value{An object with S3 class \code{"TrioSGL"}
\item{beta}{A p by \code{nlam} matrix, giving the penalized MLEs for the nlam different models, where the index corresponds to the penalty parameter \code{lambda}}
\item{lambdas}{The actual sequence of \code{lambda} values used (penalty parameter)}
\item{X.transform}{A list used in \code{predict} which gives the empirical mean and variance of the x matrix used to build the model}
}
\references{
 Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011). \emph{A Sparse-Group Lasso}. Journal of Computational and Graphical Statistics, 22, 231-245.
}
\author{Timo Stöcker, Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani\cr
Maintainer: Timo Stöcker \email{timo.stoecker@hhu.de}
}

\examples{
trios <- 4 * 10; snps <- 20; size.groups <- 4
index <- ceiling(1:snps / size.groups)
x <- floor(matrix(runif(trios * snps, min = 0, max = 3), ncol = snps, nrow = trios))
fit <- TrioSGL(x, index)
}
\keyword{trio}
\keyword{model}
\keyword{regression}
