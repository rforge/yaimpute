\name{grmsd}
\alias{grmsd}
\title{Generalized Root Mean Square Distance Between Observed and Imputed Values}

\description{
Computes the root mean square distance between predicted and corresponding
observed values in an orthoganal multivariate space. This statistic provides a way to
compare imputation (or prediction) results. While it is designed to work with
imputation, the function can be used objects that inherit from \code{\link{lm}}
or with matrices and data frames that follow the column naming convention
described in the details.
}

\usage{
grmsd(\dots,ancillaryData=NULL,vars=NULL,wts=NULL)
}

\arguments{
  \item{\dots}{objects created by any combination of
    \code{\link{yai}}, \code{\link{impute.yai}}, \code{\link{ensembleImpute}},
    \code{\link{buildConsensus}}, \code{\link{lm}} and data frames or matrices
    that follow the column naming convention described in details.}
  \item{ancillaryData}{a data frame that defines variables, passed to 
    \code{\link{impute.yai}}.}
  \item{vars}{a list of variable names you want to include; if NULL all available
    variables are included.}
  \item{wts}{A vector of weights used to compute the mean distances, see 
    details below.}
}

\details{ 
  This function is designed to compute the root mean square distance between observed
  and predicted observations over several variables at once. 
  It is similar to a Mahalanobis distance. Here are some notable characteristics.
  
  \enumerate{
    \item{In the univariate case this function returns the same value as 
      \code{\link{rmsd}} with \code{scale=TRUE}. In that case
      the root mean square difference is computed after \code{\link{scale}}
      has been called on the variable.}
    \item{Like \code{\link{rmsd}}, \code{grmsd} is zero if the imputed values are
      exactly the same as the observed values over all variables.}
    \item{Like \code{\link{rmsd}}, \code{grmsd} is 1.0 when the mean of each 
      variable is imputed in place of a near neighbor. This situation corresponds 
      to regression where the slope is zero. It indicates that the imputation 
      error is, over all, the same as it would be if the means of the varialbes 
      were imputed rather than a near neighbor (it does not signal that the means 
      were imputed).} 
    \item{Like \code{\link{rmsd}}, values of grmsd > 1.0 indicate that, on average,
      the errors in the imputation are greater than they would be if the mean
      of the corresponding variables were imputed for each observation.} 
    \item{Note that individual \code{\link{rmsd}} values can be computed even when 
      the variance of the variable is zero. In contrast, \code{grmsd} can 
      only be computed in the situation where the observed data matrix is full rank.}
  }
  
  Observed and predicted are matched using the column names. Column names
  that have "\code{.o}" are matched to those that do not. Columns that do not
  have matching observed and imputed (predicted) values are ignored. 
  
  Several objects may be passed as "\dots". Function \code{\link{impute.yai}} is
  called for any objects that were created by \code{\link{yai}};
  \code{ancillaryData} and \code{vars} are passed to \code{\link{impute.yai}}
  when it is used.
  
  When objects inherit from \code{\link{lm}}, a suitable matrix is formed using
  by calling the \code{\link{predict}} and \code{\link{resid}} functions.
  
  As noted above, the generalized distance can only be computed when the 
  observed data is full rank. Rank is determined using \code{\link{qr}} and 
  columns are removed from the analysis to create this condition if 
  necessary (with a warning). 
 
  Factors, if found, are removed (with a warning).
  
  When argument \code{wts} is defined there must be one value for each pair of
  observed and predicted variables. If the values are named (preferred), then
  the names are matched to the names of predicted variables (no \code{.o} suffix).
  The matched values are used to compute weighted average distances using
  function \code{\link{weighted.mean}}.
}

\value{
A sorted named vector of generalized root mean square distances where names are taken from the 
arguments.
}

\seealso{
\code{\link{yai}}, \code{\link{impute.yai}}, \code{\link{rmsd.yai}}, 
\code{\link{notablyDifferent}}
}


\examples{
require(yaImpute)

data(iris)
set.seed(12345)

# form some test data
refs=sample(rownames(iris),50)
x <- iris[,1:2]      # Sepal.Length Sepal.Width
y <- iris[refs,3:4]  # Petal.Length Petal.Width

# build yai objects using 2 methods
msn <- yai(x=x,y=y)
mal <- yai(x=x,y=y,method="mahalanobis")

# compute the average distances between observed and imputed (predicted)
grmsd(msn,mal,lmFit=lm(as.matrix(y) ~ ., data=x[refs,]))

# use the all variables and observations in iris
# Species is a factor and is automatically deleted with a warning
grmsd(msn,mal,ancillaryData=iris)

# here is an example using lm, and another using column
# means as predictions.

impMean <- y 
colnames(impMean) <- paste0(colnames(impMean),".o")
impMean <- cbind(impMean,y)
impMean[,3] <- mean(impMean[,3])
impMean[,4] <- mean(impMean[,4])

grmsd(msn, mal, lmFit=lm(as.matrix(y) ~ ., data=x[refs,]), impMean )

# compare to using function rmsd:
msnimp <- na.omit(impute(msn))
grmsd(msnimp[,c(1,3)])           # values match
rmsd(msnimp[,c(1,3)],scale=TRUE)

# these are multivariate cases and they don't match
# because the covariance of the two is > 0.
grmsd(msnimp)
colSums(rmsd(msnimp,scale=TRUE))/2

}


\author{
  Nicholas L. Crookston \email{ncrookston.fs@gmail.com}\cr
  John Coulston \email{jcoulston@fs.fed.us}
}

\keyword{multivariate}

