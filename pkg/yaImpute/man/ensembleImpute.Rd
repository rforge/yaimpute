\name{ensembleImpute}
\alias{ensembleImpute}
\title{Computes the mean or mode among a list of impute.yai objects}

\description{Several objects of class \code{\link{impute.yai}} are combined into a new
object by computing the mean or mode of the values. The intention is that the many would
be formed by running yai several times with \code{bootstrap=TRUE} or by varying
other options, calling impute on each to get imputations, and then calling this function.
}

\usage{
ensembleImpute(imputes, method="mean")
}

\arguments{
  \item{imputes}{a list of objects class \code{\link{impute.yai}}.}
  \item{method}{when "mean", the continuous variables are averaged using \code{mean},
        otherwise the mode is computed. Mode is always used for factors.}
}

\value{
  An object of class \code{\link{impute.yai}}
}

\seealso{
\code{\link{yai}} \code{\link{buildConsensus}} \code{\link{impute.yai}} 
}

\examples{

require (yaImpute)
data(iris)

set.seed(123) 

# form some test data, y's are defined only for reference
# observations.
refs=sample(rownames(iris),50)
x <- iris[,1:2]      # Sepal.Length Sepal.Width
y <- iris[refs,3:4]  # Petal.Length Petal.Width

reps <- replicate(20, yai(x=x,y=y,method="msn",bootstrap=TRUE,k=2),
                  simplify=FALSE)

imputes <- lapply(reps,impute.yai)

ensembleImpute(imputes)

}


\author{
  Nicholas L. Crookston \email{ncrookston.fs@gmail.com} \cr
  John Coulston \email{jcoulston@fs.fed.us} \cr
}

\keyword{multivariate}

