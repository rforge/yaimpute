\name{ensembleImpute}
\alias{ensembleImpute}
\title{Computes the mean, median, or mode among a list of impute.yai objects}

\description{Several objects of class \code{\link{impute.yai}} or \code{\link{yai}} 
are combined by computing the mean, median, or mode of separate, individual imputations.
The intention is that the members of the first argument would be formed by running 
yai several times with \code{bootstrap=TRUE} or by varying other options. 
}

\usage{
ensembleImpute(imputes, method="mean",\dots)
}

\arguments{
  \item{imputes}{a list of objects class \code{\link{impute.yai}} or \code{\link{yai}}. 
        Function \code{\link{impute.yai}} is called for those that where class is yai.}
  \item{method}{when "mean", the continuous variables are averaged using \code{mean},
        otherwise the \code{median} is used. Mode is always used for character data 
        (generally the case for factors).}
  \item{\dots}{passed to \code{\link{impute.yai}}.}
}

\value{
  A named lis, as follows:
  \itemize{
    \item{\code{ave }}{the mean, median, or mode of the imputed values for 
      all observations (and variables) found in the input list.}
    \item{\code{sd }}{the corresponding standard deviations for continuous variables
       if there are any. The columns are not reported if the standard deviation is
       zero for all observations which is typically true of "observed" values. }
    \item{\code{N }}{the number of replictions used to compute the corresponding data; 
       reported only if the number differs from the total number of replications. This 
       will be the case when \code{bootstrap}, \code{sampleVar}, or both are used in 
       \code{\link{impute.yai}}}
  }   
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

reps <- replicate(10, yai(x=x,y=y,method="msn",bootstrap=TRUE,k=2),
                  simplify=FALSE)

ensembleImpute(reps,ancillaryData=iris)

}


\author{
  Nicholas L. Crookston \email{ncrookston.fs@gmail.com} \cr
  John Coulston \email{jcoulston@fs.fed.us}
}

\keyword{multivariate}

