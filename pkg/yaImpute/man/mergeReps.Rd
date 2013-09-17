\name{mergeReps}
\alias{mergeReps}
\title{Merges a list of k-NN objects into one ensemble (EXPERIMENTAL)}

\description{Several objects of class \code{\link{yai}} are combined into a new
object forming a consensus among the several. 
}

\usage{
mergeReps(reps, noTrgs=FALSE, noRefs=FALSE, k=NULL)
}

\arguments{
  \item{reps}{a list of objects class \code{\link{yai}}.}
  \item{noTrgs}{If \code{TRUE} neighbor relationships for target observations are not merged.}
  \item{noRefs}{If \code{TRUE} neighbor relationships for reference observations are not merged.}
  \item{k}{If not specified, the minimum value of \code{k} among the objects is used.}
}

\value{
  An object of class \code{\link{yai}}
}

\seealso{
\code{\link{yai}} 
}

\examples{

require (yaImpute)
data(iris)

# form some test data, y's are defined only for reference
# observations.
refs=sample(rownames(iris),50)
x <- iris[,1:2]      # Sepal.Length Sepal.Width
y <- iris[refs,3:4]  # Petal.Length Petal.Width

set.seed(123); 

reps = list()
for (rep in 1:20)
{
  reps[[rep]] <- yai(x=x,y=y,method="msn",bootstrap=TRUE,k=2)
}

mergeReps(reps)

}


\author{
  Nicholas L. Crookston \email{ncrookston.fs@gmail.com} \cr
  Acknowledgement: This function was inspired by conversations with John Coulston.
}

\keyword{multivariate}

