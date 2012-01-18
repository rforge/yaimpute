\name{correctBias}
\alias{correctBias}
\title{Find notably distant targets}

\description{
   \emph{This is an EXPERIMENTAL FUNCTION} \cr This function changes the neighbor selections
   in a \code{\link{yai}} object such that overall bias in a \emph{variable} or
   in a combination defined by an \emph{expression} are reduced to be a user-defined 
   confidence interval.
}

\usage{
correctBias(object,trgValVariable="<some variable>")
}

\arguments{
  \item{object}{an object of class \code{\link{yai}} with \code{k > 1}.}
  \item{trgValExpression}{an \code{\link{expression}} defining how the target
    is to be computed, if \code{NULL}, the \code{trgValVariable} must be defined.}
  \item{trgValVariable}{a variable defining how the target 
    is to be computed, if \code{NULL}, the \code{trgValExpression} must be defined.} 
  \item{trgValCI}{The confidence interval of the desired target. If the target
    falls within this interval, the problem is solved. If \code{NULL}, the
    value is computed for the target using the \code{reference} observations and
    \code{nStdev}.}
  \item{nStdev}{the number of standard deviations in the target used to 
    compute the confidence interval, ignored if \code{trgValCI} is not \code{NULL}.}
  \item{trace}{if \code{TRUE}, detailed output is produced.}
}

\details{
  How this is done is being worked on.
}

\value{
 An object of class \code{\link{yai}} where \code{k = 1} and the neighbor selections
 have been changed so that the \code{trgValExpression} or \code{trgValVariable} 
 fails within the \code{trgValCI}. In addition, the \code{call} element is changed
 to show both the original call to \code{\link{yai}} and the call to this function. 
 A new list, called \code{biasParameters} is added with these tags:
 \item{trgValCI}{the call.}
 \item{curVal}{the value that was achieved.}
 \item{npasses}{the number of passes through the data taken to achieve the result.}
}

\seealso{
\code{\link{yai}}}

\examples{
data(iris)

# form some test data
refs=sample(rownames(iris),50)
x <- iris[,1:3]      # Sepal.Length Sepal.Width Petal.Length
y <- iris[refs,4:5]  # Petal.Width Species

# build an msn run, first build dummy variables for species.

sp1 <- as.integer(iris$Species=="setosa")
sp2 <- as.integer(iris$Species=="versicolor")
y2 <- data.frame(cbind(iris[,4],sp1,sp2),row.names=rownames(iris))
y2 <- y2[refs,]

names(y2) <- c("Petal.Width","Sp1","Sp2")

msn <- yai(x=x,y=y2,method="msn",k=5)

msnCorr = correctBias(msn,trgValVariable="Petal.Width")
msnCorr$biasParameters

}

\author{
  Nicholas L. Crookston \email{ncrookston.fs@gmail.com} \cr
}

\keyword{misc}
\keyword{multivariate}


