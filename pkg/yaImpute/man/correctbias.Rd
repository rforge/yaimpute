\name{correctBias}
\alias{correctBias}
\title{Correct bias by selecting different near neighbors}

\description{
   This is an \emph{EXPERIMENTAL FUNCTION} \cr Change the neighbor selections
   in a \code{\link{yai}} object such that overall bias in the average value of a \emph{variable} 
   is reduced to be within a defined confidence interval. The \emph{variable} can be an 
   expression of several variables.
}

\usage{
correctBias(object,trgVal,trgValCI=NULL,nStdev=1.5,trace=FALSE)
}

\arguments{
  \item{object}{an object of class \code{\link{yai}} with \code{k > 1}.}
  \item{trgVal}{an \code{\link{expression}} defining how the target
    is to be computed, or a character string that can be coerced into
    an expression.}
  \item{trgValCI}{The confidence interval of the desired target. If the target
    falls within this interval, the problem is solved. If \code{NULL}, the
    value is computed for the target using the \code{reference} observations, 
    \code{trgVal}, and \code{nStdev}.}
  \item{nStdev}{the number of standard deviations in the target used to 
    compute the confidence interval when one is needed, ignored if \code{trgValCI} 
    is not \code{NULL}.}
  \item{trace}{if \code{TRUE}, detailed output is produced.}
}

\details{
  First, note that the approach may change as this is an \emph{EXPERIMENTAL FUNCTION}.
  \cr\cr
  Neighbors are switched with other selections such that the bias is reduced. This 
  is done using an \emph{ad hoc} search approach. Multiple passes through the data are 
  attempted until the objective of getting the bias to be within the confidence 
  interval is met. In a given pass the following steps are followed. \cr\cr 
  First, the target value is computed among the target observations; this is accomplished 
  by creating a data frame of imputed values for all the target observations 
  and executing the \code{trgVal} expression in the \link{environment} of that data frame.
  The imputed values in this situation are based on setting \emph{k} equal to 
  the pass number (so, \code{k = 1} for the first pass).  \cr\cr
  The result of this calculation is a vector with length equal to the number 
  of target observations. The \code{trgVal} minus the mean of this vector 
  is the current bias. Individual differences between \code{trgVal} and 
  the vector members indicate the contribution to bias made by each target. \cr\cr  
  The target observations are then ordered into decreasing order of 
  distance, so that those targets that are further away from their 
  corresponding reference are first in the ordered list. Neighbors are 
  switched in this priority sequence until the bias falls to be within the 
  confidence interval. Note that the target may not be met. 
}

\value{
 An object of class \code{\link{yai}} where \code{k = 1} and the neighbor selections
 have been changed so that the \code{trgVal} 
 fails within the \code{trgValCI}. In addition, the \code{call} element is changed
 to show both the original call to \code{\link{yai}} and the call to this function. 
 A new list, called \code{biasParameters} is added with these tags:
 \item{trgValCI}{the target CI.}
 \item{curVal}{the value of the bias that was achieved.}
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

# find 5 refernece neighbors for each target
msn <- yai(x=x,y=y2,method="msn",k=5)

# check for and correct for bias. Neighbor selections will be
# changed to bring the imputed values into line with the CI.
# in this case, no changes are made (npasses returns as zero).
msnCorr = correctBias(msn,trgVal="Petal.Width")
msnCorr$biasParameters

}

\author{
  Nicholas L. Crookston \email{ncrookston.fs@gmail.com} \cr
}

\keyword{misc}
\keyword{multivariate}


