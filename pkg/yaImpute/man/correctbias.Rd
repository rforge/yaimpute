\name{correctBias}
\alias{correctBias}
\title{Correct bias by selecting different near neighbors}

\description{
   This is an \emph{EXPERIMENTAL FUNCTION}
   \cr 
   Change the neighbor selections in a \code{\link{yai}} object such that 
   bias (if any) in the average value of a \emph{expression} 
   of one or more \emph{variables} is reduced to be within a defined confidence interval. 
}

\usage{
correctBias(object,trgVal,trgValCI=NULL,nStdev=1.5,trace=FALSE)
}

\arguments{
  \item{object}{an object of class \code{\link{yai}} with \code{k > 1}.}
  \item{trgVal}{an \code{\link{expression}} defining an attribute or combination of
    attributes is to be computed, or a character string that can be coerced into
    an expression. The expression can refer to one or more X- and Y-variables 
    defined for the reference observations.}
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
  First, note that the approach and discussion presented below may change as 
  this is an \emph{EXPERIMENTAL FUNCTION}.
  \cr\cr
  Imputation as it is defined in \code{yaImpute} can yield biased
  results. Lets say that you have a collection of reference
  observations that happen to be selected in a non-biased way among a
  population. In this discussion, \emph{population} is a finite set of all
  individual sample units of interest; the reference plus target observations
  often represent this population (but this need not be true, see below). In this
  context, the term population is the same as the \emph{sampling frame}.
  If the average of a measured attribute is computed from
  this random sample, it is an unbiased estimate of the true mean. 
  \cr\cr  
  Using \code{\link{yai}}, while setting \emph{k=1}, values for each of several
  attributes are imputed from a single reference observation to a target
  observation. Once the imputation is done over all the observations, an
  average of any one measured attribute can be computed. There is no guarantee 
  that this average will be within a pre-specified confidence interval of the 
  similarly computed average among the unbiased sample that comprises the reference data.
  \cr\cr
  Experience shows that despite any lack of guarantee, the results are accurate.
  This tends to hold true when the reference data contains
  samples that cover the variation of the targets, even when they are not a
  random sample, and even if some of the reference observations are from sample
  units that are outside the target population. 
  \cr\cr
  Because there is no guarantee, and because the reference observations might
  profitably come from sample units beyond the those in the population 
  (so as to insure all kinds of targets have a matching reference), it is necessary 
  to test the imputation results for bias. If bias is found, it would be helpful to 
  do something to correct it.
  \cr\cr
  The \code{correctBias()} function is designed to check for, and correct for,
  bias by selecting alternative nearby reference observations to be imputed
  to targets that contribute to the bias. The idea is that even if
  one reference is closest to a target, its attribute of interest might be greater (or less)
  than the mean. An alternative neighbor, one that may be almost as
  close, might reduce the overall bias if it were used instead. If this is found
  to be the case, \code{correctBias()} switches the nighbor selections. 
  It makes as many switches it needs to make until the mean among the targets falls 
  within the specified confidence interval. Note that, \emph{there is no guarantee that
  the goal will be met.}
  \cr\cr
  The details of the method are: 
  \cr\cr
  1. An attribute of interest is established by naming one in the call with
  argument \code{tarVal}. Note that this can be a simple variable name enclosed
  in quotations marks or it can be an \code{\link{expression}} of one or more variables. 
  If the former, it is converted into an expression that is executed in the 
  environment of the reference observations (both the X- and Y-variables). 
  A confidence interval is computed for this value under the assumption that the 
  reference observations are an unbiased sample of the target population. This may 
  not be the case. Regardless, a confidence interval is \emph{necessary} and it can
  alternatively be supplied using \code{trgValCI}.
  \cr\cr
  2. One of several possible passes through the data are taken to find neighbor switches
  that will result in the bias being corrected. A pass includes computing the
  attribute of interest using values imputed to all the targets. This computation
  results in a vector with one element for each target observation. From this vector, the
  contribution toward the bias made by each target can be computed. The target
  observations are then ordered into increasing order of distance, so that those
  targets that are nearest to their corresponding reference are first in the
  ordered list. Reference observations that contribute to the bias are picked.
  That is, if the bias is positive, then targets that are positive deviations
  from the mean are identified and switched with the next most similar neighbors. 
  The number of possible passes is equal to \emph{k-1} where \emph{k} is
  set in the original call to \code{\link{yai}}. Note that switches are made
  among targets only, and never among reference observations that may make up the
  population. Reference observations are left to represent themselves. 
  \cr\cr
  3. The details of the \emph{FUTURE} argument \code{allTargetIds}. When computing the mean of
  the attribute, \code{correctBias()} must know which observations represent the
  population. Normally, all the target observations would be in this set. If all
  reference observations came from the population they would also be included in
  the bias calculations even though they will not be switched with more distant 
  neighbors to correct for discovered bias. In deed, this is the \emph{default}
  situation, and it takes effect when \code{allTargetIds} is left NULL. Otherwise,
  \code{allTargetIds} is a vector of row names from the reference plus target
  observations that make up the population. Bias is computed among those
  observations. 
}

\value{
 An object of class \code{\link{yai}} where \code{k = 1} and the neighbor selections
 have been changed as described above. In addition, the \code{call} element is changed
 to show both the original call to \code{\link{yai}} and the call to this function. 
 A new list, called \code{biasParameters} is added to the \code{\link{yai}} object that 
 contains with these tags:
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


