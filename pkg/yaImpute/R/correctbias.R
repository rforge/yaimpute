correctBias = function (obj,trgValExpression=NULL,trgValVariable=NULL,
                           trgValCI=NULL,nStdev=1.5,trace=FALSE)
{
  if (class(obj) != "yai") stop("obj must be of class yai")
  if (obj$k < 2) stop("k in obj must be greater than 1")
  if (is.null(obj$neiIdsTrgs)) stop("No targets, nothing to correct.")

  if (is.null(trgValExpression) && !is.null(trgValVariable)) trgValExpression <- as.expression(as.name(trgValVariable))
  if (is.null(trgValExpression)) stop ("trgValExpression and trgValVariable are both NULL")

  trgValData <- as.data.frame(cbind(obj$yRefs,obj$xRefs)) 
  trgVals <- eval(trgValExpression,trgValData)
  trgVal <- mean(trgVals)
  trgValsd <- sd(trgVals)/sqrt(nrow(trgValData))
  trgValCIc <- c(trgVal-(trgValsd*nStdev),trgVal+(trgValsd*nStdev))
  
  # if the user specified the CI, use it instead.
  if (is.null(trgValCI)) trgValCI <- trgValCIc
  
  if (trace) cat ("Target CI=", trgValCI,"\n");flush.console()
  
  newObj <- obj
  for (pass in 1:(obj$k-1))
  {
    curVals <- eval(trgValExpression,trgValData[newObj$neiIdsTrgs[,1],,drop=FALSE])
    curVal <- mean(curVals)
    if (curVal >= trgValCI[1] & curVal <= trgValCI[2]) 
    {
      if (pass == 1)
      {
        if (trace) cat ("trgValCI=",trgValCI," curVal=",curVal," -- no bias to correct.\n"); flush.console()
        pass = 0
        break
      }
      else
      {
        pass=pass-1
        if (trace) cat ("Target CI reached, pass=",pass,"\n");flush.console()
        break
      }
    }

    toHigh <- curVal > trgValCI[2]
    curBias <- if (toHigh) curVal - trgValCI[2] else curVal - trgValCI[1]
    # dV is the average contribution to the bias of each observation.
    dV <- (eval(trgValExpression,biasData[newObj$neiIdsTrgs[,pass+1],,drop=FALSE]) - curVals) / length(curVals)
    dDst <- newObj$neiDstTrgs[,pass+1] - newObj$neiDstTrgs[,1]
    if (toHigh) {
      ntoTry <- sum(dV < 0)
      dDst[dV > 0] <- .Machine$double.xmax 
    } else {
      ntoTry <- sum(dV > 0)
      dDst[dV < 0] <- .Machine$double.xmax
    }
    ord <- order(dDst)
    bsum <- cumsum(dV[ord[1:ntoTry]])
    toSwitch <- which(bsum+curBias > 0)[1]
    if (is.na(toSwitch)) toSwitch <- length(bsum)

    if (trace) cat ("trgValCI=",trgValCI," pass=",pass," curVal=",curVal," curBias=",curBias,
       " ntoTry=",ntoTry,"toSwitch=",toSwitch,"\n");flush.console()
    
    newObj$neiDstTrgs[ord[1:toSwitch],1] <- newObj$neiDstTrgs[ord[1:toSwitch],pass+1]
    newObj$neiIdsTrgs[ord[1:toSwitch],1] <- newObj$neiIdsTrgs[ord[1:toSwitch],pass+1]
  }
 
  # get rid of the k>1 columns...
  newObj$neiDstTrgs <- newObj$neiDstTrgs[,1,drop=FALSE] 
  newObj$neiIdsTrgs <- newObj$neiIdsTrgs[,1,drop=FALSE] 
  newObj$neiDstRefs <- newObj$neiDstRefs[,1,drop=FALSE] 
  newObj$neiIdsRefs <- newObj$neiIdsRefs[,1,drop=FALSE] 
  newObj$k <- 1
  newObj$call <- list(yai=newObj$call,correctBias=match.call())
  newObj$biasParameters <- list(trgValCI = trgValCI, curVal = curVal, npasses = pass)
  newObj
}


