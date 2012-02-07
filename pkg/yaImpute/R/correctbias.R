correctBias = function (object,trgVal,trgValCI=NULL,nStdev=1.5,excludeRefIds=NULL,trace=FALSE)
{
  if (class(object) != "yai") stop("object must be of class yai")

  if (missing(trgVal)) stop("trgVal must be defined.")
  trgValExpression = if (is.expression(trgVal)) trgVal else
    if (class(trgVal) == "character") parse(text=trgVal) else NULL
  if (!is.expression(trgValExpression)) stop ("trgVal can not be coerced into an expression.")

  trgVals <- NULL
  trgVal <-  NULL
  trgValsd <- NULL
  trgValCIc <- NULL
  excludeRefIds <- if (is.null(excludeRefIds)) "none" else excludeRefIds
  if (excludeRefIds == "all")
  {
    refsToCount <- NULL      
  }
  else
  {
    refsToCount <- 1:nrow(object$yRefs)
    trgValData <- as.data.frame(cbind(object$yRefs,object$xRefs))
    trgValDataCopy <- trgValData
    if (class(excludeRefIds) == "character") remove <- 
      na.omit(match(excludeRefIds,rownames(trgValDataCopy)))
    if (length(remove) > 0) 
    {
      trgValDataCopy <- trgValDataCopy[-remove,]
      refsToCount <- refsToCount[-remove]
      if (length(refsToCount) == 0) refsToCount <- NULL 
    }
    trgVals <- eval(trgValExpression,trgValDataCopy)
    trgVal <- mean(trgVals)
    trgValsd <- sd(trgVals)/sqrt(nrow(trgValDataCopy))
    trgValCIc <- c(trgVal-(trgValsd*nStdev),trgVal+(trgValsd*nStdev))
  }
  
  # if the user specified the CI, use it instead.
  if (is.null(trgValCI)) trgValCI <- trgValCIc
  if (is.null(trgValCI)) stop ("trgValIC was not supplied and can not be computed given the data.")
    
  if (trace) cat ("Target CI=", trgValCI,"\n");flush.console()
  
  newobject <- object
  pass = 0
  if (object$k-1 > 0) 
  {
    for (pass in 1:(object$k-1))
    {
      curVals <- eval(trgValExpression,trgValData[newobject$neiIdsTrgs[,1],,drop=FALSE])
      curVal <- if (is.null(trgVals)) mean(curVals) else mean(c(curVals,trgVals))
      if (curVal >= trgValCI[1] & curVal <= trgValCI[2]) 
      {
        if (trace) cat ("trgValCI=",trgValCI," curVal=",curVal)
        if (pass == 1)
        {
          pass = 0
          if (trace) cat (" -- no bias to correct.\n"); flush.console()
        }
        else
        {
          pass=pass-1
          if (trace) cat (" -- target CI reached, passes made=",pass,"\n");flush.console()
        }
        break
      }
    
      toHigh <- curVal > trgValCI[2]
      curBias <- if (toHigh) curVal - trgValCI[2] else curVal - trgValCI[1]
      # dV is the contribution to the bias of each observation.
      dV <- (eval(trgValExpression,trgValData[newobject$neiIdsTrgs[,pass+1],,drop=FALSE]) - curVals) / length(curVals)
      dDst <- newobject$neiDstTrgs[,pass+1] - newobject$neiDstTrgs[,1]
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
      
      newobject$neiDstTrgs[ord[1:toSwitch],1] <- newobject$neiDstTrgs[ord[1:toSwitch],pass+1]
      newobject$neiIdsTrgs[ord[1:toSwitch],1] <- newobject$neiIdsTrgs[ord[1:toSwitch],pass+1]
    }
  }
  else
  {
    pass = 0
    curVal = trgVal
  }
 
  # get rid of the k>1 columns...
  newobject$neiDstTrgs <- newobject$neiDstTrgs[,1,drop=FALSE] 
  newobject$neiIdsTrgs <- newobject$neiIdsTrgs[,1,drop=FALSE] 
  newobject$neiDstRefs <- newobject$neiDstRefs[,1,drop=FALSE] 
  newobject$neiIdsRefs <- newobject$neiIdsRefs[,1,drop=FALSE] 
  newobject$k <- 1
  newobject$call <- list(yai=newobject$call,correctBias=match.call())
  newobject$biasParameters <- list(trgValCI = trgValCI, curVal = curVal, npasses = pass)
  newobject
}


