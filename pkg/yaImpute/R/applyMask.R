applyMask <- function
(object, refGroups=NULL, trgGroups=NULL, method="removeWhenCommon", k=1)
{
  if (class(object) != "yai") stop("object must be of class yai")
  valid=c("removeWhenCommon","keepWhenCommon")
  if (is.na(match(method,valid))) stop (paste("method must be one of",paste(valid,collapse=", ")))
  if (is.null(refGroups) | is.null(trgGroups)) stop("refGroups and trgGroups must be defined")
  if (k >= object$k) stop("new value of k (",k,") must be less than old value (",object$k,")")

  refGrp <- refGroups[match(object$neiIdsTrgs,rownames(object$xRefs))]
  lrefGrp <- if (method=="removeWhenCommon") refGrp != trgGroups else refGrp == trgGroups
  dim(lrefGrp) <- dim(object$neiIdsTrgs) 

  # tvec is an offset in the storage of neiIdsTrgs and neiDstTrgs. 
  tvec=0:(ncol(lrefGrp)-1) * nrow(lrefGrp)
  ans <- apply(lrefGrp,1,function(x,tvec,k) tvec[x][1:k],tvec,k)
  dim (ans) = c(k,nrow(lrefGrp))
 
  tIds = vector(mode = "character", length = nrow(lrefGrp)*k)
  tDst = vector(mode = "numeric",   length = nrow(lrefGrp)*k)
  dim (tIds) = c(nrow(lrefGrp),k)
  dim (tDst) = c(nrow(lrefGrp),k)
  rownames(tIds) = rownames(object$neiIdsTrgs) 
  rownames(tDst) = rownames(object$neiDstTrgs) 
  colnames(tIds) = colnames(object$neiIdsTrgs)[1:k]
  colnames(tDst) = colnames(object$neiDstTrgs)[1:k]

  idx = 1:nrow(lrefGrp)
  for (i in 1:k)
  {
    jdx = idx+ans[i,]
    tIds[,i] = object$neiIdsTrgs[jdx]
    tDst[,i] = object$neiDstTrgs[jdx]
  }
 
  object$call = match.call()
  object$neiIdsTrgs = tIds
  object$neiDstTrgs = tDst
  object$k = k
  object$noRefs = TRUE
  object$neiIdsRefs = NULL
  object$neiDstRefs = NULL
  object
}

