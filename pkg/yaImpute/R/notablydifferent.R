# borrows a lot of code from impute 

notablyDifferent <- function (object,vars=NULL,threshold=NULL,p=.05,...)
{
   if (missing(object)) stop ("object required.")
   if (class(object)[1] != "yai")  stop ("object must be class yai")
   refIds = rownames(object$neiDstRefs)
   if (length(refIds) == 0) stop ("references are required")
   trgIds = rownames(object$neiDstTrgs)
   impObj = impute.yai(object,vars=vars,observed=TRUE,...)
   if (is.null(impObj)) stop ("no imputations found using this object")
   nuke = unlist(lapply(impObj,function (x) all(is.na(x))))
   nuke=nuke[nuke]
   if (length(nuke) > 0) impObj = impObj[,-match(names(nuke),names(impObj)),drop=FALSE]
   nuke = unlist(lapply(impObj,function (x) is.factor(x)))
   nuke = nuke[nuke]
   if (length(nuke) > 0) impObj = impObj[,-match(names(nuke),names(impObj)),drop=FALSE]
   impObj = na.omit(impObj)
   if (is.null(vars)) vars = names(impObj)
   vi = paste(unique(strsplit(vars,".o",fixed=TRUE)))
   vi = intersect(vi,names(impObj))
   notFound = setdiff(vars,names(impObj))
   if (length(notFound)>0) warning ("variables not found or had missing values: ",
                                     paste(notFound,collapse=", "))
   if (length(vi) == 0) stop("nothing to compute")
   vo = paste(vi,"o",sep=".")
   notFound = setdiff(vo,names(impObj))
   if (length(notFound)>0) warning ("variables not found or had missing values: ",
                                     paste(notFound,collapse=", "))
   vo = intersect(vo,names(impObj))
   both = intersect(paste(unique(strsplit(vo,".o",fixed=TRUE))),vi)  
   if (length(both) == 0) stop("no variables with observed and imputed values")
   vo = paste(both,"o",sep=".")
   diff = (impObj[,both]-impObj[,vo])/unlist(lapply(impObj[,vo],sd))
   rmsd = sqrt(apply(diff * diff,1,mean))
   if (is.null(threshold)) threshold = quantile(rmsd[refIds],1-p)
   ans=list(notablyDifferent.refs=sort(rmsd[rmsd[refIds]>threshold]),
        notablyDifferent.trgs=sort(rmsd[rmsd[trgIds]>threshold]),
        threshold=threshold,vars=both,rmsdS.refs=sort(rmsd[refIds]),
        rmsdS.trgs=sort(rmsd[trgIds]))
   class(ans) = "notablyDifferent"
   ans
}




plot.notablyDifferent <- function (obj,add=FALSE,...)
{
  if (missing(obj)) stop ("object required.")
  if (class(obj)[1] != "notablyDifferent")  stop ("object must be class notablyDifferent")
  all = c(obj$rmsdS.refs,obj$rmsdS.trgs)
  pch = c(rep(1,length(obj$rmsdS.refs)),rep(2,length(obj$rmsdS.trgs)))
  names(pch) = names(all)
  pch[names(obj$notablyDifferent.refs)] = 18
  pch[names(obj$notablyDifferent.trgs)] = 20
  x = 1:length(all)
  if (add) points(x=x,y=all,pch=pch,...) else 
             plot(x=x,y=all,pch=pch,
                  ylab="Scaled RMSD",xlab="Observation",...)
  abline(h=obj$threshold,...)
}
  
  









