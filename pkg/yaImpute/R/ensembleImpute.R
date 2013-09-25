ensembleImpute <- function (imputes,method="mean",...)
{
  cl = match.call()
  posM = c("mean","median")
  if (!(method %in% posM)) stop ('method="',method,'" must be one of ',paste0('"',posM,'"',collapse=", "))

  for (i in 1:length(imputes)) if (class(imputes[[i]]) == "yai") 
    imputes[[i]] = impute.yai(imputes[[i]],...)
  colns = unique(unlist(lapply(imputes,function(x) colnames(x))))
  ctf = if (method!="mean") median else mean
  rowns = sort(unique(unlist(lapply(imputes,function (x) rownames(x)))))
  ave = list()
  sd = list()
  N = list()
  methods = list()
  for (cl in colns)
  {
    one = matrix(unlist(lapply(imputes,function (x,rowns,cl) 
          x[rowns,cl],rowns,cl)), nrow=length(rowns))
    n = apply(one,1,function (x) sum(!is.na(x)))
    if (any(n != length(imputes))) N[[cl]] = n
    if (mode(one) == "character") 
    {
      ave[[cl]] = apply(one,1,function(x) 
        {
          x = na.omit(x)
          if (length(x) == 0) return(NA)
          x = table(x)
          names(x)[which.max(x)]
        })
      ave[[cl]] = as.factor(ave[[cl]])
      methods[[cl]] = "mode"
    } else {
      ave[[cl]] = apply(one,1,ctf,na.rm=TRUE)
      ave[[cl]][is.nan(ave[[cl]])] = NA
      sd [[cl]] = apply(one,1,function (x) 
        {
          x = na.omit(x)
          if (length(x) > 1) sd(x) else 0
        })
      methods[[cl]] = method
    }
  }
  ave = as.data.frame(ave)
  rownames(ave) = rowns
  ans = list(ave=ave)
  if (length(sd)>0) 
  {
    sumsgtz = unlist(lapply(sd,sum)) > 0
    if (any(sumsgtz))
    {
      sd = as.data.frame(sd[sumsgtz])
      rownames(sd) = rowns
      ans$sd = sd
    }
  }
  ans$N = if (length(N)>0) as.data.frame(N) else length(imputes)
  ans$methods = unlist(methods)
  ans  
}  

