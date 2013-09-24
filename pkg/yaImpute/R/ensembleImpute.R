ensembleImpute <- function (imputes,method="mean")
{
  cl = match.call()
  colns = unique(unlist(lapply(imputes,function(x) colnames(x))))

  ctf = if (method!="mean") median else mean
  rowns = sort(unique(unlist(lapply(imputes,function (x) rownames(x)))))
  ave = list()
  sd = list()
  for (cl in colns)
  {
    one = matrix(unlist(lapply(imputes,function (x,rowns,cl) 
          x[rowns,cl],rowns,cl)), nrow=length(rowns))
    if (mode(one) == "character") 
    {
      ave[[cl]] = apply(one,1,function(x) 
        {
          x = na.omit(x)
          if (length(x) == 0) return(NA)
          x = table(x)
          names(x)[which.max(x)]
        })
    } else {
      ave[[cl]] = apply(one,1,ctf,na.rm=TRUE)
      ave[[cl]][is.nan(ave[[cl]])] = NA
      sd [[cl]] = apply(one,1,sd,na.rm=TRUE)
    }
  }
  ave = as.data.frame(ave)
  rownames(ave) = rowns
  ans = list(ave=ave)
  if (length(sd)>0) 
  { 
    sd = as.data.frame(sd)
    rownames(sd) = rowns
    ans$sd = sd
  }
  ans  
}  

