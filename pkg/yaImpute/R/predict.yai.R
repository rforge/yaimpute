predict.yai <- 
function(x,newdata,...)
{
  if (missing(newdata)) impute.yai(x,...) else 
  { 
    al <- list(...)
    impute.yai(newtargets(x,newdata,al$k,al$ann),...)
  }
}

