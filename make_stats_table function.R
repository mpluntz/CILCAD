make_stats_table = function(data, varlist, type, filter.col=NULL, filter.row=NULL, roundcoeff=10, force.continuous=NULL, twolines = FALSE, inbrackets="N"){
  mynrows = 2*length(varlist) - 1 + (1+twolines)*sum(sapply(varlist,length))
  colnamesbase = c("Description","Variable name","Data displayed")
  if(is.null(filter.col)) mycolnames = c(colnamesbase, type$types)
  if(!is.null(filter.col)) mycolnames = c(colnamesbase, paste0(rep(type$types,2),rep(c("",filter.col$suffix),each=length(type$types))) )
  res = matrix("",mynrows,length(mycolnames))
  rownames(res) = rep("",mynrows)
  colnames(res) = mycolnames
  rowind = 0
  for(varcat in names(varlist)){
    rowind = rowind + 1
    res[rowind,1] = varcat
    rowind = rowind + 1
    for(varname in names(varlist[[varcat]])){
      datavect = data[[varname]]
      var.type = ifelse(is.logical(datavect) | all(datavect %in% c(0,1,NA)), "binary",
                        ifelse(all(datavect == round(datavect),na.rm=TRUE), "integer", "continuous"))
      if(!is.null(force.continuous)) if(varname %in% force.continuous) var.type = "continuous"
      res[rowind,1] = varlist[[varcat]][varname]
      res[rowind,2] = varname
      if(var.type == "binary"){
        res[rowind,3] = "n (%)"
      }else if(var.type == "integer"){
        res[rowind,3] = "median (quartiles)"
      }else res[rowind,3] = "mean ± st. dev."
      if(!twolines) res[rowind,3] = paste0(res[rowind,3]," [",inbrackets,"]")
      
      if(!is.null(filter.row) & varname %in% names(filter.row)){
        keepvect.row = data[[filter.row[[varname]]]] == TRUE
      } else keepvect.row  = rep(TRUE,length(datavect))
      
      for(typename in type$types){
        keepvect.cell = keepvect.row & data[[type$varname]] == typename
        filtvects = list(datavect[keepvect.cell])
        cnames = typename
        if(!is.null(filter.col)){
          filtvects[[2]] = datavect[keepvect.cell & data[[filter.col$varname]] == filter.col$value]
          cnames[2] = paste0(typename, filter.col$suffix)
        }
        for(j in 1:length(filtvects)) if(any(!is.na(filtvects[[j]]))){
          datavect.cell = filtvects[[j]]
          if(var.type == "binary"){
            num1 = sum(datavect.cell,na.rm = TRUE)
            num2 = round(1000*mean(datavect.cell,na.rm = TRUE))/10
            num2 = paste(num2,"%")
          }else if(var.type == "integer"){
            num1 = median(datavect.cell,na.rm = TRUE)
            q1 = quantile(datavect.cell,0.25,na.rm = TRUE)
            q3 = quantile(datavect.cell,0.75,na.rm = TRUE)
            num2 = paste0(q1," - ",q3)
          }else{
            num1 = round(roundcoeff*mean(datavect.cell,na.rm = TRUE))/roundcoeff
            num2 = round(roundcoeff*sd(datavect.cell,na.rm = TRUE))/roundcoeff
          }
          if(var.type == "continuous") resstr = paste0(num1," ± ",num2) else resstr = paste0(num1," (",num2,")")
          if(!twolines){
            if(inbrackets == "N") num3 = sum(!is.na(datavect.cell))
            if(inbrackets == "missing") num3 = sum(is.na(datavect.cell))
            resstr = paste0(resstr," [",num3,"]")
          }
          res[rowind,cnames[j]] = resstr
          if(twolines){
            num1 = sum(is.na(datavect.cell))
            num2 = round(100*mean(is.na(datavect.cell)))
            resstr = as.character(num1)
            if(num1 > 0){
              if(num2 > 0) resstr = paste0(resstr," (",num2," %)")
              if(num2 == 0) resstr = paste0(resstr," (< 1 %)")
            }
            resstr = paste0("[",resstr,"]")
            res[rowind+1,cnames[j]] = resstr
          }
        }
      }
      if(twolines) rowind = rowind + 2 else rowind = rowind + 1
    }
  }
  return(res)
}