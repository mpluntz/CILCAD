make_comparisions = function(data, varnames, pops1, pos2,
                             test.strategy = "Firth",nmaxit=100,
                             varnames.MRI = character(0),
                             varnames.noadjust = character(0),
                             varnames.quasipoisson = character(0),
                             adjust1 = character(0), adjust1_on = character(0),
                             adjust2 = character(0), adjust2_on = character(0)){
  ncompare = length(pops1)
  compare.names = paste(pops1,pops2,sep="_")
  pmat = matrix(NA,length(varnames),ncompare)
  row.names(pmat) = varnames
  colnames(pmat) = compare.names
  betamat = pmat
  betamat.upper = pmat
  betamat.lower = pmat
  statmat = pmat
  familymat = pmat
  formulamat = pmat
  
  for(varname in varnames) for(k in 1:ncompare){
    pop1 = pops1[k]
    pop2 = pops2[k]
    print(paste(varname,pop1,pop2))
    if(any(!is.na(data[population == pop1, get(varname)])) & any(!is.na(data[population == pop2, get(varname)]))){
      
      vars.used = c(varname, "population")
      if(!(varname %in% varnames.noadjust)){
        vars.used = c(vars.used,"sex.female")
        if(varname %in% varnames.MRI) vars.used = c(vars.used,"age.MRI") else vars.used = c(vars.used,"age")
        if(varname %in% adjust1) vars.used = c(vars.used, adjust1_on)
        if(varname %in% adjust2) vars.used = c(vars.used, adjust2_on)
      }
      data.used = data[population %in% c(pop1,pop2), ..vars.used]
      data.used[,population := ifelse(population == pop1,0,1)]
      
      for(covarname in colnames(data.used)) data.used = data.used[!is.na(data.used[[covarname]])]
      
      #type of model and test
      if(test.strategy == "Wald_scoretest") test.type = "waldtest" else test.type = "firthtest"
      if(all(data[[varname]] == 0 | data[[varname]] == 1,na.rm=TRUE)){
        family = "binomial"
        #test singularity on population
        is.singular = all(data.used[population == 0, get(varname)] == 0, na.rm=TRUE) |
          all(data.used[population == 1, get(varname)] == 0, na.rm=TRUE) |
          all(data.used[population == 0, get(varname)] == 1, na.rm=TRUE) |
          all(data.used[population == 1, get(varname)] == 1, na.rm=TRUE)
        if(is.singular & test.strategy == "Wald_scoretest"){
          test.type = "scoretest"
          print(paste("score test for",varname,pop1,pop2))
        }
      }else family = "gaussian"
      if(varname %in% varnames.quasipoisson) family = "quasipoisson"
      if(family != "binomial" & test.type == "firthtest") test.type = "waldtest"
      
      compname = compare.names[k]
      familymat[varname,compname] = family
      
      #variables in model
      if(test.type == "scoretest")  formula = paste(varname, "~ 1") else formula = paste(varname, "~ population")
      if(!(varname %in% varnames.noadjust)){
        formula = paste(formula, "+ sex.female")
        if(varname %in% varnames.MRI) formula = paste(formula, "+ age.MRI") else formula = paste(formula, "+ age")
        if(varname %in% adjust1) for(covar in adjust1_on) formula = paste(formula, "+", covar)
        if(varname %in% adjust2) for(covar in adjust2_on) formula = paste(formula, "+", covar)
      }
      formulamat[varname,compname] = formula
      
      if(test.type == "waldtest"){
        glmres = glm(formula,family=family,data=data.used)
        coeffs = summary(glmres)$coefficients
        beta.hat = coeffs['population',1]
        sd.hat = coeffs['population',2]
        betamat[varname,compname] = beta.hat
        betamat.lower[varname,compname] = beta.hat + qnorm(0.025)*sd.hat
        betamat.upper[varname,compname] = beta.hat + qnorm(0.975)*sd.hat
        statmat[varname,compname] = coeffs['population',3]
        pmat[varname,compname] = coeffs['population',4]
      }
      if(test.type == "firthtest"){
        glmres = logistf::logistf(formula,family=family,data=data.used,
                                  control=logistf::logistf.control(maxit=nmaxit),
                                  plcontrol=logistf::logistpl.control(maxit=nmaxit))
        betamat[varname,compname] = glmres$coefficients['population']
        betamat.lower[varname,compname] = glmres$ci.lower['population']
        betamat.upper[varname,compname] = glmres$ci.upper['population']
        statmat[varname,compname] = qchisq(1 - glmres$prob['population'], 1)
        pmat[varname,compname] = glmres$prob['population']
      }
      if(test.type == "scoretest"){
        glmres = glm(formula, family=family, data = data.used)
        zstat = statmod::glm.scoretest(glmres, data.used[,population])
        betamat[varname,compname] = NA
        betamat.lower[varname,compname] = NA
        betamat.upper[varname,compname] = NA
        statmat[varname,compname] = zstat
        pmat[varname,compname] = 1 - pchisq(zstat^2,df=1)
      }
    }
  }
  return(list(pmat = pmat,
              betamat = betamat,
              betamat.upper = betamat.upper,
              betamat.lower = betamat.lower,
              statmat = statmat,
              familymat = familymat,
              formulamat = formulamat))
}
