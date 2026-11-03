complete_stats_table = function(stats_table, comp_result){
  betamat = comp_result$betamat
  betamat.lower = comp_result$betamat.lower
  betamat.upper = comp_result$betamat.upper
  pmat = comp_result$pmat
  familymat = comp_result$familymat
  statmat = comp_result$statmat
  
  effectsizearray = abind::abind(beta=betamat,beta.lower=betamat.lower,beta.upper=betamat.upper,along = 3)
  familyarray = abind::abind(familymat,familymat,familymat,along = 3)
  nonlinear = !is.na(familyarray) & familyarray != "gaussian"
  effectsizearray[nonlinear] = exp(effectsizearray[nonlinear])
  effectsizearray.round = round(effectsizearray,2)
  for(e in 1:2){
    sel = !is.na(effectsizearray) & abs(effectsizearray) > 10^e
    effectsizearray.round[sel] = round(effectsizearray[sel],2-e)
  }
  chareffectsizemat = matrix("",nrow(pmat),ncol(pmat))
  dimnames(chareffectsizemat) = dimnames(pmat)
  prefixes.1 = c(gaussian = "Δ = ", quasipoisson = "r  = ", binomial = "OR = ")
  for(k in 1:ncompare) chareffectsizemat[,k] = ifelse(is.na(betamat[,k]),"",
                                                      paste0(prefixes.1[familymat[,k]],
                                                             effectsizearray.round[,k,'beta'],
                                                             " [",effectsizearray.round[,k,'beta.lower']," ; ",
                                                             effectsizearray.round[,k,'beta.upper'],"]"))
  
  roundpmat = round(pmat,3)
  roundpmat[!is.na(pmat) & pmat >= 0.1] = round(pmat[!is.na(pmat) & pmat >= 0.1],2)
  exponent = floor(log(pmat)/log(10))
  digit = round(pmat/10^exponent)
  charpmat = matrix("",nrow(pmat),ncol(pmat))
  dimnames(charpmat) = dimnames(pmat)
  for(k in 1:ncompare) charpmat[,k] = ifelse(is.na(pmat[,k]),"",
                                             ifelse(pmat[,k] >= 0.95*10^(-3),
                                                    paste0(roundpmat[,k]),
                                                    paste0(digit[,k]," × 10<sup>",exponent[,k],"</sup>")))
  
  roundstatmat = round(statmat,2)
  roundstatmat[!is.na(statmat) & abs(statmat) > 10] = round(statmat[!is.na(statmat) & abs(statmat) > 10],1)
  charstatmat = charpmat
  prefixes.2 = c(gaussian = "t = ", quasipoisson = "t  = ", binomial = "χ<sup>2</sup> = ")
  for(k in 1:ncompare) charstatmat[,k] = ifelse(is.na(statmat[,k]),"",
                                                paste0(prefixes.2[familymat[,k]],
                                                       round(statmat[,k],2)))
  charstatpmat = charpmat
  for(k in 1:ncompare) charstatpmat[,k] = ifelse(is.na(pmat[,k]),"",
                                                 paste0(charstatmat[,k]," (",charpmat[,k],")"))
  #View(charstatpmat)
  
  #Benjamini-Hochberg
  pvect = as.vector(pmat)
  pvect = sort(pvect[!is.na(pvect)])
  ntests = length(pvect)
  pvect.signif = pvect[pvect < 0.05*(1:ntests)/ntests]
  p.threshold = max(c(0,pvect.signif))
  signifmat = !is.na(pmat) & pmat <= p.threshold
  
  for(k in 1:ncompare){
    chareffectsizemat[,k] = ifelse(signifmat[,k], paste0('<b>',chareffectsizemat[,k],'</b>'), chareffectsizemat[,k])
    charstatpmat[,k] = ifelse(signifmat[,k], paste0('<b>',charstatpmat[,k],'</b>'), charstatpmat[,k])
  }
  
  nr = nrow(stats_table)
  resmat.withempty = matrix("",nr,ncompare)
  colnames(resmat.withempty) = colnames(pmat)
  for(i in 1:nr){
    varname = stats_table[i,'Variable name']
    if(varname != ""){
      resmat.withempty[i,] = chareffectsizemat[varname,]
      resmat.withempty[i+1,] = charstatpmat[varname,]
    }
  }
  
  return(cbind(stats_table,resmat.withempty))
}
