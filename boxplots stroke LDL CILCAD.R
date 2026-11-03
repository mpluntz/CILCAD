load("data.UKB_LRB_ITA.Rdata")

#test.type = "scoretest"
test.type = "Firth"
display.type = "complete"
#add.riskfactors = c("tobacco","hypertension")
add.riskfactors = character(0)

popgroups = c("CILCAD","LRB p.Arg1231Cys","UKB p.Arg1231Cys","LRB all cysteine-altering variants","UKB all cysteine-altering variants")
statingroups = c("overall","no statins")

if(test.type == "Firth"){
  resdimnames = list(popgroups,statingroups,c('N','n','beta','ci.lower','ci.upper','p','chisq'))
  resarray = array(NA,sapply(resdimnames,length))
  dimnames(resarray) = resdimnames
}


#data.UKB_LRB_ITA = data.imputed

ylim = c(30,280)
par(mfrow=c(4,2))
#par(mar=c(5.1, 4.1, 3.1, 2.1))
par(mar=c(3.1, 4.1, 2.1, 2.1))
letterindex = 1
options(warn=1)
for(popgroup in popgroups) for(statingroup in statingroups){
  print(paste(popgroup,statingroup))
  if(popgroup == "CILCAD") pop = "ITA" else pop = substr(popgroup,1,3)
  data = data.UKB_LRB_ITA[population == pop & !is.na(LDL)]
  if(!grepl("cysteine-altering variants",popgroup)) data = data[is.R1231C == 1]
  if(statingroup == "no statins") data = data[has.statin == 0]
  ntot = nrow(data)
  ncase = sum(data[,stroke == 1])
  ncontrol = sum(data[,stroke == 0])
  boxnames = paste0(c("No stroke","Stroke")," (N = ",c(ncontrol,ncase),")")
  main = paste0(statingroup," (N = ",ntot,")")
  main = paste(popgroup,main)
  
  ### test ###
  if(test.type == "scoretest") formula = "stroke ~ age + sex.male" else formula = "stroke ~ LDL + age + sex.male"
  if(statingroup == "overall") formula = paste0(formula," + has.statin")
  
  for(addvarname in add.riskfactors){
    formula = paste0(formula," + ",addvarname)
    data = data[!is.na(data[[addvarname]])]
  }
  
  if(test.type == "linear Wald"){
    lmres = lm(formula,data=data)
    pval = summary(lmres)$coefficient['LDL','Pr(>|t|)']
  }
  if(test.type == "Wald"){
    glmres = glm(formula,data=data,family="binomial")
    pval = summary(glmres)$coefficients['LDL','Pr(>|z|)']
  }
  if(test.type == "score test"){
    glmres = glm(formula,data=data,family="binomial")
    zstat = statmod::glm.scoretest(glmres, data[,LDL])
    pval = 1 - pchisq(zstat^2,df=1)
  }
  if(test.type == "Firth"){
    glmres = logistf::logistf(formula,data=data,family="binomial",
                              control=logistf::logistf.control(maxit=200),plcontrol=logistf::logistpl.control(maxit=200))
    pval = glmres$prob['LDL']
    resarray[popgroup,statingroup,] = c(nrow(data),sum(data$stroke),
                                        glmres$coefficients['LDL'],
                                        glmres$ci.lower['LDL'], glmres$ci.upper['LDL'],
                                        pval, qchisq(1-pval,df=1))
  }
  
  #### plot ###
  if(popgroup != "LRB p.Arg1231Cys"){
    boxplot(LDL ~ stroke, data = data, outline = FALSE,
            col = c("lightgreen","lightblue"), ylim = ylim, las = 1,
            xlab = "Stroke", ylab = "LDL level", names = boxnames, main = main, cex.main = 1)
    stripchart(LDL ~ stroke, data = data, vertical = TRUE, method = "jitter",
               pch = 20, col = c("darkgreen","darkblue"), add = TRUE)
    pval.char = paste('p =',round(pval,4))
    if(display.type == "complete") pval.char = paste(test.type, pval.char)
    bary = 250
    if(pval < 0.001) pvaltext = "p < 0.001" else pvaltext = paste0('p = ',round(pval,3))
    points(c(1.1,1.1,1.9,1.9), bary+c(-10,0,0,-10), type='l')
    text(1.5,bary+10, labels=pvaltext,cex=1)
    text(0.1, 315, labels=LETTERS[letterindex],xpd=NA,cex=1.5)
    letterindex = letterindex + 1
  }
}

OR20 = exp(20*resarray[,,c('beta','ci.lower','ci.upper')])
OR20 = round(OR20,3)

OR20char = matrix(paste0(OR20[,,'beta'],' [',OR20[,,'ci.lower'],' ; ',OR20[,,'ci.upper'],']'),
                  dim(OR20)[1:2],dimnames=dimnames(OR20)[1:2])

make1col = function(M){
  a = nrow(M)
  b = ncol(M)
  res = matrix(NA,a*b,1)
  rownames(res) = rep('',a*b)
  for(k in 1:a){
    rownames(res)[(k-1)*b+(1:b)] = paste(rownames(M)[k],colnames(M))
    res[(k-1)*b+(1:b),1] = M[k,]
  }
  return(res)
}
OR20charcol = make1col(OR20char)
colnames(OR20charcol) = 'OR for +20 mg/dL LDL [95 % CI]'

pcol = make1col(resarray[,,'p'])
colnames(pcol) = 'p-value'
pcol = round(pcol,3)

chisqcol = make1col(resarray[,,'chisq'])
colnames(chisqcol) = 'Chi-square statistic'
chisqcol = round(chisqcol,2)

Ncol = make1col(resarray[,,'N'])
colnames(Ncol) = 'N'
ncol = make1col(resarray[,,'n'])
colnames(ncol) = 'N strokes'

rescharmat = cbind(Ncol,ncol,OR20charcol,pcol,chisqcol)


pmat = resarray[,,'p']

#pmat = pmat.noadjust
#pmat = pmat.adjust

pvect = as.vector(pmat)
ntests = length(pvect)
qvect = sort(pvect)*ntests/(1:ntests)

for(lag in 1:(ntests-1)) qvect[1:(ntests-lag)] = pmin(qvect[1:(ntests-lag)],qvect[(1+lag):ntests])

cbind(p=sort(pvect),q=qvect)
