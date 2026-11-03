load("data.UKB_LRB_ITA.Rdata")


display.type = "basic"

#add.riskfactors = c("tobacco","hypertension")
add.riskfactors = character(0)


popgroups = c("LRB all mutations","UKB all mutations")
statingroups = c("Overall","No statin")

par(mfrow=c(2,2))
par(mar=c(5.1, 4.1, 3.1, 2.1))
letterindex = 1
for(popgroup in popgroups) for(statingroup in statingroups){
  if(popgroup == "CILCAD") pop = "ITA" else pop = substr(popgroup,1,3)
  data = data.UKB_LRB_ITA[population == pop & !is.na(LDL) & !is.na(has.lacune.CADAMRIT)]
  if(!grepl("all mutations",popgroup)) data = data[is.R1231C == 1]
  if(statingroup == "No statin") data = data[has.statin == 0]
  ntot = nrow(data)
  ncase = sum(data[,has.lacune.CADAMRIT == 1])
  ncontrol = sum(data[,has.lacune.CADAMRIT == 0])
  boxnames = paste0(c("No","Yes")," (N = ",c(ncontrol,ncase),")")
  main = paste0(statingroup," (N = ",ntot,")")
  main = paste(pop,main)
  
  ### test ###
  formula = "has.lacune.CADAMRIT ~ LDL + age.MRI + sex.male"
  
  for(addvarname in add.riskfactors){
    formula = paste0(formula," + ",addvarname)
    data = data[!is.na(data[[addvarname]])]
  }
  if(statingroup == "Overall") formula = paste0(formula," + has.statin")
  
  glmres = logistf::logistf(formula,data=data,family="binomial",control=logistf::logistf.control(maxit=100))
  pval = glmres$prob['LDL']
  
  #### plot ###
  boxplot(LDL ~ has.lacune.CADAMRIT, data = data, outline = FALSE,
          col = c("lightgreen","lightblue"), ylim = ylim, las = 1,
          xlab = "At least one lacune", ylab = "LDL level", names = boxnames, main = main, cex.main = 1)
  stripchart(LDL ~ has.lacune.CADAMRIT, data = data, vertical = TRUE, method = "jitter",
             pch = 20, col = c("darkgreen","darkblue"), add = TRUE)
  bary = 250
  points(c(1.1,1.1,1.9,1.9), bary+c(-10,0,0,-10), type='l')
  if(pval < 0.001) pvaltext = "p < 0.001" else pvaltext = paste0('p = ',round(pval,3))
  text(1.5,bary+10, labels=pvaltext,cex=1)
  text(0.1, 315, labels=LETTERS[letterindex],xpd=NA,cex=2)
  letterindex = letterindex + 1
}



# stats tables #

popgroups = c("CILCAD","LRB Arg1231Cys","UKB Arg1231Cys","LRB all mutations","UKB all mutations")
statingroups = c("Overall","No statin")
add.riskfactors.list = list(noadj=character(0),adj=c("tobacco","hypertension"))
resarrays = list()

for(addrf in c("noadj","adj")){
  add.riskfactors = add.riskfactors.list[[addrf]]
  resdimnames = list(popgroups,statingroups,c('N','n','beta','ci.lower','ci.upper','p','chisq'))
  resarray = array(NA,sapply(resdimnames,length))
  dimnames(resarray) = resdimnames
  
  for(popgroup in popgroups) for(statingroup in statingroups){
    if(popgroup == "CILCAD") pop = "ITA" else pop = substr(popgroup,1,3)
    data = data.UKB_LRB_ITA[population == pop & !is.na(LDL) & !is.na(has.lacune.CADAMRIT)]
    if(!grepl("all mutations",popgroup)) data = data[is.R1231C == 1]
    if(statingroup == "No statin") data = data[has.statin == 0]
    ntot = nrow(data)
    ncase = sum(data[,has.lacune.CADAMRIT == 1])
    ncontrol = sum(data[,has.lacune.CADAMRIT == 0])
    boxnames = paste0(c("No","Yes")," (N = ",c(ncontrol,ncase),")")
    main = paste0(statingroup," (N = ",ntot,")")
    main = paste(pop,main)
    
    ### test ###
    formula = "has.lacune.CADAMRIT ~ LDL + age.MRI + sex.male"
    
    for(addvarname in add.riskfactors){
      formula = paste0(formula," + ",addvarname)
      data = data[!is.na(data[[addvarname]])]
    }
    if(statingroup == "Overall") formula = paste0(formula," + has.statin")
    
    glmres = logistf::logistf(formula,data=data,family="binomial",control=logistf::logistf.control(maxit=100))
    pval = glmres$prob['LDL']
    resarray[popgroup,statingroup,] = c(nrow(data),sum(data$has.lacune.CADAMRIT),
                                        glmres$coefficients['LDL'],
                                        glmres$ci.lower['LDL'], glmres$ci.upper['LDL'],
                                        pval, qchisq(1-pval,df=1))
  }
  resarrays[[addrf]] = resarray
}


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


resarray = resarrays$adj

OR20 = exp(20*resarray[,,c('beta','ci.lower','ci.upper')])
OR20 = round(OR20,3)

OR20char = matrix(paste0(OR20[,,'beta'],' [',OR20[,,'ci.lower'],' ; ',OR20[,,'ci.upper'],']'),
                  dim(OR20)[1:2],dimnames=dimnames(OR20)[1:2])


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
colnames(ncol) = 'N with at least 1 lacune'

rescharmat.adj = cbind(Ncol,ncol,OR20charcol,pcol,chisqcol)
View(rescharmat.adj)


#BH
pmat.noadjust = resarrays$noadj[,,'p']
pmat.adjust = resarrays$adj[,,'p']


#pmat = pmat.noadjust
#pmat = pmat.adjust

pvect = as.vector(pmat)
ntests = length(pvect)
qvect = sort(pvect)*ntests/(1:ntests)

for(lag in 1:(ntests-1)) qvect[1:(ntests-lag)] = pmin(qvect[1:(ntests-lag)],qvect[(1+lag):ntests])

cbind(p=sort(pvect),q=qvect)