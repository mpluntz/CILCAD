load("data.UKB_LRB_ITA.Rdata")

data.UKB_LRB_ITA[,noLDL := 1*is.na(LDL)]
data.UKB_LRB_ITA[,subgroup := paste(population,ifelse(noLDL,"wo","w"),"LDL data")]

varlist = list(characteristics = c(age="Age",sex.female="Female sex",
                                   stroke = "Stroke",has.statin = "Statin use",
                                   LDL = "LDL level"))

type = list(varname = "subgroup",
            types = c("LRB w LDL data","LRB wo LDL data",
                      "UKB w LDL data","UKB wo LDL data"))

source("stats tables/make_stats_table function.R")
stats_table.noLDL = make_stats_table(data.UKB_LRB_ITA, varlist, type,
                                     roundcoeff=10, force.continuous = "age", twolines = FALSE,
                                     filter.col = list(varname = 'is.R1231C', value = 1, suffix = ' R1231C'))


## univariate regressions ##

data.UKB_LRB_ITA[,missing.LDL := is.na(LDL)]
data.UKB_LRB_ITA[,missing.lacunes := is.na(has.lacune.CADAMRIT)]


popgroups = c("CILCAD","LRB R1231C","UKB R1231C","LRB","UKB")
statingroups = c("Overall","No statin")
covars.LDL = c("sex.female","age","stroke","has.lacune.CADAMRIT","has.statin","tobacco","hypertension")

pmat.LDL = matrix(NA,2*length(popgroups),length(covars.LDL))
popstatingroups = paste(rep(popgroups,each=2),rep(statingroups,length(popgroups)))
rownames(pmat.LDL) = popstatingroups
colnames(pmat.LDL) = covars.LDL

nmissing.LDL = rep(NA,2*length(popgroups))
names(nmissing.LDL) = popstatingroups
pcmissing.LDL = nmissing.LDL

for(popgroup in popgroups) for(statingroup in statingroups){
  popstatingroup = paste(popgroup,statingroup)
  print(popstatingroup)
  print(paste(popgroup,statingroup))
  if(popgroup == "CILCAD") pop = "ITA" else pop = substr(popgroup,1,3)
  data = data.UKB_LRB_ITA[population == pop]
  if(grepl("R1231C",popgroup)) data = data[is.R1231C == 1]
  if(statingroup == "No statin") data = data[has.statin == 0]
  nmissing.LDL[popstatingroup] = sum(data[,missing.LDL])
  pcmissing.LDL[popstatingroup] = 100*round(mean(data[,missing.LDL]),4)
  if(sum(data[,missing.LDL]) > 1) for(covarname in covars.LDL) if(covarname != "has.statin" | statingroup == "Overall"){
    formula = paste("missing.LDL","~",covarname)
    glmres = logistf::logistf(formula,data=data,family="binomial")
    pmat.LDL[popstatingroup,covarname] = glmres$prob[covarname]
  }
}

min(pmat.LDL,na.rm=TRUE)
min(pmat.lacunes,na.rm=TRUE)

makeresmat = function(nmissing,pcmissing,pmat){
  nmissing.text = paste0(nmissing," (",round(pcmissing,1)," %)")
  resmat = cbind(n_missing = nmissing.text,round(pmat,3))
  return(resmat)
}

resmat.LDL = makeresmat(nmissing.LDL,pcmissing.LDL,pmat.LDL)
View(resmat.LDL)
