load("data.UKB_LRB_ITA.Rdata")

#global sur les variables de la régression des boxplots: ajustée sur âge sexe mais pas facteurs de risque
data = data.UKB_LRB_ITA[,c("is.R1231C","age","sex.male","LDL","population","has.statin","stroke")]
for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
naniar::mcar_test(data)

#global seulement les R1231C
data = data.UKB_LRB_ITA[is.R1231C == 1, c("age","sex.male","LDL","population","has.statin","stroke")]
naniar::mcar_test(data)


for(pop in c("ITA","LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[population == pop & is.R1231C == 1,c("age","sex.male","LDL","has.statin","stroke")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}
for(pop in c("LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[population == pop,c("is.R1231C","age","sex.male","LDL","has.statin","stroke")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}



data = data.UKB_LRB_ITA[,c("is.R1231C","age","sex.male","LDL","population","has.statin","stroke","tobacco","hypertension")]
for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
naniar::mcar_test(data)

data = data.UKB_LRB_ITA[,is.R1231C == 1, c("age","sex.male","LDL","population","has.statin","stroke","tobacco","hypertension")]
for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
naniar::mcar_test(data)



for(pop in c("ITA","LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[population == pop & is.R1231C == 1,c("age","sex.male","LDL","has.statin","stroke","tobacco","hypertension")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}
for(pop in c("LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[population == pop,c("is.R1231C","age","sex.male","LDL","has.statin","stroke","tobacco","hypertension")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}


data = data.UKB_LRB_ITA[,c("is.R1231C","age","sex.male","LDL","has.statin","stroke","tobacco","hypertension")]
print(naniar::mcar_test(data))

for(varname in colnames(data)[colnames(data) != "hypertension"]) print(paste(varname,mean(data[[varname]][!is.na(data$hypertension)],na.rm=TRUE)))



data = data.UKB_LRB_ITA[!is.na(age.MRI), c("is.R1231C","age.MRI","sex.male","LDL","population","has.statin","has.lacune.CADAMRIT")]
for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
naniar::mcar_test(data)


data = data.UKB_LRB_ITA[!is.na(age.MRI) & is.R1231C == 1, c("age.MRI","sex.male","LDL","population","has.statin","has.lacune.CADAMRIT")]
for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
naniar::mcar_test(data)

for(pop in c("ITA","LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[!is.na(age.MRI) & population == pop & is.R1231C == 1,c("age","sex.male","LDL","has.statin","has.lacune.CADAMRIT")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}
for(pop in c("LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[!is.na(age.MRI) & population == pop,c("is.R1231C","age","sex.male","LDL","has.statin","has.lacune.CADAMRIT")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}

for(pop in c("ITA","LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[!is.na(age.MRI) & population == pop & is.R1231C == 1,c("age","sex.male","LDL","has.statin",
                                                                                 "has.lacune.CADAMRIT","tobacco","hypertension")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}
for(pop in c("LRB","UKB")){
  print(pop)
  data = data.UKB_LRB_ITA[!is.na(age.MRI) & population == pop,c("is.R1231C","age","sex.male","LDL","has.statin",
                                                                "has.lacune.CADAMRIT","tobacco","hypertension")]
  for(varname in colnames(data)) print(paste(varname,sum(is.na(data[[varname]]))))
  print(naniar::mcar_test(data))
}



pdimnames = list(pop = c("all_pops","ITA","LRB","UKB"),
                 mut = c("R1231C","all_mutations"),
                 statin = c("overall","no_statin"),
                 rf = c("without_rf","with_rf"),
                 outcome = c("stroke","lacunes"))
parray = array(NA,sapply(pdimnames,length))
dimnames(parray) = pdimnames

for(pop in pdimnames$pop) for(mut in pdimnames$mut) for(statin in pdimnames$statin)  for(rf in pdimnames$rf)  for(outcome in pdimnames$outcome){
  print(paste(pop,mut,statin,rf,outcome))
  
  data = data.UKB_LRB_ITA
  if(pop != "all_pops") data = data[population == pop]
  if(mut == "R1231C") data = data[is.R1231C == 1]
  if(statin == "no_statin") data = data[has.statin == 0]
  if(outcome == "lacunes") data = data[!is.na(age.MRI)]
  
  varnames = c("sex.male","LDL")
  if(pop == "all_pops") varnames = c(varnames,"population")
  if(mut == "overall" & pop != "ITA") varnames = c(varnames,"is.R1231C")
  if(statin == "overall") varnames = c(varnames,"has.statin")
  if(rf == "with_rf") varnames = c(varnames,"tobacco","hypertension")
  if(outcome == "stroke") varnames = c(varnames,"stroke","age")
  if(outcome == "lacunes") varnames = c(varnames,"has.lacune.CADAMRIT","age.MRI")
  data = data[,..varnames]
  
  if(any(is.na(data))){
    for(varname in colnames(data)) if(any(is.na(data[[varname]]))) print(paste(varname,sum(is.na(data[[varname]]))))
    parray[pop,mut,statin,rf,outcome] = naniar::mcar_test(data)$p.value
  }
}


  