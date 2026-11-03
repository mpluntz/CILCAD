load("data.UKB_LRB_ITA.Rdata")

data = data.UKB_LRB_ITA
data.imputed = data

varnames.toimpute = c("LDL","hypertension","tobacco")
varnames.imputefrom = c("sex.female","stroke","population","has.statin","is.R1231C")


#V1: imputation from the closest matching age even if the age different is large
for(k in 1:nrow(data)){
  is.na.row = is.na(t(data[k,])[,1])
  varnames.toimpute.local = varnames.toimpute[is.na.row[varnames.toimpute]]
  if(length(varnames.toimpute.local) > 0){
    #print(c(data$population[k],data$id[k]))
    #print(varnames.toimpute.local)
    data.1 = data
    for(varname in varnames.imputefrom) data.1 = data.1[data.1[[varname]] == data[[varname]][k]]
    agetarget = data$age[k]
    for(varname.toimpute in varnames.toimpute.local){
      data.2 = data.1[!is.na(data.1[[varname.toimpute]])]
      minagedist = min(abs(data.2$age - agetarget))
      data.2 = data.2[abs(age - agetarget) == minagedist]
      #print(c(minagedist,nrow(data.2)))
      if(minagedist > 3){
        print(c(data$population[k],data$id[k]))
        print(varnames.toimpute.local)
      }
      data.imputed[[varname.toimpute]][k] = mean(data.2[[varname.toimpute]])
    }
  }
}


#
idmissing = data[population == "ITA" & is.na(hypertension) & stroke == 1,id]
View(data.imputed[id %in% idmissing,])


#V2: imputation only from close ages, if not possible "conservative" imputation of risk factors: equal to stroke
#one LRB R1231C has LDL imputed from non-R1231C
agediff.limit = 3
varnames.toimpute = c("LDL","hypertension","tobacco")
varnames.imputefrom = c("sex.female","stroke","population","has.statin")
data.imputed = data


for(k in 1:nrow(data)){
  is.na.row = is.na(t(data[k,])[,1])
  varnames.toimpute.local = varnames.toimpute[is.na.row[varnames.toimpute]]
  if(length(varnames.toimpute.local) > 0){
    #print(c(data$population[k],data$id[k]))
    #print(varnames.toimpute.local)
    data.1 = data
    for(varname in varnames.imputefrom) data.1 = data.1[data.1[[varname]] == data[[varname]][k]]
    data.1 = data.1[abs(age - data$age[k]) <= agediff.limit]
    for(varname.toimpute in varnames.toimpute.local){
      data.2 = data.1[!is.na(data.1[[varname.toimpute]])]
      if(nrow(data.2) > 0){
        data.imputed[[varname.toimpute]][k] = mean(data.2[[varname.toimpute]])
      }else if(varname.toimpute %in% c("hypertension","tobacco")){
        print('***')
        for(varname in varnames.imputefrom) print(paste(varname,data[[varname]][k]))
        data.imputed[[varname.toimpute]][k] = data$stroke[k]
      }
    }
  }
}
table(is.na(data.imputed$LDL))



#V2: imputation only from close ages, if not possible "conservative" imputation of risk factors: equal to stroke
#one LRB R1231C has LDL imputed from non-R1231C
agediff.limit = 3
varnames.toimpute = c("LDL","hypertension","tobacco","has.lacune.CADAMRIT")
varnames.imputefrom = c("sex.female","stroke","population","has.statin")
data.imputed = data


for(k in 1:nrow(data)){
  is.na.row = is.na(t(data[k,])[,1])
  varnames.toimpute.local = varnames.toimpute[is.na.row[varnames.toimpute]]
  if(length(varnames.toimpute.local) > 0){
    #print(c(data$population[k],data$id[k]))
    #print(varnames.toimpute.local)
    data.1 = data
    for(varname in varnames.imputefrom) data.1 = data.1[data.1[[varname]] == data[[varname]][k]]
    data.1 = data.1[abs(age - data$age[k]) <= agediff.limit]
    for(varname.toimpute in varnames.toimpute.local){
      data.2 = data.1[!is.na(data.1[[varname.toimpute]])]
      if(nrow(data.2) > 0){
        data.imputed[[varname.toimpute]][k] = mean(data.2[[varname.toimpute]])
      }else if(varname.toimpute %in% c("hypertension","tobacco")){
        print('***')
        for(varname in varnames.imputefrom) print(paste(varname,data[[varname]][k]))
        data.imputed[[varname.toimpute]][k] = data$stroke[k]
      }
    }
  }
}
table(is.na(data.imputed$LDL))




#V3: basic imputation by mean


data = data.UKB_LRB_ITA
data.imputed = data

varnames.toimpute = c("hypertension","tobacco")
varnames.imputefrom = c("sex.female","stroke","population")


#V1: imputation from the closest matching age even if the age different is large
for(k in 1:nrow(data)){
  is.na.row = is.na(t(data[k,])[,1])
  varnames.toimpute.local = varnames.toimpute[is.na.row[varnames.toimpute]]
  if(length(varnames.toimpute.local) > 0){
    #print(c(data$population[k],data$id[k]))
    #print(varnames.toimpute.local)
    data.1 = data
    for(varname in varnames.imputefrom) data.1 = data.1[data.1[[varname]] == data[[varname]][k]]
    agetarget = data$age[k]
    for(varname.toimpute in varnames.toimpute.local){
      data.2 = data.1[!is.na(data.1[[varname.toimpute]])]
      #print(c(minagedist,nrow(data.2)))
      data.imputed[[varname.toimpute]][k] = mean(data.2[[varname.toimpute]])
    }
  }
}
