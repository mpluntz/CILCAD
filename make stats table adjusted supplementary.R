load("data.UKB_LRB_ITA.Rdata")

type = list(varname = "population",
            types = c("ITA","LRB","UKB"))
varlist1 = list(Outcomes = c(stroke = "At least 1 stroke",
                             MMSE = "MMSE",
                             Rankin_at_least_2 = "mRS > 1",
                             MMSE_under_24 = "MMSE ≤ 24",
                             dementia = "Dementia"),
                CADAMRIT = c(WMHscore = "Total WMH score",
                             has.lacune.CADAMRIT = "At least one lacune",
                             has.CMB.CADAMRIT = "At least one microbleed",
                             atrophy.CADAMRIT = "Atrophy",
                             CSVD.burden.CADAMRIT = "CSVD burden score"))

varlist = list(Outcomes = varlist1$Outcomes,CADAMRIT = c(varlist1$CADAMRIT))

lapply(varlist, function(L) names(L) %in% colnames(data.UKB_LRB_ITA))

data.UKB_LRB_ITA$has.CADAMRIT = FALSE
for(varname in names(varlist1$CADAMRIT)) data.UKB_LRB_ITA$has.CADAMRIT = data.UKB_LRB_ITA$has.CADAMRIT | !is.na(data.UKB_LRB_ITA[[varname]])
data.UKB_LRB_ITA$no.statin = !data.UKB_LRB_ITA$has.statin

filter.row = list()
for(varname in names(varlist$CADAMRIT)) filter.row[[varname]] = "has.CADAMRIT"

data = data.UKB_LRB_ITA[is.R1231C == 1 & !is.na(hypertension) & !is.na(hypercholesterolemia)]

table(data[,population])
table(data[has.CADAMRIT == 1,population])


source("stats tables/make_stats_table function.R")
stats_table.R1231C = make_stats_table(data,varlist, type,
                                      roundcoeff=10, force.continuous = c("age","age.MRI"), twolines = TRUE,
                                      filter.row = filter.row)

source("stats tables/make_comparisions function.R")

ncompare = 3
varnames = unlist(lapply(varlist,names))
pops1 = c("ITA","ITA","LRB")[1:ncompare]
pops2 = c("LRB","UKB","UKB")[1:ncompare]
comp_result = make_comparisions(data, varnames,pops1,pops2,
                                varnames.MRI = names(varlist1$CADAMRIT),
                                varnames.noadjust = character(0),
                                adjust1 = c("MMSE","MMSE_under_24"), adjust1_on = "education",
                                adjust2 = varnames, adjust2_on = c("hypertension","hypercholesterolemia"))

source("stats tables/complete_stats_table function.R")

complete_stats_table.R1231C = complete_stats_table(stats_table.R1231C, comp_result)

nr = nrow(complete_stats_table.R1231C)
#keeprow = rep(TRUE,nr)
#for(i in seq_len(nr)){
#  if(complete_stats_table.R1231C[i,'Description'] == 'Main') keeprow[i] = FALSE
#  if(all(complete_stats_table.R1231C[i,] == '')) keeprow[i] = FALSE
#}
complete_stats_table.R1231C = complete_stats_table.R1231C[2:nr,]

complete_stats_table.R1231C[,'Description'] = gsub("CADAMRIT","CADA-MRIT",complete_stats_table.R1231C[,'Description'],fixed=TRUE)
colnames(complete_stats_table.R1231C) = gsub("_"," - ",colnames(complete_stats_table.R1231C),fixed=TRUE)
colnames(complete_stats_table.R1231C) = gsub("ITA","CILCAD",colnames(complete_stats_table.R1231C),fixed=TRUE)
colnames(complete_stats_table.R1231C) = gsub("Description","",colnames(complete_stats_table.R1231C),fixed=TRUE)

nc = ncol(complete_stats_table.R1231C)
kableres = knitr::kable(complete_stats_table.R1231C[,c(1,4:nc)], format = "html", escape = FALSE)
kableExtra::kable_styling(kableres)