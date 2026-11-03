load("data.UKB_LRB_ITA.Rdata")
data.UKB_LRB_ITA[,atrophy.CADAMRIT.hasBPF := ifelse(is.na(BPF), NA, atrophy.CADAMRIT)]
data.UKB_LRB_ITA[,atrophy.score.CADAMRIT.hasBPF := ifelse(is.na(BPF), NA, atrophy.score.CADAMRIT)]

type = list(varname = "population",
            types = c("ITA","LRB","UKB"))
varlist1 = list(Demography = c(age="Age", sex.female = "Female sex",
                            education = "Education level (years)"),
               Risk_factors = c(hypertension="Hypertension",
                            hypercholesterolemia = "Hypercholesterolemia",
                            cholesterol = "Cholesterol",
                            cholesterol.statin = "    In subjects taking statins",
                            cholesterol.nostatin = "    In subjects not taking statins",
                            LDL = "LDL Cholesterol",
                            LDL.statin = "    In subjects taking statins",
                            LDL.nostatin = "    In subjects not taking statins",
                            has.statin = "Statin use",
                            diabetes = "Diabetes",
                            tobacco= "Smoking",
                            alcohol = "Alcohol",
                            riskfactor.atleast1 = "At least 1 risk factor"),
               Outcomes = c(stroke = "At least 1 stroke",
                            MMSE = "MMSE",
                            Rankin_at_least_2 = "mRS > 1",
                            MMSE_under_24 = "MMSE ≤ 24",
                            dementia = "Dementia"),
               #Quantitative_MRI = c(BPF = "BPF (%)",
                            #WMHpc = "WMH volume (%)",
                            #nb.lacunes = "Lacunes (N)",
                            #nb.CMB = "Microbleeds (N)",
                            #has.lacune.quantitative = "at least one lacune",
                            #has.CMB.quantitative = "At least one microbleed"),
               CADAMRIT = c(WMHscore = "Total WMH score",
                            has.lacune.CADAMRIT = "At least one lacune",
                            has.CMB.CADAMRIT = "At least one microbleed",
                            atrophy.CADAMRIT = "Deep or superficial atrophy",
                            CSVD.burden.CADAMRIT = "CSVD burden score"))
                            #atrophy.CADAMRIT.hasBPF = "Atrophy with BPF data",
                            #atrophy.score.CADAMRIT = "Atrophy score",
                            #atrophy.score.CADAMRIT.hasBPF = "Atrophy score with BPF data"))


varlist = list(Main = c(varlist1$Demography ,varlist1$Risk_factors))

#varlist = list(Outcomes = varlist1$Outcomes, CADAMRIT = c(age.MRI = "Age at imaging", varlist1$CADAMRIT))

lapply(varlist, function(L) names(L) %in% colnames(data.UKB_LRB_ITA))

data.UKB_LRB_ITA$has.CADAMRIT = FALSE
for(varname in names(varlist1$CADAMRIT)) data.UKB_LRB_ITA$has.CADAMRIT = data.UKB_LRB_ITA$has.CADAMRIT | !is.na(data.UKB_LRB_ITA[[varname]])
data.UKB_LRB_ITA$no.statin = !data.UKB_LRB_ITA$has.statin

filter.row = list(cholesterol.statin = "has.statin", cholesterol.nostatin = "no.statin",
                  LDL.statin = "has.statin", LDL.nostatin = "no.statin")
for(varname in names(varlist$CADAMRIT)) filter.row[[varname]] = "has.CADAMRIT"


data = data.UKB_LRB_ITA[is.R1231C == 1]

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
                                varnames.noadjust = c("age","age.MRI","sex.female"),
                                adjust1 = c("MMSE","MMSE_under_24"), adjust1_on = "education")


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










