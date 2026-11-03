library(data.table)

load("CILCAD-like UKB.Rdata")

data.UKB.all = data.UKB[,.(id = f.eid,
                           population = "UKB",
                           age = age.inclusion,
                           age.MRI = age.MRI,
                           sex.male,
                           education = education.years,
                           hypertension.v1 = hypertension.atIncl, hypertension = hypertension.atIncl,
                           hypercholesterolemia = hypercholesterolemia.atIncl,
                           cholesterol = Cholesterol,
                           cholesterol.statin = ifelse(has.statin.atIncl, Cholesterol, NA),
                           cholesterol.nostatin = ifelse(has.statin.atIncl, NA, Cholesterol),
                           LDL = LDL,
                           LDL.statin = ifelse(has.statin.atIncl, LDL, NA),
                           LDL.nostatin = ifelse(has.statin.atIncl, NA, LDL),
                           has.statin = 1*has.statin.atIncl,
                           has.statin.atMRI = 1*has.statin.atMRI,
                           diabetes = diabetes.atIncl,
                           tobacco = tobacco.everoccasionally.atIncl,
                           alcohol = alcohol.daily.atIncl,
                           riskfactor.atleast1 = hypertension.atIncl | hypercholesterolemia.atIncl |
                             diabetes.atIncl | tobacco.everoccasionally.atIncl | alcohol.daily.atIncl,
                           stroke = stroke.atIncl,
                           MMSE = NA, Rankin_at_least_2 = NA, cognitive.impairment = NA,
                           dementia = allcausedementia.algorithmic.atIncl,
                           BPF = FSL.BPF,
                           WMHpc = FSL.WMH,
                           nb.lacunes, has.lacune.quantitative = ifelse(nb.lacunes == 0,0,1),
                           nb.CMB, has.CMB.quantitative = ifelse(nb.CMB == 0,0,1),
                           WMHscore = Periventricular_WMH + Deep_WMH + Superficial_WMH,
                           has.lacune.CADAMRIT = ifelse(nb.lacunes == 0,0,1),
                           has.CMB.CADAMRIT = ifelse(nb.CMB == 0,0,1),
                           atrophy.CADAMRIT = 1*(Superficial_Atrophy > 0 | Deep_Atrophy > 0),
                           atrophy.score.CADAMRIT = Superficial_Atrophy + Deep_Atrophy,
                           CSVD.burden.CADAMRIT = 1*(dPVS_CSO > 0 | dPVS_BG > 0) +
                             1*(Periventricular_WMH > 0 | Deep_WMH > 0) +
                             ifelse(nb.lacunes == 0,0,1) +
                             ifelse(nb.CMB == 0,0,1),
                           mut.type = type,
                           is.R1231C = ifelse(type == "R1231C", TRUE, FALSE),
                           is.7_34 = ifelse(type %in% c("R1231C","C1222G","R1143C","EGFR 7-34"), TRUE, FALSE))]

data.UKB.notch3 = data.UKB.all[mut.type!="no.notch3"]

for(var in colnames(data.UKB.notch3)) if(is.logical(data.UKB.notch3[[var]])){
  data.UKB.notch3[[var]] = ifelse(data.UKB.notch3[[var]] == TRUE,1,0)
}
data.UKB.notch3[is.na(BPF) & is.na(WMHpc) & is.na(nb.lacunes),age.MRI := NA]


load("LRB and ITA data/2025_12_09_Result_with_dementia.Rdata")
data.LRB_ITA = as.data.table(result)
rm(result)
data.LRB.1 = data.LRB_ITA[italy == 0]
rm(data.LRB_ITA)
colnames(data.LRB.1) = gsub(" ","_",colnames(data.LRB.1))
colnames(data.LRB.1) = gsub("<","under_",colnames(data.LRB.1))
colnames(data.LRB.1) = gsub("≥","at_least_",colnames(data.LRB.1))
colnames(data.LRB.1) = gsub("(%)","pct",colnames(data.LRB.1),fixed=TRUE)
colnames(data.LRB.1) = gsub("-",".",colnames(data.LRB.1))

data.LRB.MRI = readxl::read_excel("LRB and ITA data/2025_01_26_LRB_MRI.xlsx")
data.LRB.MRI = as.data.table(data.LRB.MRI)
data.LRB.MRI.2 = data.LRB.MRI[TimePoint == "M0",.(id = as.character(Patient_Number),
                               age.MRI = Age,
                               timepoint = TimePoint,
                               nb.CMB.T2star = Nb_microbleeds_T2Star_final,
                               nb.CMB.SWI = Nb_microbleeds_SWI_final,
                               WMH_normTIV.2 = 100*as.numeric(WMH_final)/as.numeric(ANTS_wholeMaskTIV))]
data.LRB.2 = merge(data.LRB.1, data.LRB.MRI.2,all.x = TRUE)


data.LRB = data.LRB.2[,.(id,
                         population = "LRB",
                         age = Age, #this is age_bioch2021
                         age.MRI = Age,
                         sex.male = ifelse(Sex==1,1,0),
                         education = as.numeric(education_level),
                         hypertension.v1 = as.numeric(as.character(real_hypertension)),
                         hypertension = pmax(as.numeric(as.character(real_hypertension)),
                                             1*(as.numeric(pression_art_dias) >= 90 | as.numeric(pression_art_systo) >= 140)),
                         hypercholesterolemia = as.numeric(as.character(real_hypercholesterolemia)),
                         cholesterol = as.numeric(total_chol),
                         cholesterol.statin = ifelse(statine == 1, as.numeric(total_chol), NA),
                         cholesterol.nostatin = ifelse(statine == 1, NA, as.numeric(total_chol)),
                         LDL = as.numeric(chol_ldl),
                         LDL.statin = ifelse(statine == 1, as.numeric(chol_ldl), NA),
                         LDL.nostatin = ifelse(statine == 1, NA, as.numeric(chol_ldl)),
                         has.statin = as.numeric(as.character(statine)), has.statin.atMRI = as.numeric(as.character(statine)),
                         diabetes = as.numeric(as.character(real_diabetes)),
                         tobacco = as.numeric(Smoking),
                         alcohol = as.numeric(Alcool),
                         riskfactor.atleast1 = as.numeric(as.character(real_hypertension)) |
                           as.numeric(as.character(real_hypercholesterolemia)) |
                           as.numeric(as.character(real_diabetes)) | as.numeric(Smoking) | as.numeric(Alcool),
                         stroke = as.numeric(as.character(At_least_one_stroke)),
                         MMSE = as.numeric(indexmms),
                         Rankin_at_least_2 = as.numeric(as.character(Rankin_at_least_2)),
                         cognitive.impairment = NA,
                         dementia = Dementia,
                         BPF = BPF_pct, WMHpc = WMH_normTIV,
                         nb.lacunes = nb_lacunes, has.lacune.quantitative = ifelse(nb_lacunes == 0,0,1),
                         nb.CMB = as.numeric(nb.CMB.SWI),
                         has.CMB.quantitative = ifelse(nb.CMB.SWI == 0,0,1),
                         WMHscore = as.numeric(PVH) + as.numeric(DeepWMH) + as.numeric(SuperficialWMH),
                         has.lacune.CADAMRIT = ifelse(Lacune == 0,0,1),
                         has.CMB.CADAMRIT = ifelse(CMB == 0,0,1),
                         atrophy.CADAMRIT = 1*(Superficial_atrophy != '0' | Deep_atrophy != '0'),
                         atrophy.score.CADAMRIT = as.numeric(as.character(Superficial_atrophy)) +
                           as.numeric(as.character(Deep_atrophy)),
                         CSVD.burden.CADAMRIT = 1*(PVS.CSO != 0 | PVS.BG != 0) +
                           1*(PVH != 0 | DeepWMH != 0) +
                           ifelse(Lacune == 0,0,1) +
                           ifelse(CMB == 0,0,1),
                         mut.type = substr(as.character(grp_dom_egf),5,nchar(as.character(grp_dom_egf))))]

data.LRB[substr(mut.type,1,4) == "EGFr", mut.type := paste("EGFR",substr(mut.type,5,nchar(mut.type)))]
data.LRB[,is.R1231C := ifelse(mut.type == "R1231C", TRUE, FALSE)]
data.LRB[,is.7_34 := ifelse(mut.type %in% c("R1231C","EGFR 7-34"), TRUE, FALSE)]

data.UKB_LRB = rbind(data.UKB.notch3, data.LRB)

save(data.UKB_LRB,file="data.UKB_LRB.Rdata")


### separately obtain ITA data ###

data.ITA = readxl::read_excel("LRB and ITA data/2026_09_Cilento_data.xlsx",sheet="Cilento data (53) NEW")
data.ITA = as.data.table(data.ITA)
data.ITA[,ID:=as.character(ID)]

for(varname in colnames(data.ITA)){
  data.ITA[[varname]][data.ITA[[varname]] == "NA"] = NA
  data.ITA[[varname]] = as.numeric(data.ITA[[varname]])
}

colnames(data.ITA) = gsub(" ","_",colnames(data.ITA),fixed=TRUE)
colnames(data.ITA) = gsub("-","_",colnames(data.ITA),fixed=TRUE)
colnames(data.ITA) = gsub("(","",colnames(data.ITA),fixed=TRUE)
colnames(data.ITA) = gsub(")","",colnames(data.ITA),fixed=TRUE)

data.ITA.2 = data.ITA[,.(id = ID, population = "ITA",
                         age = age, age.MRI = age_Imaging,
                         sex.male = sex,
                         education = Education_Y,
                         hypertension.v1 = Hypertension, hypertension = Hypertension,
                         hypercholesterolemia = Hypercholesterolemia,
                         cholesterol = TC,
                         cholesterol.statin = ifelse(Statins == 1, TC, NA),
                         cholesterol.nostatin = ifelse(Statins == 1, NA, TC),
                         LDL = LDL,
                         LDL.statin = ifelse(Statins == 1, LDL, NA),
                         LDL.nostatin = ifelse(Statins == 1, NA, LDL),
                         has.statin = Statins, has.statin.atMRI = Statins,
                         diabetes = Diabetes, tobacco = Smoking, alcohol = Alchool,
                         riskfactor.atleast1 = 1*(Hypertension | Hypercholesterolemia |
                                                    Diabetes | Smoking | Alchool),
                         stroke = Stroke, MMSE = MMSE, Rankin_at_least_2 = 1*(RANKIN >= 2),
                         cognitive.impairment = Cognitive_impairment, dementia = Dementia,
                         BPF = BPF, WMHpc = normWMH,
                         nb.lacunes = Lacunes_MRI, has.lacune.quantitative = At_least_one_lacune_MRI,
                         nb.CMB = Number_of_microbleeds_MRI, has.CMB.quantitative = At_least_one_microbleed_MRI,
                         WMHscore = WMH_score,
                         has.lacune.CADAMRIT = At_least_one_lacune_CADA_MRIT,
                         has.CMB.CADAMRIT = At_least_one_microbleed_CADA_MRIT,
                         atrophy.CADAMRIT = Atrophy_CADA_MRIT,
                         atrophy.score.CADAMRIT = Atrophy_score,
                         CSVD.burden.CADAMRIT = 1*(Dilated_perivascular_spaces_dPVS_Centrum_semiovale_CSO != 0 | Dilated_perivascular_spaces_dPVS_Basal_ganglia != 0) +
                           1*(Peri_WMH != 0 | Deep_WMH != 0) +
                           At_least_one_lacune_CADA_MRIT +
                           At_least_one_microbleed_CADA_MRIT,
                         mut.type = "R1231C", is.R1231C = TRUE, is.7_34 = TRUE)]


data.UKB_LRB_ITA = rbind(data.UKB_LRB, data.ITA.2)

data.UKB_LRB_ITA[,MMSE_under_24 := 1*(MMSE <= 24)]
data.UKB_LRB_ITA[,MMSE_under_24_or_ci := ifelse(population=="ITA",cognitive.impairment, MMSE_under_24)]

data.UKB_LRB_ITA[,sex.female := 1 - sex.male]


data.UKB_LRB_ITA[,sex.female := 1 - sex.male]

save(data.UKB_LRB_ITA,file="data.UKB_LRB_ITA.Rdata")