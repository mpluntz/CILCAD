data.ITA = readxl::read_excel("LRB and ITA data/2026_09_Cilento_data.xlsx",sheet="CILCAD_stroke_data")
data.ITA = as.data.table(data.ITA)
data.ITA[,age_stroke := as.numeric(age_stroke)]
data.ITA[,stroke := as.numeric(stroke)]

data.ITA[,age_info := age_stroke]
data.ITA[,age_stroke := ifelse(stroke == 1, age_info, NA)]


### NAIVE FORMULA ###
age = 80
table(data.ITA[age_info >= age | age_stroke < age, (stroke & age_stroke < age)])
mean(data.ITA[age_info >= age | age_stroke < age, (stroke & age_stroke < age)])

age = min(data.ITA[,age_info],na.rm=TRUE):max(data.ITA[,age_info],na.rm=TRUE)
penetrance = numeric(length(age))
for(k in seq_along(age)) penetrance[k] = mean(data.ITA[age_info >= age[k] | age_stroke < age[k],
                                                       (stroke & age_stroke < age[k])])
plot(age,1 - penetrance,type='l',las = 1, ylim = c(0,1))
### NAIVE FORMULA END ###



km_fit = survival::survfit(Surv(age_info, stroke) ~ 1, data=data.ITA)

plot(km_fit,xlab='Age',ylab='Probability of absence of stroke',las=1,col='darkblue')
legend("bottomleft",c("Kaplan-Meier estimator","95% CI bounds"),lty=c(1,2),col='darkblue')

km_fit_bysex = survival::survfit(Surv(age_info, stroke) ~ SEX, data=data.ITA)
plot(km_fit_bysex,xlab='Age',ylab='Probability of absence of stroke',las=1,col='darkblue')
legend("bottomleft",c("Kaplan-Meier estimator","95% CI bounds"),lty=c(1,2),col='darkblue')

km_fit_dt = as.data.table(summary(km_fit,data.frame=TRUE))
max(km_fit_dt[time < 60,1-surv])
max(km_fit_dt[time < 80,1-surv])


