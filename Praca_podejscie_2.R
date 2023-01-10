library(survival)
library(survminer)

dane <- read.csv("clinical.tsv", sep = "\t")
dataa <- na.omit(dane)
dataa <- dataa[!duplicated(dataa$case_id), ]

view(dane[1:20, 100:ncol(dane)])

for(i in 1:nrow(dataa)){
  if(dataa[i,16] == 'Dead'){
    dataa[i,16] = 0
  }
  else{
    dataa[i,16] = 1
  }
}

dataa$vital_status = as.numeric(dataa$vital_status)

for(i in 1:nrow(dataa)){
  if(dataa[i,16] == 1){
    dataa[i,10] = dataa[i,50]
  }
}

dataa <- subset(dataa, dataa$days_to_death != '\'--')

for(i in 1:nrow(dataa)){
  if(dataa[i,16] == 1){
    dataa[i,10] = dataa[i,50]
  }
}

dataa$days_to_death = as.numeric(dataa$days_to_death)
dataa$age_at_index = as.numeric(dataa$age_at_index)
dataa <- subset(dataa, dataa$days_to_death >= 0)
dataa <- subset(dataa, dataa$ajcc_pathologic_stage != '\'--')

dataa <- dataa[, !sapply(dataa, function(x) any(x == '\'--'))]




#KAPLAN

data_kaplan = dataa
data_kaplan[, sapply(data_kaplan, is.character)] <- lapply(data_kaplan[, sapply(data_kaplan, is.character)], as.factor)

mayer = survfit(Surv(days_to_death, vital_status)~ajcc_pathologic_stage, data = data_kaplan)
summary(mayer)

plot(mayer, xlab = 'Czas(w dniach)', ylab = 'Prawdopodobienstwo smierci')
ggsurvplot(mayer, xlab = 'Czas(w dniach)', ylab = 'Prawdopodobienstwo przezycia', color = 'ajcc_pathologic_stage')

  #kaplan - przezywalnosc od rodzaju terapii
  kaplan_terapia <- survfit(Surv(days_to_death, vital_status) ~ treatment_type, data=data_kaplan)
  autoplot(kaplan_terapia)
  ggsurvplot(kaplan_terapia, xlab = 'Czas(w dniach)', ylab = 'Prawdopodobienstwo przezycia', palette = 'treatment_type')
  
  #kaplan - przezywalnosc od plci
  kaplan_plec <- survfit(Surv(days_to_death, vital_status) ~ gender, data=data_kaplan)
  autoplot(kaplan_plec)
  ggsurvplot(kaplan_plec, xlab = 'Czas(w dniach)', ylab = 'Prawdopodobienstwo przezycia', palette = 'vital_status')
  
  #kaplan - czy bylo priorytetowe leczenie ###USUNIETE W TYM PRZYPADKU WIERSZE Z NOT REPORTED
  data_kaplan_notreported = subset(data_kaplan, data_kaplan$prior_treatment != 'Not Reported')
  kaplan_prio <- survfit(Surv(days_to_death, vital_status) ~ prior_treatment, data=data_kaplan_notreported)
  autoplot(kaplan_prio)
  ggsurvplot(kaplan_prio, xlab = 'Czas(w dniach)', ylab = 'Prawdopodobienstwo przezycia', palette = 'prior_treatment')
  
  #kaplan - diagnoza ### DANE Z 'NA' DO WYWALENIA, ALE ICH NIE MA
  kaplan_diagnoza <- survfit(Surv(days_to_death, vital_status) ~ primary_diagnosis, data=data_kaplan_diagnoza)
  ggsurvplot(kaplan_diagnoza, xlab = 'Czas(w dniach)', ylab = 'Prawdopodobienstwo przezycia', color = 'primary_diagnosis')
  
  
#COX

library(survival)
library(survminer)
library(ggfortify)

data_cox = dataa
data_cox[, sapply(data_cox, is.character)] <- lapply(data_cox[, sapply(data_cox, is.character)], as.factor)
cox = coxph(Surv(days_to_death, vital_status) ~ age_at_index + ajcc_pathologic_stage + treatment_type, data = data_cox)
fit = survfit(cox, data = data_cox)
summary(cox)
ggsurvplot(fit)
autoplot(fit)

  #duzo wykresow
  test22 = aareg(Surv(days_to_death, vital_status) ~ age_at_index + ajcc_pathologic_stage + treatment_type, data = data_cox)
  autoplot(test22)

  

#RANDOM FOREST
  
library(survival)
library(caret)  
library(randomForestSRC)

data_las = dataa

data_las[, sapply(data_las, is.character)] <- lapply(data_las[, sapply(data_las, is.character)], as.factor)

set.seed(2505)


ID = seq(from = 1, to = nrow(data_las), by = 1)
new.data = cbind(ID,data_las)
ID.train = sample(new.data$ID,723)
ID.test = ID[-ID.train]
train_data = new.data[ID.train,]
test_data = new.data[ID.test,]

#train = sample(1:1084, 723, replace = FALSE)
#train_data = data_las[train, ]
#test_data = data_las[-train, ]

las = rfsrc(Surv(days_to_death, vital_status) ~ ajcc_pathologic_stage + age_at_index + treatment_type , data = train_data, ntree = 1500, mtry = 2, importance = TRUE)
prediction = predict(las, test_data)
plot.survival(prediction)




