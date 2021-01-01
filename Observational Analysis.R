#This script allows one to run univariate and multivariate regression on UK biobank data
#Author: Originally Robyn Wootton October 2018. Modified by Daniel Jones October 2019
#######################################################################################################################
#Contents
#1. Read in the data
#2. Tidy up data - name variables appropriately
#3. Exclude participants who have withdrawn consent
#4. Read in SNP and other smoking variables
#5. Split into ever and never smokers
#6. Get the demographics per genotype of the SNP
#7. Test Assosiation of genetic score with suicide in smokers
#8. Test Assosiation of genetic score with suicide in non-smokers 
#9. Plot the results

##########################################################################
#1. Read in data
##########################################################################
rm(list=ls(all=TRUE)) #empties your R environment

##Connect to IEU server using finder-go-connect to server
#Read in ICD code data
library(readstata13)
read <- read.dta13("20170613_sleep_icd10.dta")
str(read)

#Identify the ICD codes using UK Biobank variables [main diagosis = 41202/ secondary diagnosis = 41204]
diagnosis_main<-read[ , grepl( "41202" , names( read ) ) ]
str(diagnosis_main)
table(diagnosis_main$s_41202_0_10, exclude=NULL)
diagnosis_second<-read[ , grepl( "41204" , names( read ) ) ]
str(diagnosis_second)

#Stick all of the datasets together respectively
data<-data.frame(diagnosis_main, diagnosis_second)
str(data)

#Read in phenotype data
library(readstata13)
read2 <- read.dta13("20160609_7152_allphenotype.dta")
str(read2)

# ##########################################################################
# #2. Tidy up data - name variables appropriately
# ##########################################################################

###For ICD coded cases 
#Set all blank values as missing
data[c(2:816)][data[c(2:816)] == ""] <-NA
str(data)

#Restrict all of the ICD 10 codes to the first two numbers only
for (i in 2:816){
	data[, i] <- substr(data[, i], 0, 5)
 }
#Save dataset
write.csv(data, "ICD10_codes.csv", quote=F, row.names=F)
data<-read.csv("ICD10_codes.csv", header=T)

#Now sort data for CVD
CVD <- data
for (i in 2:816){
	CVD[, i] <- ifelse(CVD[, i]!="I20" & CVD[, i]!="I21" & CVD[, i]!="I63", 0, 1)
}
table(CVD $s_41202_0_10) #N=575 
CVD$cvd_sum <- rowSums(CVD[c(2:816)], na.rm=T)
table(CVD$cvd_sum)
#make all values above 1 into 1
CVD$cvd[CVD$cvd_sum>1] <-1
table(CVD$cvd)

#Now depression
DEP <- data
for (i in 2:816){
	DEP[, i] <- ifelse(DEP[, i]!="F32" & DEP[, i]!="F33", 0, 1)
}
table(DEP $s_41202_0_10) #N=54 
DEP$dep_sum <- rowSums(DEP[c(2:816)], na.rm=T)
table(DEP$dep_sum)
#make all values above 1 into 1
DEP$dep[DEP$dep_sum>1] <-1
table(DEP$dep) 

#Now anxiety
ANXIETY <- data
for (i in 2:816){
	ANXIETY[, i] <- ifelse(ANXIETY[, i]!="F41", 0, 1)
}
table(ANXIETY $s_41202_0_10) #N=26 
ANXIETY$anxiety_sum <- rowSums(ANXIETY[c(2:816)], na.rm=T)
table(ANXIETY$anxiety_sum)
#make all values above 1 into 1
ANXIETY$anxiety[ANXIETY$anxiety_sum>1] <-1
table(ANXIETY$anxiety) 

#Now HTN
HTN <- data
for (i in 2:816){
	HTN[, i] <- ifelse(HTN[, i]!="I10", 0, 1)
}
table(HTN $s_41202_0_10) #N=32 
HTN$htn_sum <- rowSums(HTN[c(2:816)], na.rm=T)
table(HTN$htn_sum)
#make all values above 1 into 1
HTN$htn[HTN$htn_sum>1] <-1
table(HTN$htn) 

#Now obesity
OBS <- data
for (i in 2:816){
	OBS[, i] <- ifelse(OBS[, i]!="E66", 0, 1)
}
table(OBS $s_41202_0_10) #N=29 
OBS$obs_sum <- rowSums(OBS[c(2:816)], na.rm=T)
table(OBS$obs_sum)
#make all values above 1 into 1
OBS$obs[OBS$obs_sum>1] <-1
table(OBS$obs)

#Subset into the relevant columns
CVD2 <- subset(CVD, select=c("read.n_eid", "cvd"))
str(CVD2)
#N=502621
DEP <- subset(DEP, select=c("read.n_eid", "dep"))
str(DEP)
#N=502621
ANXIETY <- subset(ANXIETY, select=c("read.n_eid", "anxiety"))
str(ANXIETY)
#N=502621
HTN <- subset(HTN, select=c("read.n_eid", "htn"))
str(HTN)
#N=502621
OBS <- subset(OBS, select=c("read.n_eid", "obs"))
str(OBS)
#N=502621

#Merge ICD data together
final <- merge(CVD2, DEP, by="read.n_eid", all=T)
final<- merge(final, ANXIETY, by="read.n_eid", all=T)
final <- merge(final, HTN, by="read.n_eid", all=T)
final<- merge(final, OBS, by="read.n_eid", all=T)
str(final) #N=502621

#Save dataset
setwd("/Volumes/UK_Biobank_Project_9142/Working Folder Daniel")
write.csv(final, "UKBB_ICD10variables.csv", quote=F, row.names=F)
CVD_Dep <- read.csv("UKBB_ICD10variables.csv", header=T)

###For exposure and mediator variables
#Identify required exposure variables by UK biobank variables
ID<-read2[ , grepl( "n_eid" , names( read2 ) ) ]
str(ID)
#Seen a psychiatrist for nerves, anxiety, tension, depression = 2100
seen_psych<-read2[ , grepl( "2100" , names( read2 ) ) ]
str(seen_psych)
#Seen GP for nerves, anxiety, tension, depression = 2100
seen_gp<-read2[ , grepl( "2090" , names( read2 ) ) ]
str(seen_gp)
#Neuroticism - Tense/Highly strung? = 1990
neuro<-read2[ , grepl( "1990" , names( read2 ) ) ]
str(neuro)
#Worry - Do you have nerves? = 2010
worry<-read2[ , grepl( "2010" , names( read2 ) ) ]
str(worry)
#Depressed affect - Miserable? = 1930
lowmood<-read2[ , grepl( "1930" , names( read2 ) ) ]
str(lowmood)
#Education - Qual = 6138
ea<-read2[ , grepl( "6138" , names ( read2 ) ) ]
str(ea)

#Create a phenotype file
data2<-data.frame(ID, seen_psych, seen_gp, neuro, worry, lowmood, ea)
str(data2)

#Save dataset
write.csv(data2, "Phenotypes.csv", quote=F, row.names=F)
data2<-read.csv("Phenotypes.csv", header=T)

#Set all blank values as missing
data2[c(2:175)][data2[c(2:175)] == ""] <-NA
str(data2)

#Rename columns
colnames(data2)[colnames(data2)=="n_2100_0_0"]<-"Seen_Psych"
colnames(data2)[colnames(data2)=="n_2090_0_0"]<-"Seen_GP"
colnames(data2)[colnames(data2)=="n_1990_0_0"]<-"Neuroticism"
colnames(data2)[colnames(data2)=="n_2010_0_0"]<-"Worry"
colnames(data2)[colnames(data2)=="n_1930_0_0"]<-"Low_Mood"
colnames(data2)[colnames(data2)=="n_6138_0_0"]<-"Educational_Attainment"

#Select these columns out
data2 = data2[c("ID", "Seen_Psych", "Seen_GP", "Neuroticism", "Worry", "Low_Mood", "Educational_Attainment",)]
str(data2)
write.csv(data2, "PhenotypesII.csv", quote=F, row.names=F)

#Finally
rm(list=ls(all=TRUE))
data1 <- read.csv("UKBB_ICD10variables.csv", header=T)
data2<-read.csv("PhenotypesII.csv", header=T)

#########################################################################
#3. Exclude participants who have withdrawn consent
#########################################################################
#Read in the exclusion list
ex<-read.csv("w9142_20181016.csv", header=F)
str(ex)
ex$exclude<-1

#Merge these into the data files
df<-merge(data1, ex, by.x="read.n_eid", by.y="V1", all.x=T)
df<-merge(df, data2, by.x="read.n_eid", by.y="ID", all=T)
str(df) #N= 510073
dat<-subset(df, is.na(df$exclude))
str(dat) #N=509,995

######################################################################################
#4. Read in Smoking data and perform exclusion
######################################################################################
#For smoking data
data<-read.csv("Feb2018_smoking_mortality.csv", header=T)
str(data)

#Merge with the previous data
final <- merge(dat, data, by.x="read.n_eid", by.y="ID", all.x=T)
head(final)

#Apply the genetic exclusions
table(final $snp) #0=151220, 1=149076, 2=36754

#Exclude people who are not european ancestry using IEU guide (DOI: 10.5523/bris.3074krb6t2frj29yh2b03x3wxj)
#Read in the genetic data from BMI - here is the 500k with exclusions applied and we can subset on this
scores<-read.csv("UKBiobank_BMI_PRS.csv", header=T) 
str(scores) #Check N=337,115

#Read in the ID linker file
library(readstata13)
mergelist <- read.dta13("matching_id.dta")
str(mergelist)

data2 <-merge(final, mergelist, by.x="read.n_eid", by.y="appid9142", all.x=T)
str(data2)

df<-merge(data2, scores, by.x="appid8786", by.y="IID", all.x=T)
str(df) #N=509995

#Now subset only those who have a score for the 500k
data3<-subset(df, complete.cases(df$SCORE500))
str(data3) #N = 337, 058
table(data3$cvd) #Checking UKBB population prevalence = 6.07%
table(data3$dep) #UKBB prevalence = 2.86%
table(data3$anxiety) #UKBB prevalence for anxiety = 1.42%
table(data3$htn)
table(data3$obs)
table(data3$Seen_Psych) #UKBB prevalence = 11.51%
table(data3$Seen_GP) #UKBB prevalence = 33.79%
table(data3$Neuroticism) #UKBB prevalence = 16.63%
table(data3$Worry) #UKBB prevalence = 20.72%
table(data3$Low_Mood) #UKBB prevalence = 42.05%

#Save final data file
write.csv(data3, "FullData.csv", quote=F, row.names=F)
data3<-read.csv("FullData.csv", header=T)

######################################################################################
#5. Sort by smoking status
######################################################################################
#Smoking status
ever<-subset(data3, data3 $nevev=="Ever")
str(ever) #N= 151, 822

never<-subset(data3, data3 $nevev =="Never")
str(never) #N= 184044

current<-subset(data3, data3 $smoking_status =="Current")
str(current) #N= 33357

former <- subset(data3, data3 $smoking_status =="Previous")
str(former) #N= 118465

#Descriptive stats
table(data3$sex) #Female = 181325, Male = 155725
table(data3$ageatstart)
mean(data3$ageatstart, na.rm=T) #56.9
sd(data3$ageatstart, na.rm=T) #8.00
mean(data3$SES, na.rm=T) #-1.58
sd(data3$SES, na.rm=T) #2.928618

######################################################################################
#6. Run  the observational associations (controlling for age, sex and SES)
######################################################################################
## First flip the levels of the variables the right way
gen<-data3

#Never/Ever smoking
gen$nevev<-as.factor(gen$nevev)
levels(gen$nevev)
gen$nevev<-as.numeric(gen$nevev)
table(gen$nevev)
gen$nevev[gen$nevev==2]<-0
table(gen$nevev)
gen$nevev<-factor(gen$nevev, levels = c(0,1),labels = c("Never", "Ever"))
table(gen$nevev)

#Current/Previous smoking
gen$smoking_status<-as.factor(gen$smoking_status)
levels(gen$smoking_status)
table(gen$smoking_status)
gen$smoking_status <-as.numeric(gen$smoking_status)
table(gen$smoking_status)
gen$smoking_status[gen$smoking_status ==2]<-0
gen$smoking_status[gen$smoking_status ==1]<-2
gen$smoking_status[gen$smoking_status ==3]<-1
table(gen$smoking_status)
gen$smoking_status <-factor(gen$smoking_status, levels = c(0,1,2),labels = c("Never", "Previous", "Current"))
table(gen$smoking_status)

#Ever seen a psych?
gen$Seen_Psych <-as.factor(gen$Seen_Psych)
levels(gen$Seen_Psych)
table(gen$Seen_Psych)
gen$Seen_Psych <-as.numeric(gen$Seen_Psych)
table(gen$Seen_Psych)
gen$Seen_Psych[gen$Seen_Psych ==1]<-0
gen$Seen_Psych[gen$Seen_Psych ==2]<-0
gen$Seen_Psych[gen$Seen_Psych ==3]<-0
gen$Seen_Psych[gen$Seen_Psych ==4]<-1
table(gen$Seen_Psych)
gen$Seen_Psych <-factor(gen$Seen_Psych, levels = c(0,1),labels = c("No", "Yes"))
table(gen$Seen_Psych)

#Ever seen a GP?
gen$Seen_GP <-as.factor(gen$Seen_GP)
levels(gen$Seen_GP)
table(gen$Seen_GP)
gen$Seen_GP <-as.numeric(gen$Seen_GP)
table(gen$Seen_GP)
gen$Seen_GP[gen$Seen_GP ==1]<-0
gen$Seen_GP[gen$Seen_GP ==2]<-0
gen$Seen_GP[gen$Seen_GP ==3]<-0
gen$Seen_GP[gen$Seen_GP ==4]<-1
table(gen$Seen_GP)
gen$Seen_GP <-factor(gen$Seen_GP, levels = c(0,1),labels = c("No", "Yes"))
table(gen$Seen_GP)

#Low mood
levels(gen$Low_Mood)
gen$Low_Mood <-as.numeric(gen$Low_Mood)
table(gen$Low_Mood)
gen$Low_Mood[gen$Low_Mood ==1]<-0
gen$Low_Mood[gen$Low_Mood ==2]<-0
gen$Low_Mood[gen$Low_Mood ==3]<-0
gen$Low_Mood[gen$Low_Mood ==4]<-1
table(gen$Low_Mood)
gen$Low_Mood <-factor(gen$Low_Mood, levels = c(0,1),labels = c("No", "Yes"))
table(gen$Low_Mood)
#No=195322 Yes=141736

#Tense/nervous
levels(gen$Neuroticism)
gen$Neuroticism <-as.numeric(gen$Neuroticism)
table(gen$Neuroticism)
gen$Neuroticism[gen$Neuroticism ==1]<-0
gen$Neuroticism[gen$Neuroticism ==2]<-0
gen$Neuroticism[gen$Neuroticism ==3]<-0
gen$Neuroticism[gen$Neuroticism ==4]<-1
table(gen$Neuroticism)
gen$Neuroticism <-factor(gen$Neuroticism, levels = c(0,1),labels = c("No", "Yes"))
table(gen$Neuroticism)
#No=280991 Yes=56067

#Worry
levels(gen$Worry)
gen$Worry <-as.numeric(gen$Worry)
table(gen$Worry)
gen$Worry[gen$Worry ==1]<-0
gen$Worry[gen$Worry ==2]<-0
gen$Worry[gen$Worry ==3]<-0
gen$Worry[gen$Worry ==4]<-1
table(gen$Worry)
gen$Worry <-factor(gen$Worry, levels = c(0,1),labels = c("No", "Yes"))
table(gen$Worry)
#No=267576 Yes=69842

#EA
gen$Educational_Attainment<-as.factor(gen$Educational_Attainment)
levels(gen$Educational_Attainment)
gen$Educational_Attainment<-as.numeric(gen$Educational_Attainment)
table(gen$Educational_Attainment)
gen$Educational_Attainment[gen$Educational_Attainment ==1]<-13
gen$Educational_Attainment[gen$Educational_Attainment ==2]<-20
gen$Educational_Attainment[gen$Educational_Attainment ==3]<-10
gen$Educational_Attainment[gen$Educational_Attainment ==5]<-19
gen$Educational_Attainment[gen$Educational_Attainment ==6]<-10
gen$Educational_Attainment[gen$Educational_Attainment ==7]<-15
gen$Educational_Attainment[gen$Educational_Attainment ==4]<-7
gen$Educational_Attainment[gen$Educational_Attainment ==8]<-NA
table(gen$Educational_Attainment)

#7=57099, 10=92248, 13=38440, 15=17290, 19=22105, 20=106748
mean(gen$Educational_Attainment, na.rm=T) #13.9
sd(gen$Educational_Attainment, na.rm=T) #5.1

#Create dummy codes
gen$Educational_Attainment[gen$Educational_Attainment ==7]<-0
gen$Educational_Attainment[gen$Educational_Attainment ==10]<-0
gen$Educational_Attainment[gen$Educational_Attainment ==13]<-1
gen$Educational_Attainment[gen$Educational_Attainment ==15]<-1
gen$Educational_Attainment[gen$Educational_Attainment ==19]<-2
gen$Educational_Attainment[gen$Educational_Attainment ==20]<-2
table(gen$Educational_Attainment)

#create ever again with these changes
ever<-subset(gen, gen$nevev=="Ever")
str(ever) #N= 151822

#Run regression models
#Education on CVD
model1<- glm(cvd ~ Educational_Attainment + ageatstart + sex + SES, data=gen, family=binomial(link="logit"))
mod1<-summary(model1)
or1<-mod1$coefficients[2,1]
se1<-mod1$coefficients[2,2]
p1<-mod1$coefficients[2,4]
n1<-mod1$df[2]
row1<-cbind(or1, se1, p1, n1)
row1

#Education on Seen Psych?
model2<- glm(Seen_Psych ~ Educational_Attainment + ageatstart + sex + SES, data=gen, family=binomial(link="logit"))
mod2<-summary(model2)
or2<-mod2$coefficients[2,1]
se2<-mod2$coefficients[2,2]
p2<-mod2$coefficients[2,4]
n2<-mod2$df[2]
row2<-cbind(or2, se2, p2, n2)
row2

#Education on Seen GP?
model3<- glm(Seen_GP ~ Educational_Attainment + ageatstart + sex + SES, data=gen, family=binomial(link="logit"))
mod3<-summary(model3)
or3<-mod3$coefficients[2,1]
se3<-mod3$coefficients[2,2]
p3<-mod3$coefficients[2,4]
n3<-mod3$df[2]
row3<-cbind(or3, se3, p3, n3)
row3

#Education on Depression?
model4<- glm(dep ~ Educational_Attainment + ageatstart + sex + SES, data=gen, family=binomial(link="logit"))
mod4<-summary(model4)
or4<-mod4$coefficients[2,1]
se4<-mod4$coefficients[2,2]
p4<-mod4$coefficients[2,4]
n4<-mod4$df[2]
row4<-cbind(or4, se4, p4, n4)
row4

#Education on Anxiety?
model5<- glm(anxiety ~ Educational_Attainment + ageatstart + sex + SES, data=gen, family=binomial(link="logit"))
mod5<-summary(model5)
or5<-mod5$coefficients[2,1]
se5<-mod5$coefficients[2,2]
p5<-mod5$coefficients[2,4]
n5<-mod5$df[2]
row5<-cbind(or5, se5, p5, n5)
row5

#Seeing a Psych on CVD?
model6<- glm(cvd ~ Seen_Psych + Educational_Attainment + ageatstart + sex + smoking_status + htn + obs + SES , data=gen, family=binomial(link="logit"))
mod6<-summary(model6)
or6<-mod6$coefficients[2,1]
se6<-mod6$coefficients[2,2]
p6<-mod6$coefficients[2,4]
n6<-mod6$df[2]
row6<-cbind(or6, se6, p6, n6)
row6

#Seeing a GP on CVD?
model7<- glm(cvd ~ Seen_GP + Educational_Attainment + ageatstart + sex + smoking_status + htn + obs + SES , data=gen, family=binomial(link="logit"))
mod7<-summary(model7)
or7<-mod7$coefficients[2,1]
se7<-mod7$coefficients[2,2]
p7<-mod7$coefficients[2,4]
n7<-mod7$df[2]
row7<-cbind(or7, se7, p7, n7)
row7

#Depression on CVD?
model8<- glm(cvd ~ dep + anxiety + Educational_Attainment + ageatstart + sex + smoking_status + htn + obs + SES , data=gen, family=binomial(link="logit"))
mod8<-summary(model8)
or8<-mod8$coefficients[2,1]
se8<-mod8$coefficients[2,2]
p8<-mod8$coefficients[2,4]
n8<-mod8$df[2]
row8<-cbind(or8, se8, p8, n8)
row8

#Anxiety on CVD?
model9<- glm(cvd ~ anxiety + dep + Educational_Attainment + ageatstart + sex + smoking_status + htn + obs + SES, data=gen, family=binomial(link="logit"))
mod9<-summary(model9)
or9<-mod9$coefficients[2,1]
se9<-mod9$coefficients[2,2]
p9<-mod9$coefficients[2,4]
n9<-mod9$df[2]
row9<-cbind(or9, se9, p9, n9)
row9

#Education on Low mood?
model10<- glm(Low_Mood ~ edu + ageatstart + sex + SES, data=gen, family=binomial(link="logit"))
mod10<-summary(model10)
or10<-mod10$coefficients[2,1]
se10<-mod10$coefficients[2,2]
p10<-mod10$coefficients[2,4]
n10<-mod10$df[2]
row10<-cbind(or10, se10, p10, n10)
row10

#Education on Worry?
model12<- glm(Worry ~ edu + ageatstart + sex + SES, data=gen, family=binomial(link="logit"))
mod12<-summary(model12)
or12<-mod12$coefficients[2,1]
se12<-mod12$coefficients[2,2]
p12<-mod12$coefficients[2,4]
n12<-mod12$df[2]
row12<-cbind(or12, se12, p12, n12)
row12

#Low mood on CVD, adj for worry
model13<- glm(cvd ~ Low_Mood + Worry + edu + ageatstart + sex + smoking_status + htn + obs + SES, data=gen, family=binomial(link="logit"))
mod13<-summary(model13)
or13<-mod13$coefficients[2,1]
se13<-mod13$coefficients[2,2]
p13<-mod13$coefficients[2,4]
n13<-mod13$df[2]
row13<-cbind(or13, se13, p13, n13)
row13

#Worry on CVD, adj for low mood
model14<- glm(cvd ~ Worry + Low_Mood + edu + ageatstart + sex + smoking_status + htn + obs + SES, data=gen, family=binomial(link="logit"))
mod14<-summary(model14)
or14<-mod14$coefficients[2,1]
se14<-mod14$coefficients[2,2]
p14<-mod14$coefficients[2,4]
n14<-mod14$df[2]
row14<-cbind(or14, se14, p14, n14)
row14

###CSI & MDD models
#Merge in CSI
csi<-read.table("lifetimesmoking_phen_revision1.txt", header=T)
str(csi)
csi$csi <- as.numeric(csi$csi)
gen_csi<-merge(gen, csi, by.x="FID", by.y="IID", all.x=T)

#MDD on LSI (Linear regression)
model15<- lm(csi ~ dep + ageatstart + Educational_Attainment + sex +SES, data=gen_csi)
mod15<-summary(model15)
or15<-mod15$coefficients[2,1]
se15<-mod15$coefficients[2,2]
p15<-mod15$coefficients[2,4]
n15<-mod15$df[2]
row15<-cbind(or15, se15, p15, n15)
row15

#LSI on CVD
model16<- glm(cvd ~ csi + dep + ageatstart + sex + htn + obs + SES, data=gen_csi, family=binomial(link="logit"))
mod16<-summary(model16)
or16<-mod16$coefficients[2,1]
se16<-mod16$coefficients[2,2]
p16<-mod16$coefficients[2,4]
n16<-mod16$df[2]
row16<-cbind(or16, se16, p16, n16)
row16

#stick all of the log odds together
results<-rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12, row13, row14, row15, row16)
results<-data.frame(results)
exposure<-c("Education","Education","Education", "Education", "Education","Seen_Psych","Seen_GP", "Depression","Anxiety","Education", "Education","Education", "Low_Mood", "Worry", "Neuroticism", "Neuroticism-Low_Mood", "Neuroticism-Worry")
outcome<-c("CVD", "Seen_Psych", "Seen_GP", "Depression", "Anxiety", "CVD","CVD", "CVD", "CVD","Low_Mood", "Neuroticism", "Worry", "CVD", "CVD", "CVD", "CVD", "CVD")
results2<-cbind(exposure, outcome, results)
names(results2)<-c("exposure", "outcome", "b", "se", "pval", "N")
str(results2)

#convert to odds ratios
library(TwoSampleMR)
or<-generate_odds_ratios(results2)
or

#save the dataset
write.csv(or, "Results.csv", quote=F, row.names=F)
