#This script enables one to use correctly formatted MVMR beta files built in an MVMR script and apply MR Egger methods and Q tests. Author: Daniel Jones September 2019
#######################################################################################################################


######################################################################################################################################
#1. Load packages
######################################################################################################################################
rm(list=ls(all=TRUE)) #empties your R environment

#Install Packages
install.packages("MendelianRandomization")
library(MendelianRandomization)
install.packages("data.table")
library(data.table)

######################################################################################################################################
#2. Perform MR Egger and Q tests
######################################################################################################################################

##For EA and MH on CVD
#Read in file
mv1<-fread("EA-MH-CVD-mvmr_betas.csv")

#Define vectors
gy<-c(mv1$CVD.beta)
gy_se<-c(mv1$CVD.SE)
gx<-c(mv1$EA.beta)
gx_se<-c(mv1$EA.SE)
gz<-c(mv1$MH.beta)
gz_se<-c(mv1$MH.SE)

#For MVMR Egger
mv1<-lm(gy~gx+gz, weights=gy_se^-2)
summary(mv1)

#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.0010381  0.0005319   1.952   0.0515    
#EA           -0.5094456  0.0464266 -10.973   <2e-16 ***
#MH            0.3703038  0.2206429   1.678   0.0938
#Residual standard error: 1.111 on 570 degrees of freedom
#Multiple R-squared:  0.1759,	Adjusted R-squared:  0.173 
#F-statistic:   60.81 on 2 and 570 DF,  p-value: < 2.2e-16
                       
#For Q-stats
#Create weights
fitted<-predict(mv1)
weight = gy_se^2+((gz^2)*(gz_se^2))+((gx^2)*(gx_se^2))

#Perform Q stat
Q<- sum(((gy-fitted)^2)/weight)
str(Q)
pchisq(Q, 570, lower.tail=FALSE) #This gives a p-value. 1269 is degrees of freedom
#Q=703
#P-value=0.0001130101

##For EA and MDD
#Read in file
mv2<-fread("EA-MDD-mvmr_betas.csv")

#Define vectors
gy<-c(mv2$CVD.beta)
gy_se<-c(mv2$CVD.SE)
gx<-c(mv2$EA.beta)
gx_se<-c(mv2$EA.SE)
gz<-c(mv2$MDD.beta)
gz_se<-c(mv2$MDD.SE)

#For MVMR Egger
mv2<-lm(gy~gx+gz, weights=gy_se^-2)
summary(mv2)

#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.0005548  0.0003571   1.553  0.12056    
#EA           -0.4539457  0.0329343 -13.783  < 2e-16 ***
#MDD           0.0774048  0.0254967   3.036  0.00245
#Residual standard error: 1.11 on 1237 degrees of freedom
#Multiple R-squared:  0.1602,	Adjusted R-squared:  0.1588 
#F-statistic:   117.9 on 2 and 1237 DF,  p-value: < 2.2e-16
                      
#For Q-stats
#Create weights
fitted<-predict(mv2)
weight = gy_se^2+((gz^2)*(gz_se^2))+((gx^2)*(gx_se^2))

#Perform Q stat
Q<- sum(((gy-fitted)^2)/weight)
str(Q)
pchisq(Q, 1237, lower.tail=FALSE) #This gives a p-value. 1237 is degrees of freedom
#Q=1524
#P-value=3.40e-08

##For EA and Anxiety
#Read in file
mv3<-fread("EA-Anxiety-mvmr_betas.csv")

#Define vectors
gy<-c(mv3$CVD.beta)
gy_se<-c(mv3$CVD.SE)
gx<-c(mv3$EA.beta)
gx_se<-c(mv3$EA.SE)
gz<-c(mv3$Anx.beta)
gz_se<-c(mv3$anx.SE)

#For MVMR Egger
mv3<-lm(gy~gx+gz, weights=gy_se^-2)
summary(mv3)

#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.0004047  0.0003537   1.144    0.253    
#EA           -0.4853889  0.0327542 -14.819   <2e-16 ***
#Anx           0.0046622  0.0065431   0.713    0.476
#Residual standard error: 1.11 on 1237 degrees of freedom
#Multiple R-squared:  0.1561,	Adjusted R-squared:  0.1548 
#F-statistic:   113.4 on 2 and 1226 DF,  p-value: < 2.2e-16
           
#For Q-stats
#Create weights
fitted<-predict(mv3)
weight = gy_se^2+((gz^2)*(gz_se^2))+((gx^2)*(gx_se^2))

#Perform Q stat
Q<- sum(((gy-fitted)^2)/weight)
str(Q)
pchisq(Q, 1226, lower.tail=FALSE) #This gives a p-value. 1226 is degrees of freedom
#Q=0
#P-value=1

##For MDD and LSI
#Read in file
mv4<-fread("MDD-LSI-mvmr_betas.csv")

#Define vectors
gy<-c(mv4$CVD.beta)
gy_se<-c(mv4$CVD.SE)
gx<-c(mv4$MDD.beta)
gx_se<-c(mv4$MDD.SE)
gz<-c(mv4$LSI.beta)
gz_se<-c(mv4$LSI.SE)

#For MVMR Egger
mv4<-lm(gy~gx+gz, weights=gy_se^-2)
summary(mv4)

#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.002360   0.001018   2.317   0.0219    
#MDD           0.050647   0.052440   0.966   0.3358
#LSI           0.625956   0.125719   4.979   1.84e-06 ***
#Residual standard error: 1.159 on 141 degrees of freedom
#Multiple R-squared:  0.219,	Adjusted R-squared:  0.2079 
#F-statistic:   19.76 on 2 and 141 DF,  p-value: 2.713e-08
                 
#For Q-stats
#Create weights
fitted<-predict(mv4)
weight = gy_se^2+((gz^2)*(gz_se^2))+((gx^2)*(gx_se^2))

#Perform Q stat
Q<- sum(((gy-fitted)^2)/weight)
str(Q)
pchisq(Q, 141, lower.tail=FALSE) #This gives a p-value. 141 is degrees of freedom
#Q=189
#P-value=0.00423

##For EA and MDD on LSI
#Read in file
mv5<-fread("MDD-EA-LSI-mvmr_betas.csv")

#Define vectors
gy<-c(mv5$LSI.beta)
gy_se<-c(mv5$LSI.SE)
gx<-c(mv5$MDD.beta)
gx_se<-c(mv5$MDD.SE)
gz<-c(mv5$EA.beta)
gz_se<-c(mv5$EA.SE)

#For MVMR Egger
mv5<-lm(gy~gx+gz, weights=gy_se^-2)
summary(mv5)

#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -2.043e-05  8.509e-05  -0.240     0.81    
#MDD           5.255e-02  6.041e-03   8.699   <2e-16
#EA           -2.213e-01  7.812e-03 -28.332   <2e-16 
#Residual standard error: 1.706 on 1084 degrees of freedom
#Multiple R-squared:  0.4962,	Adjusted R-squared:  0.4953 
#F-statistic:  533.9 on 2 and 1084 DF,  p-value: < 2.2e-16

#For Q-stats
#Create weights
fitted<-predict(mv5)
weight = gy_se^2+((gz^2)*(gz_se^2))+((gx^2)*(gx_se^2))

#Perform Q stat
Q<- sum(((gy-fitted)^2)/weight)
str(Q)
pchisq(Q, 1084, lower.tail=FALSE) #This gives a p-value. 1084 is degrees of freedom
#Q=3129
#P-value=3.027558e-197
