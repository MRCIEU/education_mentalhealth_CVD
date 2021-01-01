#Analysis conducted by Daniel Jones
#Based on script by Robyn Wootton
#April 2019
######################################################################################################################
#Contents
#1. Load packages
#2. Read in xposure data
#3. Read in outcome data
#4. Harmonise data
#5. Run 2 sample MR
#6. Steiger filtering
#7. Plot results
#8. Regression dilution I2 GX
#9. Simex correction 

######################################################################################################################################
#1. Load packages
######################################################################################################################################
rm(list=ls(all=TRUE))

install.packages("devtools")
install_github("MRCIEU/TwoSampleMR")
install.packages("ggplot2")
install.packages("knitr")
install_github('qingyuanzhao/mr.raps')
install.packages("Cairo")
install.packages("curl")

library(devtools)
library(TwoSampleMR)
library(ggplot2)
library(knitr)
library(mr.raps)
library(Cairo)
library(curl)

######################################################################################################################################
#2. Read in instruments
######################################################################################################################################

#Educational Attainment
ea_exp_dat <- read_exposure_data(
       filename = "EA-Lead_SNPs.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "Pval",
   )
ea_exp_dat$exposure<-"Educational Attainment"

#Anxiety
anx_exp_dat <- read_exposure_data(
       filename = "Anxiety-Instrument.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Effect",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "P",
   )
anx_exp_dat$exposure<-"Anxiety"
anx_exp_dat<-clump_data(anx_exp_dat, clump_r2 = 0.01)

#Depression
dep_exp_dat <- read_exposure_data(
       filename = "MDD-Instrument.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "OR",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       pval_col = "P",
       eaf_col= "EAF"
   )
dep_exp_dat$exposure<-"Depression"

#Mental Health
mh_exp_dat <- read_exposure_data(
       filename = "MH-LeadSNPs.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "P",
   )
mh_exp_dat$exposure<-"Mental Health"

#Lifetime smoking
lsi_exp_dat <- read_exposure_data(
       filename = "CSI_genome-wideSNPs_MRBase.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "beta.exposure",
       se_col = "se.exposure",
       effect_allele_col = "effect_allele.exposure",
       other_allele_col = "other_allele.exposure",
       eaf_col = "eaf.exposure",
       pval_col = "pval.exposure",
   )
lsi_exp_dat$exposure<-"Lifetime Smoking"

#Worry
worry_exp_dat <- read_exposure_data(
       filename = "Worry_LeadSNPs.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "Pval",
   )
worry_exp_dat$exposure<-"Worry"
worry_exp_dat<-clump_data(worry_exp_dat)

#Neuroticism (No EAF available)
neuro_exp_dat <- read_exposure_data(
       filename = "Neuro-Leads.txt",
       sep = ",",
       snp_col = "rsID",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       pval_col = "Pval",
   )
neuro_exp_dat$exposure<-"Neuroticism"

#Depressed affect
dep2_exp_dat <- read_exposure_data(
       filename = "Neuro-DepressionLeads.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Beta",
       se_col = "SE",
       eaf_col= "EAF",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       pval_col = "Pval",
   )
dep2_exp_dat$exposure<-"Depressed Affect"
dep2_exp_dat<-clump_data(dep2_exp_dat)

######################################################################################################################################
#3. Read in relevant outcome data
######################################################################################################################################

#Educational Attainment
ea_outcome_dat <- read_outcome_data(
       filename = "GWAS_EA_excl23andMe.txt",
       sep = ",",
       snp_col = "MarkerName",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "Pval",
   )
ea_outcome_dat$outcome<-"Educational Attainment"

#Anxiety
anx_outcome_dat <- read_outcome_data(
       filename = "Anxiety(Sum).txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Effect",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "P",
   )
anx_outcome_dat$outcome<-"Anxiety"

#Depression
dep_outcome_dat <- read_outcome_data(
       filename = "MDD(Sum).txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "OR",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       pval_col = "P",
   )
dep_outcome_dat$outcome<-"Depression"

#CVD
cvd_outcome_dat <- read_outcome_data(
         filename = "cad1000genomes.txt",
         sep = ",",
         snp_col = "SNP",
         beta_col = "beta",
         se_col = "se",
         effect_allele_col = "A1",
         other_allele_col = "A2",
         eaf_col = "eaf",
         pval_col = "pval",
              )
cvd_outcome_dat $outcome<-"CVD"

#Mental Health
mh_outcome_dat <- read_outcome_data(
       filename = "MH_Summary.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "P",
   )
mh_outcome_dat$outcome<-"Mental Health"

#Lifetime smoking (LSI)
lsi_outcome_dat <- read_outcome_data(
       filename = "LSI_Summary.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "BETA",
       se_col = "SE",
       effect_allele_col = "EFFECT_ALLELE",
       other_allele_col = "OTHER_ALLELE",
       eaf_col = "EAF",
       pval_col = "P",
   )
lsi_outcome_dat$outcome<-"Lifetime Smoking"

#Worry
worry_outcome_dat <- read_outcome_data(
       filename = "Worry.txt",
       sep = ",",
       snp_col = "RSID",
       beta_col = "BETA",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "P",
   )
worry_outcome_dat$outcome<-"Worry"

#Neuroticism (No EAF available)
neuro_outcome_dat <- read_outcome_data(
       filename = "NEUROTICISM-full.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       pval_col = "Pval",
   )
neuro_outcome_dat$outcome<-"Neuroticism"

#Neuroticism-Depression
dep2_outcome_dat <- read_outcome_data(
       filename = "NeuroDepressionSumStats.txt",
       sep = ",",
       snp_col = "RSID",
       beta_col = "BETA",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "P",
   )
dep2_outcome_dat$outcome<-"Neuroticism-Depression"

######################################################################################################################################
#4. Harmonise data
######################################################################################################################################

#EA-CVD
dat1 <- harmonise_data( 
  	exposure_dat = ea_exp_dat,
  	outcome_dat = cvd_outcome_dat,
  	action = 2
  )

#EA-Dep
dat2 <- harmonise_data( 
  	exposure_dat = ea_exp_dat,
  	outcome_dat = dep_outcome_dat,
  	action = 2
  )

#Dep-CVD
dat3 <- harmonise_data( 
  	exposure_dat = dep_exp_dat,
  	outcome_dat = cvd_outcome_dat,
  	action = 2
  )


#EA-Anxiety
dat4 <- harmonise_data( 
  	exposure_dat = ea_exp_dat,
  	outcome_dat = anx_outcome_dat,
  	action = 2
  )

#Anxiety-CVD
dat5 <- harmonise_data( 
  	exposure_dat = anx_exp_dat,
  	outcome_dat = cvd_outcome_dat,
  	action = 2
  )

#EA-MH
dat6 <- harmonise_data( 
  	exposure_dat = ea_exp_dat,
  	outcome_dat = mh_outcome_dat,
  	action = 2
  )

#MH-CVD
dat7 <- harmonise_data( 
  	exposure_dat = mh_exp_dat,
  	outcome_dat = cvd_outcome_dat,
  	action = 2
  )

#EA-Worry
dat8 <- harmonise_data( 
  	exposure_dat = ea_exp_dat,
  	outcome_dat = worry_outcome_dat,
  	action = 2
  )
  
#Worry-CVD
dat9 <- harmonise_data( 
  	exposure_dat = worry_exp_dat,
  	outcome_dat = cvd_outcome_dat,
  	action = 2
  )  

#EA-Depressed Affect
dat10 <- harmonise_data( 
  	exposure_dat = ea_exp_dat,
  	outcome_dat = dep2_outcome_dat,
  	action = 2
  )
  
#Depressed Affect-CVD
dat11 <- harmonise_data( 
  	exposure_dat = dep2_exp_dat,
  	outcome_dat = cvd_outcome_dat,
  	action = 2
  )  

#DEP-LSI
dat12 <- harmonise_data( 
  	exposure_dat = dep_exp_dat,
  	outcome_dat = lsi_outcome_dat,
  	action = 2
  )  

#LSI-CVD
dat13 <- harmonise_data( 
  	exposure_dat = lsi_exp_dat,
  	outcome_dat = cvd_outcome_dat,
  	action = 2
  )  

######################################################################################################################################
#5. Run MR
######################################################################################################################################
#EA-CVD
mr_het1 <- mr_heterogeneity(dat1)
str(mr_het1)
mr_ruck1 <- mr_rucker(dat1)
res1 <- mr(dat1, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int1 <- mr_pleiotropy_test(dat1)
res1

or1<-generate_odds_ratios(res1)
or1


#EA-Dep
mr_het2 <- mr_heterogeneity(dat2)
str(mr_het2)
mr_ruck2 <- mr_rucker(dat2)
res2 <- mr(dat2, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int2 <- mr_pleiotropy_test(dat2)
res2

or2<-generate_odds_ratios(res2)
or2

#Dep-CVD
mr_het3 <- mr_heterogeneity(dat3)
str(mr_het3)
mr_ruck3 <- mr_rucker(dat3)
res3 <- mr(dat3, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int3 <- mr_pleiotropy_test(dat3)
res3

or3<-generate_odds_ratios(res3)
or3

#EA-Anxiety
mr_het4 <- mr_heterogeneity(dat4)
str(mr_het4)
mr_ruck4 <- mr_rucker(dat4)
res4 <- mr(dat4, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int4 <- mr_pleiotropy_test(dat4)
res4

or4<-generate_odds_ratios(res4)
or4

#Anxiety-CVD
mr_het5 <- mr_heterogeneity(dat5)
str(mr_het5)
mr_ruck5 <- mr_rucker(dat5)
res5 <- mr(dat5, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int5 <- mr_pleiotropy_test(dat5)
res5

or5<-generate_odds_ratios(res5)
or5

#EA-MH
mr_het6 <- mr_heterogeneity(dat6)
str(mr_het6)
mr_ruck6 <- mr_rucker(dat6)
res6<- mr(dat6, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int6 <- mr_pleiotropy_test(dat6)
res6

or6<-generate_odds_ratios(res6)
or6


#MH-CVD
mr_het7 <- mr_heterogeneity(dat7)
str(mr_het7)
mr_ruck7 <- mr_rucker(dat7)
res7<- mr(dat7, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int7 <- mr_pleiotropy_test(dat7)
res7

or7<-generate_odds_ratios(res7)
or7

#merge results together and save
results <- rbind (or1, or2, or3, or4, or5, or6, or7)
write.csv(results, "Two-sample MR results - Education, MDD, Anxiety, MH & CVD.csv", quote=F, row.names=F)

#EA-Worry
mr_het8 <- mr_heterogeneity(dat8)
str(mr_het8)
mr_ruck8 <- mr_rucker(dat8)
res8<- mr(dat8, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int8 <- mr_pleiotropy_test(dat8)
res8

or8<-generate_odds_ratios(res8)
or8

#Worry-CVD
mr_het9 <- mr_heterogeneity(dat9)
str(mr_het9)
mr_ruck9 <- mr_rucker(dat9)
res9<- mr(dat9, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int9 <- mr_pleiotropy_test(dat9)
res9

or9<-generate_odds_ratios(res9)
or9

#EA-Depressed Affect
mr_het10 <- mr_heterogeneity(dat10)
str(mr_het10)
mr_ruck10 <- mr_rucker(dat10)
res10<- mr(dat10, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int10 <- mr_pleiotropy_test(dat10)
res10

or10<-generate_odds_ratios(res10)
or10


#Depressed Affect-CVD
mr_het11 <- mr_heterogeneity(dat11)
str(mr_het11)
mr_ruck11 <- mr_rucker(dat11)
res11<- mr(dat11, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int11 <- mr_pleiotropy_test(dat11)
res11

or11<-generate_odds_ratios(res11)
or11

#DEP-LSI
mr_het12 <- mr_heterogeneity(dat12)
str(mr_het12)
mr_ruck12 <- mr_rucker(dat12)
res12<- mr(dat12, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int12 <- mr_pleiotropy_test(dat12)
res12

or12<-generate_odds_ratios(res12)
or12

#LSI-CVD
mr_het13 <- mr_heterogeneity(dat13)
str(mr_het13)
mr_ruck13 <- mr_rucker(dat13)
res13<- mr(dat13, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int13 <- mr_pleiotropy_test(dat13)
res13

or13<-generate_odds_ratios(res13)
or13

#merge results together and save
results <- rbind (or12, or13)
write.csv(results, "Two-sample MR results - Smoking & CVD.csv", quote=F, row.names=F)

######################################################################################################################################
#6. Steiger Filtering
######################################################################################################################################

##Ensure relevant data harmonised

#For EA & CVD
dat1$samplesize.exposure<-1131881 
dat1$units.exposure<-"SD"
dat1$units.outcome <- "log odds"
dat1$samplesize.outcome<-184305 
dat1$ncase.outcome<-60801
dat1$ncontrol.outcome<-123504
dat1$prevalence.outcome<-0.076 
steiger <- steiger_filtering(dat1)
table(steiger$steiger_dir)

FALSE  TRUE 
   60   1206 

#For EA and Depression
dat2$samplesize.exposure<-766343 
dat2$units.exposure<-"SD"
dat2$units.outcome <- "log odds"
dat2$samplesize.outcome<-807553 
dat2$ncase.outcome<-135458
dat2$ncontrol.outcome<-344901
dat2$prevalence.outcome<-0.033 
dat2$eaf.outcome <- dat2$eaf.exposure
steiger <- steiger_filtering(dat2)
table(steiger$steiger_dir)

FALSE  TRUE 
   43   1228 

#For depression and CVD
dat3$samplesize.exposure<-480359 
dat3$units.exposure<-"log odds"
dat3$units.outcome <- "log odds"
dat3$samplesize.outcome<-184305
dat3$ncase.exposure<-135458
dat3$ncontrol.exposure<-344901
dat3$ncase.outcome<-60801
dat3$ncontrol.outcome<-123504
dat3$prevalence.exposure<-0.033
dat3$prevalence.outcome<-0.076
steiger <- steiger_filtering(dat3)
table(steiger$steiger_dir)

FALSE  TRUE 
   0   40 

#For EA and Anxiety
dat4$samplesize.exposure<-766343 
dat4$units.exposure<-"SD"
dat4$units.outcome <- "log odds"
dat4$samplesize.outcome<-21806 
dat4$ncase.outcome<-7061
dat4$ncontrol.outcome<-14745
dat4$prevalence.outcome<-0.03
steiger <- steiger_filtering(dat4)
table(steiger$steiger_dir)

FALSE  TRUE 
   505   702 

#For Anxiety and CVD
dat5$samplesize.exposure<-21806 
dat5$units.exposure<-"log odds"
dat5$units.outcome <- "log odds"
dat5$samplesize.outcome<-184305
dat5$ncase.exposure<-7061
dat5$ncontrol.exposure<-14745
dat5$ncase.outcome<-60801
dat5$ncontrol.outcome<-123504
dat5$prevalence.exposure<-0.03
dat5$prevalence.outcome<-0.076
steiger <- steiger_filtering(dat5)
table(steiger$steiger_dir)

FALSE  TRUE 
    0   95 

#For EA and Mental Health
dat6$samplesize.exposure<-766343 
dat6$units.exposure<-"SD"
dat6$units.outcome <- "log odds"
dat6$samplesize.outcome<-322580
dat6$ncase.outcome<-113769
dat6$ncontrol.outcome<-208811
dat6$prevalence.outcome<-0.3 
steiger <- steiger_filtering(dat6)
table(steiger$steiger_dir)

FALSE  TRUE 
   0    595 

#For mental health and CVD
dat7$samplesize.exposure<-322580 
dat7$units.exposure<-"SD"
dat7$units.outcome <- "log odds"
dat7$samplesize.outcome<-184305
dat7$ncase.outcome<-60801
dat7$ncontrol.outcome<-123504
dat7$prevalence.outcome<-0.076 
steiger <- steiger_filtering(dat7)
table(steiger$steiger_dir)

FALSE  TRUE 
   1    12   
   
#For depression and LSI
dat12$samplesize.exposure<-480359 
dat12$units.exposure<-"log odds"
dat12$units.outcome <- "SD"
dat12$samplesize.outcome<-462690
dat12$ncase.exposure<-135458
dat12$ncontrol.exposure<-344901
dat12$prevalence.exposure<-0.033
dat12$eaf.exposure <- dat12$eaf.outcome 
steiger <- steiger_filtering(dat12)
table(steiger$steiger_dir)

FALSE  TRUE 
   0   36  

#For LSI and CVD
dat13$samplesize.exposure<-462690 
dat13$units.exposure<-"SD"
dat13$units.outcome <- "log odds"
dat13$samplesize.outcome<-184305
dat13$ncase.outcome<-60801
dat13$ncontrol.outcome<-123504
dat13$prevalence.outcome<-0.076 
steiger <- steiger_filtering(dat13)
table(steiger$steiger_dir)

FALSE  TRUE 
  10   116  


######################################################################################################################################
#7. Plot results
######################################################################################################################################
    
#For EA and CVD
res1_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res1_loo <- mr_leaveoneout(dat1)
p1 <- mr_scatter_plot(res1, dat1)
ggsave(p1[[1]], file="EA&CVD_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res3_single)
ggsave(p2[[1]], file="EA&CVD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res3_loo)
ggsave(p3[[1]], file="EA&CVD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res3_single)
ggsave(p4[[1]], file="EA&CVD_funnel.png", width=7, height=7)
mr_report(dat1, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "EA&CVD_Report")

#For EA and depression
res2_single <- mr_singlesnp(dat2, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res2_loo <- mr_leaveoneout(dat2)
p1 <- mr_scatter_plot(res2, dat2)
ggsave(p1[[1]], file="EA&MDDscatter.png", width=7, height=7)
p2 <- mr_forest_plot(res2_single)
ggsave(p2[[1]], file="EA&MDD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res2_loo)
ggsave(p3[[1]], file="EA&MDD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res2_single)
ggsave(p4[[1]], file="EA&MDD_funnel.png", width=7, height=7)
mr_report(dat2, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "EA&MDD_Report")
    
#For depression and CVD
res3_single <- mr_singlesnp(dat3, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res3_loo <- mr_leaveoneout(dat3)
p1 <- mr_scatter_plot(res3, dat3)
ggsave(p1[[1]], file="MDD&CVD_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res3_single)
ggsave(p2[[1]], file="MDD&CVD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res3_loo)
ggsave(p3[[1]], file="MDD&CVD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res3_single)
ggsave(p4[[1]], file="MDD&CVD_funnel.png", width=7, height=7)
mr_report(dat3, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "MDD&CVD_Report")
    
#For EA and anxiety
res4_single <- mr_singlesnp(dat4, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res4_loo <- mr_leaveoneout(dat4)
p1 <- mr_scatter_plot(res4, dat4)
ggsave(p1[[1]], file="EA&Anxscatter.png", width=7, height=7)
p2 <- mr_forest_plot(res4_single)
ggsave(p2[[1]], file="EA&Anx_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res4_loo)
ggsave(p3[[1]], file="EA&ANx_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res4_single)
ggsave(p4[[1]], file="EA&Anx_funnel.png", width=7, height=7)
mr_report(dat4, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "EA&Anxiety_Report")
    
#For anxiety and CVD
res5_single <- mr_singlesnp(dat5, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res5_loo <- mr_leaveoneout(dat5)
p1 <- mr_scatter_plot(res5, dat5)
ggsave(p1[[1]], file="Anx&CVD_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res5_single)
ggsave(p2[[1]], file="Anx&CVD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res5_loo)
ggsave(p3[[1]], file="Anx&CVD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res5_single)
ggsave(p4[[1]], file="Anx&CVD_funnel.png", width=7, height=7)
mr_report(dat5, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "Anx&CVD_Report")  
    
#For EA & MH
res6_single <- mr_singlesnp(dat6, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res6_loo <- mr_leaveoneout(dat6)
p1 <- mr_scatter_plot(res6, dat6)
ggsave(p1[[1]], file="EA&MH_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res6_single)
ggsave(p2[[1]], file="EA&MH_single.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res6_loo)
ggsave(p3[[1]], file="EA&MH_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res6_single)
ggsave(p4[[1]], file="EA&MH_funnel.png", width=7, height=7)
mr_report(dat6, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "EA&MH_Report")
    
#For MH & CVD
res7_single <- mr_singlesnp(dat7, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res7_loo <- mr_leaveoneout(dat7)
p1 <- mr_scatter_plot(res7, dat7)
ggsave(p1[[1]], file="MH&CVD_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res7_single)
ggsave(p2[[1]], file="MH&CVD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res7_loo)
ggsave(p3[[1]], file="MH&CVD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res7_single)
ggsave(p4[[1]], file="MH&CVD_funnel.png", width=7, height=7)
mr_report(dat7, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "MH&CVD_Report")        

#For Dep & LSI
res12_single <- mr_singlesnp(dat12, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res12_loo <- mr_leaveoneout(dat12)
p1 <- mr_scatter_plot(res12, dat12)
ggsave(p1[[1]], file="DEP&LSI_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res12_single)
ggsave(p2[[1]], file="DEP&LSI_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res12_loo)
ggsave(p3[[1]], file="DEP&LSI_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res12_single)
ggsave(p4[[1]], file="DEP&LSI_funnel.png", width=7, height=7)
mr_report(dat12, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "DEP&LSI_Report")     

#For LSI & CVD
res13_single <- mr_singlesnp(dat13, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res13_loo <- mr_leaveoneout(dat13)
p1 <- mr_scatter_plot(res13, dat13)
ggsave(p1[[1]], file="LSI&CVD_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res13_single)
ggsave(p2[[1]], file="LSI&CVD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res13_loo)
ggsave(p3[[1]], file="LSI&CVD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res13_single)
ggsave(p4[[1]], file="LSI&CVD_funnel.png", width=7, height=7)
mr_report(dat13, output_path="Report", output_type = "html",
    author = "Daniel Jones", study = "LSI&CVD_Report")  
                            
######################################################################################################################################
#9. Regression Dilution
######################################################################################################################################
#Isq >0.9 = ok to do MR Egger
#0.3<Isq<0.9 = perform SIMEX correction
#Isq<0.3 do not report MR Egger at all

# I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
  }

#calculate Isq wieghted and unweighted in a loop
I2<-c()

###For EA &CVD

#Update ouctome data values and then harmonise data
cvd_outcome_dat.proxies<-1
cvd_outcome_dat.rsq<-0.8
cvd_outcome_dat.align_alleles<-1
cvd_outcome_dat.maf_threshold<-0.3
cvd_outcome_dat.palindromes<-1
dat14 <- harmonise_data(ea_exp_dat, cvd_outcome_dat, action = 1)

#Rename required columns
for(i in 1:length(dat14)){
   	Vars<-dat14[i]
   }
dat14$BetaXG<-dat14$beta.exposure
dat14$seBetaXG<-dat14$se.exposure
dat14$seBetaYG<-dat14$se.outcome
BetaXG   = dat14$BetaXG
seBetaXG = dat14$seBetaXG 
seBetaYG = dat14$seBetaYG
BXG             = abs(BetaXG)

# Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Educational_Attainment", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="EA&CVD_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For EA & MDD
dep_outcome_dat.proxies<-1
dep_outcome_dat.rsq<-0.8
dep_outcome_dat.align_alleles<-1
dep_outcome_dat.maf_threshold<-0.3
dep_outcome_dat.palindromes<-1
dat15 <- harmonise_data(ea_exp_dat, dep_outcome_dat, action = 1)
dat15$eaf.exposure <- dat15$eaf.outcome 

for(i in 1:length(dat15)){
   	Vars<-dat15[i]
   }
dat15$BetaXG<-dat15$beta.exposure
dat15$seBetaXG<-dat15$se.exposure
dat15$seBetaYG<-dat15$se.outcome
BetaXG   = dat15$BetaXG
seBetaXG = dat15$seBetaXG 
seBetaYG = dat15$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Educational_Attainment", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="EA&MDD_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For EA and anxiety
#Update ouctome data values and then harmonise data
anx_outcome_dat.proxies<-1
anx_outcome_dat.rsq<-0.8
anx_outcome_dat.align_alleles<-1
anx_outcome_dat.maf_threshold<-0.3
anx_outcome_dat.palindromes<-1
dat16 <- harmonise_data(ea_exp_dat, anx_outcome_dat, action = 1)

#Rename required columns
for(i in 1:length(dat16)){
   	Vars<-dat16[i]
   }
dat16$BetaXG<-dat16$beta.exposure
dat16$seBetaXG<-dat16$se.exposure
dat16$seBetaYG<-dat16$se.outcome
BetaXG   = dat16$BetaXG
seBetaXG = dat16$seBetaXG 
seBetaYG = dat16$seBetaYG
BXG             = abs(BetaXG)

# Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-c()
I2<-rbind(I2, output)
colnames(I2) <- c("Educational_Attainment", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="EducationalAttainment&Anxiety_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For EA and MH
#Update ouctome data values and then harmonise data
mh_outcome_dat.proxies<-1
mh_outcome_dat.rsq<-0.8
mh_outcome_dat.align_alleles<-1
mh_outcome_dat.maf_threshold<-0.3
mh_outcome_dat.palindromes<-1
dat17 <- harmonise_data(ea_exp_dat, mh_outcome_dat, action = 1)

for(i in 1:length(dat17)){
   	Vars<-dat17[i]
   }
dat17$BetaXG<-dat17$beta.exposure
dat17$seBetaXG<-dat17$se.exposure
dat17$seBetaYG<-dat17$se.outcome
BetaXG   = dat17$BetaXG
seBetaXG = dat17$seBetaXG 
seBetaYG = dat17$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Educational Attainment", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Educational_Attainment&MH_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For MDD and CVD
dat18 <- harmonise_data(dep_exp_dat, cvd_outcome_dat, action = 1)
dat18$eaf.exposure <- dat18$eaf.outcome

for(i in 1:length(dat18)){
   	Vars<-dat18[i]
   }
dat18$BetaXG<-dat18$beta.exposure
dat18$seBetaXG<-dat18$se.exposure
dat18$seBetaYG<-dat18$se.outcome
BetaXG   = dat18$BetaXG
seBetaXG = dat18$seBetaXG 
seBetaYG = dat18$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Depression", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="MDD&CVD_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For anxiety and CVD
dat19 <- harmonise_data(anx_exp_dat, cvd_outcome_dat, action = 1)

for(i in 1:length(dat19)){
   	Vars<-dat19[i]
   }
dat19$BetaXG<-dat19$beta.exposure
dat19$seBetaXG<-dat19$se.exposure
dat19$seBetaYG<-dat19$se.outcome
BetaXG   = dat19$BetaXG
seBetaXG = dat19$seBetaXG 
seBetaYG = dat19$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Anxiety", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Anxiety&CVD_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For Mental Health and CVD
dat20 <- harmonise_data(mh_exp_dat, cvd_outcome_dat, action = 1)

for(i in 1:length(dat20)){
   	Vars<-dat20[i]
   }
dat20$BetaXG<-dat20$beta.exposure
dat20$seBetaXG<-dat20$se.exposure
dat20$seBetaYG<-dat20$se.outcome
BetaXG   = dat20$BetaXG
seBetaXG = dat20$seBetaXG 
seBetaYG = dat20$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Mental Health", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="MH&CVD_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For MDD and LSI
lsi_outcome_dat.proxies<-1
lsi_outcome_dat.rsq<-0.8
lsi_outcome_dat.align_alleles<-1
lsi_outcome_dat.maf_threshold<-0.3
lsi_outcome_dat.palindromes<-1
dat21 <- harmonise_data(dep_exp_dat, lsi_outcome_dat, action = 1)

for(i in 1:length(dat21)){
   	Vars<-dat21[i]
   }
dat21$BetaXG<-dat21$beta.exposure
dat21$seBetaXG<-dat21$se.exposure
dat21$seBetaYG<-dat21$se.outcome
BetaXG   = dat21$BetaXG
seBetaXG = dat21$seBetaXG 
seBetaYG = dat21$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Depression", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="MDD&LSI_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For LSI and CVD
cvd_outcome_dat.proxies<-1
cvd_outcome_dat.rsq<-0.8
cvd_outcome_dat.align_alleles<-1
cvd_outcome_dat.maf_threshold<-0.3
cvd_outcome_dat.palindromes<-1
dat22 <- harmonise_data(lsi_exp_dat, cvd_outcome_dat, action = 1)

for(i in 1:length(dat22)){
   	Vars<-dat22[i]
   }
dat22$BetaXG<-dat22$beta.exposure
dat22$seBetaXG<-dat22$se.exposure
dat22$seBetaYG<-dat22$se.outcome
BetaXG   = dat22$BetaXG
seBetaXG = dat22$seBetaXG 
seBetaYG = dat22$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Lifetime Smoking", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="LSI&CVD_regression_dilution_isq_weighted.csv", row.names = FALSE)

######################################################################################################################################
#10.Simex corrections
######################################################################################################################################

#Install Packages
install.packages("simex")
library(simex)

###For EA and CVD
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat14)){
   	Vars<-dat14[i]
   }
dat14$BetaXG<-dat14$beta.exposure
dat14$BetaYG<-dat14$beta.outcome
dat14$seBetaXG<-dat14$se.exposure
dat14$seBetaYG<-dat14$se.outcome
BetaXG <- dat14$BetaXG
BetaYG <- dat14$BetaYG
seBetaXG <- dat14$seBetaXG
seBetaYG <- dat14$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Educational_Attainment", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Educational_Attainment", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Educational_Attainment&CVD_mreggersimex_weighted.csv", row.names = FALSE)

###For EA & MDD
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat15)){
   	Vars<-dat15[i]
   }
dat15$BetaXG<-dat15$beta.exposure
dat15$BetaYG<-dat15$beta.outcome
dat15$seBetaXG<-dat15$se.exposure
dat15$seBetaYG<-dat15$se.outcome
BetaXG <- dat15$BetaXG
BetaYG <- dat15$BetaYG
seBetaXG <- dat15$seBetaXG
seBetaYG <- dat15$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Educational Attainment", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Depression", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="EA&MDD_mreggersimex_weighted.csv", row.names = FALSE)

###For EA & Anxiety
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat16)){
   	Vars<-dat16[i]
   }
dat16$BetaXG<-dat16$beta.exposure
dat16$BetaYG<-dat16$beta.outcome
dat16$seBetaXG<-dat16$se.exposure
dat16$seBetaYG<-dat16$se.outcome
BetaXG <- dat16$BetaXG
BetaYG <- dat16$BetaYG
seBetaXG <- dat16$seBetaXG
seBetaYG <- dat16$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Educational Attainment", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Educational Attainment", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="EA&Depression_mreggersimex_weighted.csv", row.names = FALSE)

###For EA and Mental Health
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat17)){
   	Vars<-dat17[i]
   }
dat17$BetaXG<-dat17$beta.exposure
dat17$BetaYG<-dat17$beta.outcome
dat17$seBetaXG<-dat17$se.exposure
dat17$seBetaYG<-dat17$se.outcome
BetaXG <- dat17$BetaXG
BetaYG <- dat17$BetaYG
seBetaXG <- dat17$seBetaXG
seBetaYG <- dat17$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Educational Attainment", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Educational Attainment", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Educational_Attainment&MH_mreggersimex_weighted.csv", row.names = FALSE)

###For Anxiety and CVD
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat18)){
   	Vars<-dat18[i]
   }
dat18$BetaXG<-dat18$beta.exposure
dat18$BetaYG<-dat18$beta.outcome
dat18$seBetaXG<-dat18$se.exposure
dat18$seBetaYG<-dat18$se.outcome
BetaXG <- dat18$BetaXG
BetaYG <- dat18$BetaYG
seBetaXG <- dat18$seBetaXG
seBetaYG <- dat18$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 
# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Anxiety", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Anxiety", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Anxiety&CVD_mreggersimex_weighted.csv", row.names = FALSE)

###For MDD & LSI
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat19)){
   	Vars<-dat19[i]
   }
dat19$BetaXG<-dat19$beta.exposure
dat19$BetaYG<-dat19$beta.outcome
dat19$seBetaXG<-dat19$se.exposure
dat19$seBetaYG<-dat19$se.outcome
BetaXG <- dat19$BetaXG
BetaYG <- dat19$BetaYG
seBetaXG <- dat19$seBetaXG
seBetaYG <- dat19$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Depression", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Anxiety", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="LSI&MDD_mreggersimex_weighted.csv", row.names = FALSE)

###For LSI & CVD
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat20)){
   	Vars<-dat20[i]
   }
dat20$BetaXG<-dat20$beta.exposure
dat20$BetaYG<-dat20$beta.outcome
dat20$seBetaXG<-dat20$se.exposure
dat20$seBetaYG<-dat20$se.outcome
BetaXG <- dat20$BetaXG
BetaYG <- dat20$BetaYG
seBetaXG <- dat20$seBetaXG
seBetaYG <- dat20$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Smoking", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("LSI", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="LSI&CVD_mreggersimex_weighted.csv", row.names = FALSE)
