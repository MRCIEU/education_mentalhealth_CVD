#Script for analysing the Lifetime smoking summary data in two-sample MR bi-directionally
#UPDATED for the newest version of MRBase R package
#For more information on commands/additional options see: https://mrcieu.github.io/TwoSampleMR/
#Created by Robyn Wootton September 2018 - adapted for MVMR by Hannah Sallis June 2019
#Updated by Daniel Jones August 2019
#######################################################################################################################


######################################################################################################################################
#1. Load packages
######################################################################################################################################
rm(list=ls(all=TRUE)) #empties your R environment

#Load packages - you will have to install the packages first time you run the script
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)
library(knitr)

######################################################################################################################################
#2. Read exposure data
######################################################################################################################################
#The following code is for reading in your own data
#If instead, you want to use data from MR Base, see https://mrcieu.github.io/TwoSampleMR/

# Make sure headings are formatted like so: 
# "Phenotype 	SNP			CHR 	BP 			effect_allele 	other_allele 		eaf   beta 	  se 		  pval"
# "Neuroticism 	rs4653651 	1 		225862060 	A 				G 				0.32	-0.091 	0.0202 	6.443e-06"

# THE EXPOSURE DATASETS MUST CONTAIN THE SNPS FROM BOTH INSTRUMENTS

###Read in initial lead SNPs
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

#Lifetime Smoking
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

###Create additional instruments for each exposure that pick out the lead SNPs from the relevant other instruments
#E.g. For using EA with mental health
mh_exp_dat <- read_exposure_data(
       filename = "DepressionLeadSNPs.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "Beta",
       se_col = "SE",
       effect_allele_col = "A1",
       other_allele_col = "A2",
       eaf_col = "EAF",
       pval_col = "Pval",
   )
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
ea_outcome_dat <- ea_outcome_dat[ea_outcome_dat$SNP %in% mh_exp_dat$SNP,] #This picks out the lead SNPs from your other exposure data
write.csv(ea_outcome_dat, "EA(MDD)-MVMR.csv", quote=F, row.names=F)

##Read in these additional instruments 
#For EA and Mental Health  
ea2_exp_dat <- read_exposure_data("EA(MH)-MVMR.csv", sep=",")
ea2_exp_dat$exposure <- "Educational Attainment"
  
mh2_exp_dat <- read_exposure_data("MH-MVMR.csv", sep=",")
mh2_exp_dat$exposure <- "Mental Health"

#Read in new instruments for EA and Anxiety
ea3_exp_dat <- read_exposure_data("EA(Anx)-MVMR.csv", sep=",")
ea3_exp_dat$exposure <- "Educational Attainment"
  
anx2_exp_dat <- read_exposure_data("Anxiety-MVMR.csv", sep=",")
anx2_exp_dat$exposure <- "Anxiety"

#Read in new instruments for EA and MDD
ea4_exp_dat <- read_exposure_data("EA(MDD)-MVMR.csv", sep=",")
ea4_exp_dat$exposure <- "Educational Attainment"

dep2_exp_dat <- read_exposure_data("MDD-MVMR.csv", sep=",")
dep2_exp_dat$exposure <- "Depression"

#Read in new instruments for LSI and Depression analysis
lsi2_exp_dat <- read_exposure_data("LSI(MDD)-MVMR.csv", sep=",")
lsi2_exp_dat$exposure <- "Lifetime Smoking"

lsi3_exp_dat <- read_exposure_data("LSI(EA)-MVMR.csv", sep=",")
lsi3_exp_dat$exposure <- "Lifetime Smoking"

dep2_exp_dat <- read_exposure_data("MDD-MVMR.csv", sep=",")
dep2_exp_dat$exposure <- "Depression"

dep3_exp_dat <- read_exposure_data("MDD(LSI)-MVMR.csv", sep=",")
dep3_exp_dat$exposure <- "Depression"

ea4_exp_dat <- read_exposure_data("EA(MDD)-MVMR.csv", sep=",")
ea4_exp_dat$exposure <- "Educational Attainment"

ea5_exp_dat <- read_exposure_data("EA(LSI)-MVMR.csv", sep=",")
ea5_exp_dat$exposure <- "Educational Attainment"

##Append the dataframes
#For EA and mental health
ea2_exp_dat <- subset(ea2_exp_dat, !(ea2_exp_dat$SNP %in% ea_exp_dat$SNP))
ea2_exp_dat<-rbind(ea_exp_dat,ea2_exp_dat)
ea2_exp_dat$id.exposure<-"lHOAgH" #There must be only one ID per exposure for harmonisation
mh2_exp_dat <- subset(mh2_exp_dat, !(mh2_exp_dat$SNP %in% mh_exp_dat$SNP))
mh2_exp_dat<-rbind(mh_exp_dat,mh2_exp_dat)
mh2_exp_dat$id.exposure<-"0LMMTd"
ea2_exp_dat <- ea2_exp_dat[ea2_exp_dat$SNP %in% mh2_exp_dat$SNP,]
mh2_exp_dat <- mh2_exp_dat[mh2_exp_dat$SNP %in% ea2_exp_dat$SNP,]
exposure_dat <- rbind(ea2_exp_dat,mh2_exp_dat)

#For EA and anxiety
ea3_exp_dat <- subset(ea3_exp_dat, !(ea3_exp_dat$SNP %in% ea_exp_dat$SNP))
ea3_exp_dat<-rbind(ea_exp_dat,ea3_exp_dat)
ea3_exp_dat$id.exposure<-"lHOAgH" #There must be only one ID per exposure for harmonisation
anx2_exp_dat <- subset(anx2_exp_dat, !(anx2_exp_dat$SNP %in% anx_exp_dat$SNP))
anx2_exp_dat<-rbind(anx_exp_dat,anx2_exp_dat)
anx2_exp_dat$id.exposure<-"8YmxF0"
ea3_exp_dat <- ea3_exp_dat[ea3_exp_dat$SNP %in% anx2_exp_dat$SNP,]
anx2_exp_dat <- anx2_exp_dat[anx2_exp_dat$SNP %in% ea3_exp_dat$SNP,]
exposure2_dat <- rbind(ea3_exp_dat,anx2_exp_dat)

#For EA and depression
ea4_exp_dat <- subset(ea4_exp_dat, !(ea4_exp_dat$SNP %in% ea_exp_dat$SNP))
ea4_exp_dat<-rbind(ea_exp_dat,ea4_exp_dat)
ea4_exp_dat$id.exposure<-"lHOAgH" #There must be only one ID per exposure for harmonisation
dep2_exp_dat <- subset(dep2_exp_dat, !(dep2_exp_dat$SNP %in% dep_exp_dat$SNP))
dep2_exp_dat<-rbind(dep_exp_dat,dep2_exp_dat)
dep2_exp_dat$id.exposure<-"Vz2RFT"
ea4_exp_dat <- ea4_exp_dat[ea4_exp_dat$SNP %in% dep2_exp_dat$SNP,]
dep2_exp_dat <- dep2_exp_dat[dep2_exp_dat$SNP %in% ea4_exp_dat$SNP,]
dep2_exp_dat$eaf.exposure<-ea4_exp_dat$eaf.exposure
exposure3_dat <- rbind(ea4_exp_dat,dep2_exp_dat)

#For EA & Depression on Smoking
ea4_exp_dat <- subset(ea4_exp_dat, !(ea4_exp_dat$SNP %in% ea_exp_dat$SNP))
ea4_exp_dat<-rbind(ea_exp_dat,ea4_exp_dat)
ea4_exp_dat$id.exposure<-"lHOAgH" #There must be only one ID per exposure for harmonisation
dep2_exp_dat <- subset(dep2_exp_dat, !(dep2_exp_dat$SNP %in% dep_exp_dat$SNP))
dep2_exp_dat<-rbind(dep_exp_dat,dep2_exp_dat)
dep2_exp_dat$id.exposure<-"Vz2RFT"
ea4_exp_dat <- ea4_exp_dat[ea4_exp_dat$SNP %in% dep2_exp_dat$SNP,]
dep2_exp_dat <- dep2_exp_dat[dep2_exp_dat$SNP %in% ea4_exp_dat$SNP,]
exposure4_dat <- rbind(dep2_exp_dat, ea4_exp_dat)

#For EA, Depression and Smoking on CVD
dep2_exp_dat <- subset(dep2_exp_dat, !(dep2_exp_dat$SNP %in% dep_exp_dat$SNP))
dep2_exp_dat<-rbind(dep_exp_dat,dep2_exp_dat)
dep3_exp_dat <- subset(dep3_exp_dat, !(dep3_exp_dat$SNP %in% dep2_exp_dat$SNP))
dep3_exp_dat<-rbind(dep3_exp_dat,dep2_exp_dat)
dep3_exp_dat$id.exposure<-"Vz2RFT"
lsi2_exp_dat <- subset(lsi2_exp_dat, !(lsi2_exp_dat$SNP %in% lsi_exp_dat$SNP))
lsi2_exp_dat<-rbind(lsi_exp_dat,lsi2_exp_dat)
lsi3_exp_dat <- subset(lsi3_exp_dat, !(lsi3_exp_dat$SNP %in% lsi2_exp_dat$SNP))
lsi3_exp_dat<-rbind(lsi2_exp_dat,lsi3_exp_dat)
lsi3_exp_dat$id.exposure<-"6rJOAy"
dep3_exp_dat <- dep3_exp_dat[dep3_exp_dat$SNP %in% lsi3_exp_dat$SNP,]
dep3_exp_dat$eaf.exposure<-lsi3_exp_dat$eaf.exposure
exposure5_dat <- rbind(lsi3_exp_dat, dep3_exp_dat)

######################################################################################################################################
#3. Read in outcome data
######################################################################################################################################

#These must be formatted by restricting SNPs to SNPs in the relevant instruments
##For EA & Mental Health on CVD
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
cvd_outcome_dat <- cvd_outcome_dat[cvd_outcome_dat$SNP %in% exposure_dat$SNP,]

#restrict exposure dataset to snps in outcome
exp_subset <- exposure_dat[exposure_dat$SNP %in% cvd_outcome_dat$SNP,]

##For EA and anxiety on CVD
cvd2_outcome_dat <- read_outcome_data(
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
cvd2_outcome_dat $outcome<-"CVD"
cvd2_outcome_dat <- cvd2_outcome_dat[cvd2_outcome_dat$SNP %in% exposure2_dat$SNP,]

#restrict exposure dataset to snps in outcome
exp2_subset <- exposure2_dat[exposure2_dat$SNP %in% cvd2_outcome_dat$SNP,]

##For EA and depression on CVD
cvd3_outcome_dat <- read_outcome_data(
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
cvd3_outcome_dat $outcome<-"CVD"
cvd3_outcome_dat <- cvd3_outcome_dat[cvd3_outcome_dat$SNP %in% exposure3_dat$SNP,]

#restrict exposure dataset to snps in outcome
exp3_subset <- exposure3_dat[exposure3_dat$SNP %in% cvd3_outcome_dat$SNP,]


##For EA & MDD on Smoking
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
lsi_outcome_dat <- lsi_outcome_dat[lsi_outcome_dat$SNP %in% exposure4_dat$SNP,]

#restrict exposure dataset to snps in outcome
exp4_subset <- exposure4_dat[exposure4_dat$SNP %in% lsi_outcome_dat$SNP,]

##For EA & MDD & smoking on cvd
cvd5_outcome_dat <- read_outcome_data(
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
cvd5_outcome_dat $outcome<-"CVD"
cvd5_outcome_dat <- cvd5_outcome_dat[cvd5_outcome_dat$SNP %in% exposure5_dat$SNP,]

#restrict exposure dataset to snps in outcome
exp5_subset <- exposure5_dat[exposure5_dat$SNP %in% cvd5_outcome_dat$SNP,]


######################################################################################################################################
#4. Harmonise data
######################################################################################################################################
#For EA and MH
mvdat1 <- mv_harmonise_data(exp_subset, cvd_outcome_dat)

#For EA and anxiety
mvdat2 <- mv_harmonise_data(exp2_subset, cvd2_outcome_dat)

#For EA and depression
mvdat3 <- mv_harmonise_data(exp3_subset, cvd3_outcome_dat)

#For EA & depression ON LSI
mvdat4 <- mv_harmonise_data(exp4_subset, lsi_outcome_dat)

#For depression and LSI ON CVD
mvdat5 <- mv_harmonise_data(exp5_subset, cvd5_outcome_dat)

######################################################################################################################################
#5. Run MVMR
######################################################################################################################################
#res <- mv_multiple(mvdat)
#res <- mv_multiple(mvdat,instrument_specific = TRUE)
#res <- mv_multiple(mvdat, pval_threshold = 1)

res1 <- mv_residual(mvdat1, pval_threshold = 1, instrument_specific = TRUE, plots = TRUE)
res1$result
#exposure 		 	outcome nsnp    b          se         pval
#1 EA    		 	  CVD	 573 -0.4595952  0.04843343 2.327788e-21
#2 MH			      CVD	 573  0.3950711. 0.24378664 1.051118e-01

Tab1 <- cbind(mvdat1[["exposure_beta"]], mvdat1[["exposure_se"]], mvdat1[["exposure_pval"]],mvdat1[["outcome_beta"]], mvdat1[["outcome_se"]], mvdat1[["outcome_pval"]])
write.table(Tab1,file="EA-MH-CVD-mvmr_betas.csv",sep=",")

res2 <- mv_residual(mvdat2, pval_threshold = 1, instrument_specific = TRUE, plots = TRUE)
res2$result
#exposure 		 	outcome nsnp    b           se         pval
#1 EA    		 	  CVD	 1146 -0.4764533371 0.032882663 1.411673e-47
#2 Anxiety		      CVD	 1146 -0.0003637509 0.009701754 9.700917e-01

Tab2 <- cbind(mvdat2[["exposure_beta"]], mvdat2[["exposure_se"]], mvdat2[["exposure_pval"]],mvdat2[["outcome_beta"]], mvdat2[["outcome_se"]], mvdat2[["outcome_pval"]])
write.table(Tab2,file="EA-Anxiety-mvmr_betas.csv",sep=",")

res3 <- mv_residual(mvdat3, pval_threshold = 1, instrument_specific = TRUE, plots = TRUE)
res3$result
#exposure 		 	outcome nsnp    b          se         pval
#1 EA    		 	  CVD	 1206 -0.36206739 0.03368275 5.968537e-27
#2 Depression	      CVD	 1206  0.09452519 0.02930712 1.258268e-03

Tab3 <- cbind(mvdat3[["exposure_beta"]], mvdat3[["exposure_se"]], mvdat3[["exposure_pval"]],mvdat3[["outcome_beta"]], mvdat3[["outcome_se"]], mvdat3[["outcome_pval"]])
write.table(Tab3,file="EA-MDD-mvmr_betas.csv",sep=",")

res4 <- mv_residual(mvdat4, pval_threshold = 1, instrument_specific = TRUE, plots = TRUE)
res4$result
#exposure 		 	outcome nsnp    b          se         pval
#1 EA    		 	  LSI	1087 -0.20290128 0.007405947  2.971750e-165
#2 MDD			      LSI	1087  0.03565828 0.005707266  4.160796e-10

Tab4 <- cbind(mvdat4[["exposure_beta"]], mvdat4[["exposure_se"]], mvdat4[["exposure_pval"]],mvdat4[["outcome_beta"]], mvdat4[["outcome_se"]], mvdat4[["outcome_pval"]])
write.table(Tab4,file="MDD-EA-LSI-mvmr_betas.csv",sep=",")

res5 <- mv_residual(mvdat5, pval_threshold = 1, instrument_specific = TRUE, plots = TRUE)
res5$result
#exposure 		 	outcome nsnp    b          se         pval
#1 MDD    		 	  CVD	 160  0.07483449 0.02861218 0.0089103078
#2 LSI			      CVD	 160  0.37619924 0.09845675 0.0001329336

Tab5 <- cbind(mvdat5[["exposure_beta"]], mvdat4[["exposure_se"]], mvdat4[["exposure_pval"]],mvdat4[["outcome_beta"]], mvdat4[["outcome_se"]], mvdat4[["outcome_pval"]])
write.table(Tab5,file="LSI-MDD-mvmr_betas.csv",sep=",")