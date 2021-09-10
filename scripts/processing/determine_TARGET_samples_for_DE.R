###################################
# Author: Komal S Rathi
# Date: 08/03/2018
# Function: Determine which samples do not have clinical data from 
# additional samples and remove them before adding to DE
###################################

library(reshape2)
library(plyr)
library(dplyr)

setwd('~/Projects/toil-rnaseq-20k/data/metadata/target/clinical')
sample.info <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/target/to_process/data/TARGET_samples_toadd.txt', stringsAsFactors = F)
sample.info$patient_barcode <- gsub('-[0-9]{2}[A-Z]{1}-[0-9]{2}[A-Z]{1}','',sample.info$barcode)

########################################################################
# TARGET
# ALL (n = 149 patients)
all1 <- read.csv('ALL_Phase1/TARGET_ALL_Phase_I_ClinicalData_20170525.csv', check.names = F, stringsAsFactors = F)
all2 <- read.csv('ALL_Phase2/TARGET_ALL_Dicentric_ClinicalData_20170525.csv', check.names = F, stringsAsFactors = F)
all2.1 <- read.csv('ALL_Phase2/TARGET_ALL_Phase_II_ClinicalData_20170525.csv', check.names = F, stringsAsFactors = F)
all2.2 <- read.csv('ALL_Phase2/TARGET_ALL_Xenografts_ClinicalData_20170525.csv', check.names = F, stringsAsFactors = F)
all2.3 <- read.csv('ALL_Phase2/TARGET_ALL_Validation_ClinicalData_20170525.csv', check.names = F, stringsAsFactors = F)
all1 <- all1[,colnames(all2.1)]
all2 <- all2[,colnames(all2.1)]
all2.2 <- all2.2[,colnames(all2.1)]
all2.3 <- all2.3[,colnames(all2.1)]
all <- unique(rbind(all1, all2, all2.1, all2.2, all2.3))
all <- unique(all$`TARGET USI`)
all <- all[all %in% sample.info$patient_barcode] # 338

# AML (n = 186 patients)
aml1 <- read.csv('AML/TARGET_AML_Discovery_ClinicalData_20170525.csv',check.names = F, stringsAsFactors = F)
aml2 <- read.csv('AML/TARGET_AML_Validation_ClinicalData_20170525.csv', check.names = F, stringsAsFactors = F)
aml <- rbind(aml1, aml2)
aml <- unique(aml$`TARGET USI`)
aml <- aml[aml %in% sample.info$patient_barcode] # 151

# NBL 
nbl <- read.csv('NBL/TARGET_NBL_ClinicalData_20151124.csv', check.names = F, stringsAsFactors = F)
nbl <- unique(nbl$`TARGET USI`)
nbl <- nbl[nbl %in% sample.info$patient_barcode] # 7

# OS 
os.1 <- read.csv('OS/TARGET_OS_ClinicalData_Discovery_20150729.csv', check.names = F, stringsAsFactors = F)
os.2 <- read.csv('OS/TARGET_OS_ClinicalData_Validation_20160401.csv', check.names = F, stringsAsFactors = F)
os <- unique(rbind(os.1, os.2))
os <- unique(os$`TARGET USI`)
os <- os[os %in% sample.info$patient_barcode] # 87

# RT
rt.1 <- read.csv('RT/TARGET_RT_ClinicalData_Discovery_20150710_public.csv', stringsAsFactors = F, check.names = F)
rt.2 <- read.csv('RT/TARGET_RT_ClinicalData_Validation_20150710_public.csv', stringsAsFactors = F, check.names = F)
rt <- unique(rbind(rt.1, rt.2))
rt <- unique(rt$`TARGET USI`)
rt <- rt[rt %in% sample.info$patient_barcode] # 65

# WT
wt.1 <- read.csv('WT/TARGET_WT_ClinicalData_Discovery_20160714_public.csv', stringsAsFactors = F, check.names = F)
wt.2 <- read.csv('WT/TARGET_WT_ClinicalData_Validation_20160714_public.csv', stringsAsFactors = F, check.names = F)
wt <- unique(rbind(wt.1, wt.2))
wt <- unique(wt$`TARGET USI`)
wt <- wt[wt %in% sample.info$patient_barcode] # 4

# CCSK
ccsk <- read.csv('CCSK/TARGET_CCSK_Discovery_ClinicalData_20170525.csv', stringsAsFactors = F, check.names = F)
ccsk <- unique(ccsk$`TARGET USI`)
ccsk <- ccsk[ccsk %in% sample.info$patient_barcode] # 13

# merge
total <- unique(c(all, aml, ccsk, nbl, os, rt, wt))
neg <- sample.info[-which(sample.info$patient_barcode %in% total),]

sample.info <- sample.info[-which(sample.info$analysis_id %in% neg$analysis_id),]
write.table(sample.info, file = '~/Projects/toil-rnaseq-20k/data/metadata/target/to_process/data/TARGET_samples_toadd_v2.txt', quote = F, sep = "\t", row.names = F)
