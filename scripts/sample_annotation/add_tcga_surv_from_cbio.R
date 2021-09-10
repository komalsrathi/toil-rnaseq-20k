############################################
# Author: Komal S Rathi
# Function: This script adds surv data 
# from cbio for TCGA patients
# Date: 08/06/2018
# For now don't use it until model is changed
############################################

setwd('~/Projects/toil-rnaseq-20k/')

library(dplyr)
clin <- read.delim('data/clinical/clinical_info.txt', stringsAsFactors = F)
clin$age_normal_tissue_collected_in_years[is.na(clin$age_normal_tissue_collected_in_years)] <- ''

from.cbio <- read.delim('data/clinical/tcga_survival_data.tsv', stringsAsFactors = F)
from.cbio$STABLE_ID_LONG <- NULL
from.cbio$CANCER_STUDY_IDENTIFIER <- NULL
from.cbio <- from.cbio[which(from.cbio$STABLE_ID %in% clin$patient_barcode),]
from.cbio$DFS_STATUS <- gsub('/.*','', from.cbio$DFS_STATUS)
from.cbio[from.cbio == "unavailable"] <- NA

from.cbio <- from.cbio %>% group_by(STABLE_ID) %>% 
  summarise(OS_MONTHS = floor(as.numeric(max(OS_MONTHS))*30), 
            DFS_MONTHS = floor(as.numeric(max(DFS_MONTHS))*30),
            OS_STATUS = toString(unique(na.omit(OS_STATUS))),
            DFS_STATUS = toString(unique(na.omit(DFS_STATUS)))) %>%
  as.data.frame()

# don't include discrepancies for now
from.cbio <- from.cbio[grep(',', from.cbio$OS_STATUS, invert = T),]
from.cbio <- from.cbio[grep(',', from.cbio$DFS_STATUS, invert = T),]
from.cbio <- from.cbio[-which(from.cbio$OS_MONTHS < 0 | from.cbio$DFS_MONTHS < 0),]
from.cbio[is.na(from.cbio)] <- 'unavailable'
from.cbio[from.cbio == ""] <- 'unavailable'

# add to clinical file
# clin <- clin[,c('patient_barcode','vital_status','event_free_status','event_free_survival_time_in_days','overall_survival_time_in_days')]
# clin <- clin[which(clin$patient_barcode %in% from.cbio$STABLE_ID),]
res <- merge(clin, from.cbio, by.x = 'patient_barcode', by.y = 'STABLE_ID', all.x = T)
tc <- colnames(from.cbio[,2:5])
res[which(res$group == "normals"),tc] <- ''
apply(res, 2, FUN = function(x) length(x[is.na(x)])) # what are NAs in the data
res[is.na(res)] <- 'unavailable'

res$vital_status <- ifelse(res$vital_status == "unavailable", res$OS_STATUS, res$vital_status)
res$event_free_survival_time_in_days <- ifelse(res$event_free_survival_time_in_days == "unavailable", res$DFS_MONTHS, res$event_free_survival_time_in_days)
res$overall_survival_time_in_days <- ifelse(res$overall_survival_time_in_days == "unavailable", res$OS_MONTHS, res$overall_survival_time_in_days)
res$event_free_status <- ifelse(res$event_free_status == "unavailable", res$DFS_STATUS, res$event_free_status)
res$OS_MONTHS <- NULL
res$DFS_MONTHS <- NULL
res$OS_STATUS <- NULL
res$DFS_STATUS <- NULL
write.table(res, file = '~/Projects/toil-rnaseq-20k/data/clinical/clinical_info.txt', quote = F, sep = "\t", row.names = F)
