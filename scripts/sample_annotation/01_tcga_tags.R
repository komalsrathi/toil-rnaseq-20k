###################################
# Author: Komal S Rathi
# Date: 08/03/2018
# Function: Create tags for TCGA 
# Note: event free survival later**
###################################

library(reshape2)
library(plyr)
library(dplyr)

# load target clinical data
load('~/Projects/toil-rnaseq-20k/data/clinical/targetclin.RData')
sample.info <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/sample_info.txt', stringsAsFactors = F)

# read list of clinical files from GDC
l <- list.files(path = '~/Projects/toil-rnaseq-20k/data/metadata/tcga/clinical/gdc_download_20161214_200436/', pattern = ".txt", recursive = T, full.names = T)
l <- l[grep('patient',l)]

tcga.clin <- data.frame()
for(i in 1:length(l)){
  x <- read.delim(file = l[i], skip = 1, check.names = F, stringsAsFactors = F)
  x <- x[-1,]
  if(i == 1){
    tcga.clin <- x
  } else {
    tcga.clin <- plyr::rbind.fill(tcga.clin, x) 
  }
}
tcga.clin <- tcga.clin[which(tcga.clin$bcr_patient_barcode %in% sample.info$patient_barcode),]
cols <- c('pathologic_stage','tumor_type','igcccg_stage','hypertension','diabetes','hpv_status_by_p16_testing','hpv_status_by_ish_testing')
cols.standard <- c('bcr_patient_barcode','gender','race','ethnicity',
                   'age_at_initial_pathologic_diagnosis',
                   'vital_status','days_to_death')
tcga.clin <- tcga.clin[,c(cols.standard, cols)]

# fix days
tcga.clin$age_at_initial_pathologic_diagnosis <- as.numeric(as.character(tcga.clin$age_at_initial_pathologic_diagnosis))
tcga.clin$age_at_initial_pathologic_diagnosis <- floor(tcga.clin$age_at_initial_pathologic_diagnosis*365)
tcga.clin$days_to_death <- as.numeric(as.character(tcga.clin$days_to_death))
tcga.clin$days_to_death <- ifelse(is.na(tcga.clin$days_to_death), NA, tcga.clin$age_at_initial_pathologic_diagnosis + tcga.clin$days_to_death)

# only select the below as the others are mostly NAs
# clean up
tcga.clin$pathologic_stage <- ifelse(tcga.clin$pathologic_stage %in% c("I/II NOS","IS"), paste0('Stage ', tcga.clin$pathologic_stage), tcga.clin$pathologic_stage)
clean.up <- function(x){
  x[grep("[[]",x)] <- 'unavailable'
  x[is.na(x)] <- 'unavailable'
  if(length(grep('Stage|TCGA', x)) == 0){
    x <- tolower(x)
  }
  return(x)
}
tcga.clean <- as.data.frame(apply(tcga.clin, 2, clean.up), stringsAsFactors = F)
tcga.clean$hpv_status <- tcga.clean$hpv_status_by_p16_testing
tcga.clean$hpv_status <- ifelse(tcga.clean$hpv_status == "unavailable", tcga.clean$hpv_status_by_ish_testing, tcga.clean$hpv_status)
tcga.clean$hpv_status_by_ish_testing <- NULL
tcga.clean$hpv_status_by_p16_testing <- NULL
tcga.clean$pathologic_stage <- gsub('Stage ','', tcga.clean$pathologic_stage)
tcga.clean$pathologic_stage <- sub(' .*','',tcga.clean$pathologic_stage)
tcga.clean$tumor_type <- sub(' ','',tolower(tcga.clean$tumor_type))

# create tags column
tcga.clean <- merge(sample.info[,c('analysis_id','patient_barcode','disease_subtype')], 
                    tcga.clean, by.x = 'patient_barcode', by.y = 'bcr_patient_barcode')
tags <- melt(tcga.clean, measure.vars = c('pathologic_stage','tumor_type','igcccg_stage','hypertension','diabetes','hpv_status'),
          na.rm = TRUE, variable.name = 'tags')
tags$tags <- paste0(tags$tags,'=',tags$value)
tags$value <- NULL

tcga.clin.sub <- tags %>% 
  group_by(patient_barcode, analysis_id, disease_subtype, gender, 
           race, ethnicity, age_at_initial_pathologic_diagnosis, 
           vital_status, days_to_death) %>%
  summarise(tags = toString(tags)) %>% as.data.frame()
tcga.clin.sub$tags <- gsub(', ','; ',tcga.clin.sub$tags)
tcga.clin.sub$patient_barcode <- NULL
colnames(tcga.clin.sub)[c(1,6,8)] <- c('sample_id','age_at_diagnosis_in_days','overall_survival_time_in_days')
tcga.clin.sub$event_free_survival_time_in_days <- 'unavailable'
tcga.clin.sub <- tcga.clin.sub[,colnames(target.clin)] 

# save data
save(tcga.clin.sub, file = '~/Projects/toil-rnaseq-20k/data/clinical/tcgaclin.RData')
