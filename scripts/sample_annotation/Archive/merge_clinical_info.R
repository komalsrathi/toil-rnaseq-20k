###################################
# Author: Komal S Rathi
# Date: 08/03/2018
# Function: Merge all clinical data 
###################################

library(reshape2)

# merge with target data (Left from here....)
load('~/Projects/toil-rnaseq-20k/data/clinical/targetclin.RData')
load('~/Projects/toil-rnaseq-20k/data/clinical/tcgaclin.RData')
total <- rbind(target.clin, tcga.clin.sub)
sample.info <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/sample_info.txt', stringsAsFactors = F)
clin <- merge(sample.info, total, by.x = c("analysis_id","disease_subtype"), by.y = c("sample_id","disease_subtype"), all.x = T)

# now fill in missing info
# what are all the columns common for all studies
cc <- c('analysis_id',
        'patient_barcode',
        'sample_barcode',
        'group',
        'study',
        'definition',
        'library_type',
        'platform',
        'center',
        'gender',
        'race',
        'ethnicity')
tc <- c('disease','disease_name',
        'disease_subtype',
        'age_at_diagnosis_in_days',
        'event_free_survival_time_in_days',
        'vital_status',
        'overall_survival_time_in_days',
        'tags')

# common columns
fill.cc <- function(x){
  x[is.na(x)] <- 'unavailable'
  return(x)
}

clin[,cc] <- apply(clin[,cc], 2, fill.cc)
clin[which(clin$group == "Tumors"),tc] <- apply(clin[which(clin$group == "Tumors"),tc], 2, fill.cc)
clin[which(clin$group == "Normals"),tc] <- ''
clin$vital_status[clin$group == "Normals"] <- 'dead'

# prepare GTEx
gtex.clin <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/gtex/gtex_clinical.txt', stringsAsFactors = F)
gtex.clin$Age <- NULL

# now merge GTEx with clinical 
clin <- merge(clin, gtex.clin, by = 'analysis_id', all.x = TRUE)
clin$age_normal_tissue_collected_in_years[clin$group == "Tumors"] <- ''
clin$gender <- ifelse(clin$group == "Normals", clin$Gender, clin$gender)
clin$Gender <- NULL

# reorder columns but you can skip this step
clin <- clin[,c(1:16,22,17,18,20,19,21)]
colnames(clin)[1] <- 'sample_id'
colnames(clin)[6] <- 'study_id'

# prepare CBTTC
cbttc <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/cbttc/cbttc_metadata.txt', stringsAsFactors = F)
cbttc[is.na(cbttc)] <- ''

# merge CBTTC with clin
clin <- rbind(clin, cbttc)

# merge PNOC with clin
pnoc <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/pnoc/pnoc_data_clinic_info.txt', stringsAsFactors = F)
pnoc[is.na(pnoc)] <- ''
clin <- rbind(clin, pnoc)

# convert some values to lower case
clin$group <- tolower(clin$group)

# add event_free_status
# clin$event_free_status <- ifelse(clin$group == "normals", '', 'unavailable')

# write out for diseasexpress
cols <- c("analysis_id","patient_barcode","sample_barcode","group","study","disease","disease_name","disease_subtype","tissue","definition","library_type","platform","center","gender","race","ethnicity","age_normal_tissue_collected_in_years","age_at_diagnosis_in_days","event_free_survival_time_in_days","overall_survival_time_in_days","vital_status","tags")
colnames(clin)[1] <- 'analysis_id'
colnames(clin)[6] <- 'study'
clin <- clin[,cols]
clin$tags[clin$tags == "unavailable"] <- ""
write.table(clin, file = '~/Projects/toil-rnaseq-20k/data/clinical/clinical_info.txt', quote = F, sep = "\t", row.names = F)

