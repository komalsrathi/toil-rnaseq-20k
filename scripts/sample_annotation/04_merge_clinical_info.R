###################################
# Author: Komal S Rathi
# Date: 08/20/2018
# Function: Merge all clinical data
###################################

library(reshape2)

# merge with target data (Left from here....)
load('~/Projects/toil-rnaseq-20k/data/clinical/targetclin.RData')
load('~/Projects/toil-rnaseq-20k/data/clinical/tcgaclin.RData')
load('~/Projects/toil-rnaseq-20k/data/clinical/targetclin_additional.RData')
target.clin <- rbind(target.clin, target.clin2)
rm(target.clin2)

total <- rbind(target.clin, tcga.clin.sub)
sample.info <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/sample_info.txt', stringsAsFactors = F)
clin <- merge(sample.info, total, by.x = 'analysis_id', by.y = 'sample_id', all.x = TRUE)
clin$disease_subtype.x <- ifelse(clin$disease_subtype.x == "", clin$disease_subtype.y, clin$disease_subtype.x)
clin$disease_subtype.y <- NULL
colnames(clin)[colnames(clin) == "disease_subtype.x"] <- "disease_subtype"
clin$disease_subtype[clin$disease == "SCLC"] <- "Small Cell Lung Cancer"
# clin <- merge(sample.info, total, by.x = c("analysis_id","disease_subtype"), by.y = c("sample_id","disease_subtype"), all.x = T)

# now fill in missing info
# what are all the columns common for all studies
cc <- c('analysis_id',
        'patient_barcode',
        'sample_barcode',
        'group',
        'study',
        'definition',
        'gender',
        'race',
        'ethnicity')
tc <- c('disease',
        'disease_name',
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
# clin <- clin[,c(1:16,22,17,18,20,19,21)]
# colnames(clin)[colnames(clin) == "analysis_id"] <- 'sample_id'
# colnames(clin)[colnames(clin) == "study"] <- 'study_id'

# prepare CBTTC
c1 <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/cbttc/cbttc_metadata.txt', stringsAsFactors = F)
colnames(c1)[colnames(c1) == "sample_id"] <- 'analysis_id'
colnames(c1)[colnames(c1) == "study_id"] <- 'study'
c2 <- read.csv('~/Projects/toil-rnaseq-20k/data/metadata/cbttc/cbttc_metadata_v2.csv', stringsAsFactors = F)
c3 <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/cbttc/cbttc_metadata_v3.txt', stringsAsFactors = F)
cbttc <- rbind(c1, c2, c3)
cbttc[is.na(cbttc)] <- ''
cbttc$subtissue <- ''
cbttc <- cbttc[,colnames(clin)]

# merge CBTTC with clin
clin <- rbind(clin, cbttc)

# prepare PNOC
pnoc <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/pnoc/pnoc_data_clinic_info.txt', stringsAsFactors = F)
pnoc[is.na(pnoc)] <- ''
pnoc$subtissue <- ''
colnames(pnoc)[1] <- 'analysis_id'
colnames(pnoc)[5] <- 'study'
pnoc <- pnoc[,colnames(clin)]

# merge with PNOC
clin <- rbind(clin, pnoc)

# convert some values to lower case
clin$group <- tolower(clin$group)
clin$gender <- tolower(clin$gender)
clin$race <- tolower(clin$race)
clin$ethnicity <- tolower(clin$ethnicity)
clin$vital_status <- tolower(clin$vital_status)

# remove unwanted 
clin <- clin[which(clin$study != "MB_TaylorLab"),]
clin <- clin[which(clin$definition != "Control Analyte"),]
dups <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/target_dups.txt')
clin <- clin[-which(clin$analysis_id %in% dups$x),]
clin <- unique(clin)

# write out for diseasexpress
cols <- c("analysis_id","patient_barcode","sample_barcode","group","study","disease","disease_name","disease_subtype","tissue","subtissue","definition","gender","race","ethnicity","age_normal_tissue_collected_in_years","age_at_diagnosis_in_days","event_free_survival_time_in_days","overall_survival_time_in_days","vital_status","tags")
clin <- clin[,cols]
clin$tags[clin$tags == "unavailable"] <- ""
clin[is.na(clin)] <- ""
write.table(clin, file = '~/Projects/toil-rnaseq-20k/data/clinical/final_clinical/clinical_info.txt', quote = F, sep = "\t", row.names = F)

