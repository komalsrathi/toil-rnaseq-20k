#######################################
# Author: Komal S Rathi
# Function: Creates the input table for 
# coxReg.R
#######################################

setwd('~/Projects/toil-rnaseq-20k/data/metadata/target/TARGET_clinical_data/')
dat <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata_filtered/filtered_datasets.txt', stringsAsFactors = F)
dat$USI <- gsub('-[0-9]{2}[A-Z]{1}-[0-9]{2}[A-Z]{1}$|-[0-9]{2}[A-Z]{1}-[0-9]{2}[A-Z]{1}-[A-Z0-9]{4}-[0-9]{2}$','',dat$barcode)
disease <- unique(dat[which(dat$study == 'TARGET'),'disease'])
cols.to.keep <- c('TARGET.USI','Vital.Status','Overall.Survival.Time.in.Days')
files <- list.files(path = '.', pattern = '.csv', recursive = T)
y <- grep(paste(disease, collapse = "|"), files)

target <- data.frame()
for(i in 1:length(files)){
  if(i %in% y){
    x <- read.csv(files[i], stringsAsFactors = F)
    x <- x[,cols.to.keep]
    if(i == 1){
      target <- x
    }
    if(i > 1){
      target <- rbind(target, x)
    }
  }
}

target$Vital.Status <- sub(' .*','',tolower(target$Vital.Status))
target <- target[which(target$Vital.Status == "alive" | target$Vital.Status == "dead"),]
colnames(target) <- c('USI','eventVar', 'TimeVar')
target$eventVar <- ifelse(target$eventVar == "dead", 1, 0)

# TCGA
l <- list.files(path = '~/Projects/toil-rnaseq-20k/data/metadata/tcga/clinical/gdc_download_20161214_200436/', pattern = ".txt", recursive = T, full.names = T)
l <- l[grep('patient',l)]

tcga <- data.frame()
for(i in 1:length(l)){
  x <- read.delim(file = l[i], skip = 1)
  x <- x[-1,]
  x <- x[,c("bcr_patient_barcode", "days_to_last_followup", "days_to_death", "vital_status")]
  if(i == 1){
    tcga <- x
  } else {
    tcga <- rbind(tcga, x) 
  }
}
tcga <- tcga[which(tcga$vital_status == "Dead" | tcga$vital_status == "Alive"),]
tcga$vital_status <- ifelse(tcga$vital_status == "Dead", 1, 0)
tcga$days_to_last_followup[tcga$days_to_last_followup=="[Not Available]"] <- NA
tcga$days_to_death[tcga$days_to_death=="[Not Applicable]"] <- NA
tcga$TimeVar <- apply(tcga[,2:3], MARGIN = 1, FUN = function(x) max(x, na.rm = T))
tcga <- tcga[-which(tcga$TimeVar<0),]
colnames(tcga)[4] <- 'eventVar'
colnames(tcga)[1] <- 'USI'
tcga <- tcga[,c('USI', 'eventVar', 'TimeVar')]

total <- rbind(target,tcga)
results <- unique(merge(dat, total, by = 'USI'))
results <- results[-which(results$TimeVar == 3123),]
results <- results[!is.na(results$TimeVar),]
write.table(results, file = '~/Projects/toil-rnaseq-20k/data/metadata_filtered/filtered_datasets_with_clinical.txt', quote = F, sep = "\t", row.names = F)
