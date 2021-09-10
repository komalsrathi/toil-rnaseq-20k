###################################
# Author: Komal S Rathi
# Date: 08/03/2018
# Function: Create tags for TARGET
###################################

library(reshape2)
library(plyr)
library(dplyr)

setwd('~/Projects/toil-rnaseq-20k/data/metadata/target/clinical')
sample.info <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/sample_info.txt', stringsAsFactors = F)

# this is to merge duplicates from ALL 
toString_all <- function(x){
  print(x)
  if(length(grep('[A-Za-z]', x)) > 0){
    print('First')
    y.str <- toString(na.omit(unique(x)))
    z.str <- unlist(strsplit(y.str, ', '))
    if(length(z.str[which.max(nchar(z.str))]) > 0){
      z.str[which.max(nchar(z.str))]
    } else {
      ''
    }
  } else if(all(is.na(x))){
    print('Second')
    ''
  } else {
    print('Third')
    as.character(max(as.numeric(na.omit(unique(x)))))
  }
}

########################################################################
# TARGET
# ALL (n = 149 patients)
all1 <- read.csv('ALL_Phase1/TARGET_ALL_ClinicalData_Phase1_20160714.csv', check.names = F, stringsAsFactors = F)
all2 <- read.csv('ALL_Phase2/TARGET_ALL_ClinicalData_Phase2_Discovery_20160714.csv', check.names = F, stringsAsFactors = F)
all2.1 <- read.csv('ALL_Phase2/TARGET_ALL_ClinicalData_Phase2_Validation_20160714.csv', check.names = F, stringsAsFactors = F)
all2.2 <- read.csv('ALL_Phase2/TARGET_ALL_ClinicalData_Phase2_Xenografts_20160714.csv', check.names = F, stringsAsFactors = F)
all1 <- all1[,colnames(all2)]
all2.1 <- all2.1[,colnames(all2)]
all2.2 <- all2.2[,colnames(all2)]
all <- unique(rbind(all1, all2, all2.1, all2.2))
all <- unique(all[,c(1:9,14)])
clean.up <- function(x){
  if(length(grep('Dead|Male|White',x, ignore.case = T)) > 0){
    x <- tools::toTitleCase(tolower(x))
  }
  x <- trimws(x)
  x[x %in% c("Not Reported", "Unknown", "", "None", ".", "Not reported")] <- "unavailable"
  return(x)
}
all <- as.data.frame(apply(all, 2, clean.up), stringsAsFactors = F)
dups <- all[which(duplicated(all$`TARGET USI`)),'TARGET USI'] 
dups <- all[which(all$`TARGET USI` %in% dups),]
dups[dups == "unavailable"] <- NA
dups <- dups %>% group_by(`TARGET USI`) %>% 
  summarise_all(toString_all) %>% 
  as.data.frame(stringsAsFactors= F)
all <- all[-which(all$`TARGET USI` %in% dups$`TARGET USI`),]
all <- rbind(all, dups)
all <- all[which(all$`TARGET USI` %in% sample.info$patient_barcode),]
all <- all[,colSums(is.na(all))<nrow(all)] 
all.colstokeep <- 'CNS Status at Diagnosis'
all$disease_subtype <- 'B-Precursor'
rm(all1, all2, all2.1, all2.2)

# AML (n = 186 patients)
aml <- read.csv('AML/TARGET_AML_ClinicalData_20160714.csv',check.names = F, stringsAsFactors = F)
aml <- aml[,colSums(is.na(aml))<nrow(aml)]
aml <- aml[which(aml$`TARGET USI` %in% sample.info$patient_barcode),]
aml$disease_subtype <- aml$`FAB Category`
aml.colstokeep <- 'Risk group'
aml[95,'ISCN'] <- "46,XY,9ph,t(16;17)(p13.1;p13.3)"
aml <- as.data.frame(apply(aml, 2, clean.up), stringsAsFactors = F)

# NBL (n = 155)
nbl <- read.csv('NBL/TARGET_NBL_ClinicalData_20151124.csv', check.names = F, stringsAsFactors = F)
nbl <- nbl[,colSums(is.na(nbl))<nrow(nbl)]
nbl <- nbl[which(nbl$`TARGET USI` %in% sample.info$patient_barcode),]
nbl.colstokeep <- c('INSS Stage', 'MYCN status', 'COG Risk Group')
nbl$`MYCN status`[nbl$`MYCN status` == "Not Amplified"] <- 'Non-amplified'
nbl$disease_subtype <- nbl$`Diagnostic Category`
nbl <- as.data.frame(apply(nbl, 2, clean.up), stringsAsFactors = F)

# OS not required
# os.1 <- read.csv('OS/TARGET_OS_ClinicalData_Discovery_20150729.csv', check.names = F, stringsAsFactors = F)
# os.2 <- read.csv('OS/TARGET_OS_ClinicalData_Validation_20160401.csv', check.names = F, stringsAsFactors = F)
# os.3 <- unique(rbind(os.1, os.2))
# os.3$Comment[os.3$Comment==""] <- NA
# os.3 <- unique(os.3)
# os.3 <- os.3[,colSums(is.na(os.3))<nrow(os.3)]
# colnames(os.3)[7] <- "Event Free Survival Time in Days"
# os.3 <- os.3[which(os.3$`TARGET USI` %in% sample.info$patient_barcode),]
# rm(os.1, os.2)

# RT
rt.1 <- read.csv('RT/TARGET_RT_ClinicalData_Discovery_20150710_public.csv', stringsAsFactors = F, check.names = F)
rt.2 <- read.csv('RT/TARGET_RT_ClinicalData_Validation_20150710_public.csv', stringsAsFactors = F, check.names = F)
rt <- unique(rbind(rt.1, rt.2))
rt <- rt[which(rt$`TARGET USI` %in% sample.info$patient_barcode),]
rt.colstokeep <- 'Stage'
rt$disease_subtype <- 'Rhabdoid Tumor'
rt <- as.data.frame(apply(rt, 2, clean.up), stringsAsFactors = F)
rm(rt.1, rt.2)

# WT
wt.1 <- read.csv('WT/TARGET_WT_ClinicalData_Discovery_20160714_public.csv', stringsAsFactors = F, check.names = F)
wt.2 <- read.csv('WT/TARGET_WT_ClinicalData_Validation_20160714_public.csv', stringsAsFactors = F, check.names = F)
wt <- unique(rbind(wt.1, wt.2))
wt <- wt[which(wt$`TARGET USI` %in% sample.info$patient_barcode),]
wt$disease_subtype <- wt$`Histology Classification of Primary Tumor`
wt.colstokeep <- c('Stage')
wt <- as.data.frame(apply(wt, 2, clean.up), stringsAsFactors = F)
rm(wt.1, wt.2)

# add cancer types
all$`Cancer Type` <- "ALL"
aml$`Cancer Type` <- "AML"
nbl$`Cancer Type` <- "NBL"
rt$`Cancer Type` <- "RT"
wt$`Cancer Type` <- "WT"

all$`Cancer Type Detail` <- "Acute Lymphoblastic Leukemia"
aml$`Cancer Type Detail` <- "Acute Myeloid Leukemia"
nbl$`Cancer Type Detail` <- "Neuroblastoma"
rt$`Cancer Type Detail` <- "Rhabdoid Tumor"
wt$`Cancer Type Detail` <- "Wilms Tumor"

# get common columns 
common.cols <- data.frame(cols = c(colnames(all),colnames(aml),colnames(nbl),colnames(rt),colnames(wt)), stringsAsFactors = F)
common.cols$cancer <- c(rep('ALL',length(colnames(all))), 
              rep('AML',length(colnames(aml))),
              rep('NBL',length(colnames(nbl))),
              rep('RT',length(colnames(rt))),
              rep('WT',length(colnames(wt))))
common.cols.ct <- plyr::count(common.cols$cols)
common.cols <- merge(common.cols, common.cols.ct, by.x = 'cols', by.y ='x')
common.cols <- unique(common.cols[which(common.cols$freq == 5),'cols'])

# now subset the columns but keep the ones you want to
aml <- aml[,which(colnames(aml) %in% c(common.cols, aml.colstokeep))]
colnames(aml)[10] <- "risk"
aml[grep('TARGET-21', aml$`TARGET USI`),'Cancer Type'] <- "AML-IF"
aml[grep('TARGET-21', aml$`TARGET USI`),'Cancer Type Detail'] <- "Acute Myeloid Leukemia - Induction Failure"
aml$risk <- tolower(aml$risk)

all <- all[,which(colnames(all) %in% c(common.cols,all.colstokeep))]
colnames(all)[10] <- 'cns_stage'
all$cns_stage <- gsub('CNS ','',all$cns_stage)

nbl <- nbl[,which(colnames(nbl) %in% c(common.cols, nbl.colstokeep))]
colnames(nbl)[c(10:12)] <- c('stage','mycn_status','risk')
nbl$stage <- gsub('.* ','',nbl$stage)
nbl$risk <- gsub(' .*','', tolower(nbl$risk))
nbl$mycn_status <- tolower(nbl$mycn_status)

wt <- wt[,which(colnames(wt) %in% c(common.cols, wt.colstokeep))]
colnames(wt)[10] <- "stage"

rt <- rt[,which(colnames(rt) %in% c(common.cols, rt.colstokeep))]
colnames(rt)[10] <- 'stage'

# merge everything
allsets <- c('all','aml','nbl','rt','wt')
for(i in 1:length(allsets)){
  if(i==1){
    raw <- get(allsets[i])
  }
  if(i>1){
    tmp <- get(allsets[i])
    raw <- plyr::rbind.fill(raw, tmp)
  }
}
raw <- raw[,colSums(is.na(raw))<nrow(raw)]
target.clin <- raw
target.clin$`First Event` <- NULL
rm(list=setdiff(ls(), c("target.clin","sample.info")))
target.clin <- merge(sample.info[,c('analysis_id','patient_barcode','definition')], target.clin, by.x = 'patient_barcode', by.y = 'TARGET USI')

# now tag information
tags <- melt(target.clin, measure.vars = c('cns_stage','mycn_status','risk','stage'), na.rm = TRUE, variable.name = 'tags')
tags$tags <- paste0(tags$tags,'=',tags$value)
tags$value <- NULL
target.clin <- ddply(.data = tags, 
                 .variables = .(patient_barcode, analysis_id, definition,
                                Gender, Race, Ethnicity,
                                `Age at Diagnosis in Days`,
                                `Event Free Survival Time in Days`,
                                `Vital Status`,`Overall Survival Time in Days`,
                                `Cancer Type`,`Cancer Type Detail`, disease_subtype), 
                 .fun = summarize, 
                 tags = toString(tags))
target.clin$tags <- gsub(', ','; ',target.clin$tags)
colnames(target.clin) <- gsub(' ','_',tolower(colnames(target.clin)))
target.clin$patient_barcode <- NULL
target.clin$definition <- NULL
target.clin$cancer_type <- NULL
target.clin$cancer_type_detail <- NULL
colnames(target.clin)[1] <- 'sample_id'

rm(sample.info, tags)
save(target.clin, file = '~/Projects/toil-rnaseq-20k/data/clinical/targetclin.RData')
