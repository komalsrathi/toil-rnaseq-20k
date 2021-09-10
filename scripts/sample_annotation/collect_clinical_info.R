############################################
# Author: Komal S Rathi
# Function: This script creates
# clinical files for TCGA and TARGET
############################################

library(reshape2)
setwd('~/Projects/toil-rnaseq-20k/data/metadata/target/clinical')
sample.info <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/sample_info.txt', stringsAsFactors = F)

########################################################################
# TARGET
# ALL (149 patients)
all1 <- read.csv('ALL_Phase1/TARGET_ALL_ClinicalData_Phase1_20160714.csv',check.names = F, stringsAsFactors = F)
all2 <- read.csv('ALL_Phase2/TARGET_ALL_ClinicalData_Phase2_Discovery_20160714.csv',check.names = F, stringsAsFactors = F)
all2.1 <- read.csv('ALL_Phase2/TARGET_ALL_ClinicalData_Phase2_Validation_20160714.csv',check.names = F, stringsAsFactors = F)
all2.2 <- read.csv('ALL_Phase2/TARGET_ALL_ClinicalData_Phase2_Xenografts_20160714.csv',check.names = F, stringsAsFactors = F)
all1 <- all1[,colnames(all2)]
all2.1 <- all2.1[,colnames(all2)]
all2.2 <- all2.2[,colnames(all2)]
all.3 <- rbind(all1, all2, all2.1, all2.2)
all.3$`ALL Molecular Subtype`[all.3$`ALL Molecular Subtype` == "None of the above"] <- NA
all.3 <- unique(all.3)
colnames(all.3) <- sub('[/\n]','',colnames(all.3))
tmp <- all.3[which(duplicated(all.3$`TARGET USI`)),'TARGET USI'] 
tmp <- all.3[which(all.3$`TARGET USI` %in% tmp),]
# write.table(tmp, file = 'ALL_duplicates.txt', quote = F, sep = "\t", row.names = F)

# manually merge duplicates
# so after merging/removing duplicates
tmp <- read.delim('ALL_duplicates_2.txt', check.names = F, stringsAsFactors = F)
all.3 <- all.3[-which(all.3$`TARGET USI` %in% tmp$`TARGET USI`),]
all.4 <- rbind(all.3, tmp)
all.4 <- all.4[which(all.4$`TARGET USI` %in% sample.info$patient_barcode),]
all.4 <- all.4[,colSums(is.na(all.4))<nrow(all.4)] 
all.4$`Cell of Origin` <- "B-Precursor"
all.colstokeep <- 'CNS Status at Diagnosis'
rm(all1, all2, all2.1, all2.2, all.3, tmp)

# AML (186 patients)
aml <- read.csv('AML/TARGET_AML_ClinicalData_20160714.csv',check.names = F, stringsAsFactors = F)
aml <- aml[,colSums(is.na(aml))<nrow(aml)]
aml <- aml[which(aml$`TARGET USI` %in% sample.info$patient_barcode),]
aml.colstokeep <- 'Risk group'

# NBL
nbl <- read.csv('NBL/TARGET_NBL_ClinicalData_20151124.csv', check.names = F, stringsAsFactors = F)
nbl <- nbl[,colSums(is.na(nbl))<nrow(nbl)]
nbl <- nbl[which(nbl$`TARGET USI` %in% sample.info$patient_barcode),]
nbl.colstokeep <- c('INSS Stage', 'MYCN status', 'COG Risk Group')

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
rt.3 <- unique(rbind(rt.1, rt.2))
colnames(rt.3)[16] <- "Histology Classification of Primary Tumor"
rt.3 <- rt.3[which(rt.3$`TARGET USI` %in% sample.info$patient_barcode),]
rt.colstokeep <- 'Stage'
rm(rt.1, rt.2)

# WT
wt.1 <- read.csv('WT/TARGET_WT_ClinicalData_Discovery_20160714_public.csv', stringsAsFactors = F, check.names = F)
wt.2 <- read.csv('WT/TARGET_WT_ClinicalData_Validation_20160714_public.csv', stringsAsFactors = F, check.names = F)
wt.3 <- unique(rbind(wt.1, wt.2))
wt.3 <- wt.3[which(wt.3$`TARGET USI` %in% sample.info$patient_barcode),]
wt.colstokeep <- c('Stage','Histology Classification of Primary Tumor')
rm(wt.1, wt.2)

# add cancer types
all.4$`Cancer Type` <- "ALL"
aml$`Cancer Type` <- "AML"
nbl$`Cancer Type` <- "NBL"
rt.3$`Cancer Type` <- "RT"
wt.3$`Cancer Type` <- "WT"

all.4$`Cancer Type Detail` <- "Childhood Acute Lymphoblastic Leukemia"
aml$`Cancer Type Detail` <- "Childhood Acute Myeloid Leukemia"
nbl$`Cancer Type Detail` <- "Neuroblastoma"
rt.3$`Cancer Type Detail` <- "Rhabdoid Tumor"
wt.3$`Cancer Type Detail` <- "Wilms Tumor"

# cols frequency
d <- data.frame(cols = c(colnames(all.4),colnames(aml),colnames(nbl),colnames(rt.3),colnames(wt.3)), stringsAsFactors = F)
d$cancer <- c(rep('ALL',length(colnames(all.4))), 
              rep('AML',length(colnames(aml))),
              rep('NBL',length(colnames(nbl))),
              rep('RT',length(colnames(rt.3))),
              rep('WT',length(colnames(wt.3))))
d.ct <- plyr::count(d$cols)
d <- merge(d, d.ct, by.x = 'cols', by.y ='x')
common.cols <- unique(d[which(d$freq == 5),'cols'])
common.cols <- common.cols[-grep('Comment|Protocol|Year',common.cols)]

# now subset the columns but keep the ones you want to
aml <- aml[,which(colnames(aml) %in% c(common.cols,aml.colstokeep))]
colnames(aml)[10] <- "AML Risk group"
aml[which(aml$`TARGET USI` %in% sample.info[which(sample.info$study== 'TARGET' & sample.info$disease == "AML-IF"),'patient_barcode']),'Cancer Type'] <- "AML-IF"
aml[which(aml$`TARGET USI` %in% sample.info[which(sample.info$study== 'TARGET' & sample.info$disease == "AML-IF"),'patient_barcode']),'Cancer Type Detail'] <- "Childhood Acute Myeloid Leukemia (Induction Failure)"

all.4 <- all.4[,which(colnames(all.4) %in% c(common.cols,all.colstokeep))]
colnames(all.4)[10] <- 'ALL CNS Status at Diagnosis'

nbl <- nbl[,which(colnames(nbl) %in% c(common.cols, nbl.colstokeep))]
colnames(nbl)[c(10:12)] <- c('NBL INSS Stage','NBL MYCN Status','NBL COG Risk Group')

wt.3 <- wt.3[,which(colnames(wt.3) %in% c(common.cols, wt.colstokeep))]
colnames(wt.3)[10:11] <- c("WT Stage","WT Subtype")

rt.3 <- rt.3[,which(colnames(rt.3) %in% c(common.cols, rt.colstokeep))]
colnames(rt.3)[10] <- c('RT Stage')

# merge everything
allsets <- c('all.4','aml','nbl','rt.3','wt.3')
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
raw[is.na(raw)] <- "Not Applicable"
target.clin <- raw
colnames(target.clin)[grep(paste(common.cols,collapse = "|"),colnames(target.clin))] <- paste0('x ',colnames(target.clin)[grep(paste(common.cols,collapse = "|"),colnames(target.clin))])
target.clin$`x First Event` <- NULL

rm(list=setdiff(ls(), c("target.clin","sample.info")))
target.clin <- merge(sample.info[,c('analysis_id','patient_barcode','definition')], target.clin, by.x = 'patient_barcode', by.y = 'x TARGET USI')
write.table(target.clin, file = '~/Projects/toil-rnaseq-20k/data/clinical/target_clinical_info.txt', quote = F, row.names = F, sep = "\t")

colnames(target.clin)[1:3] <- paste0('x ',colnames(target.clin)[1:3])
cols <- grep('x',colnames(target.clin), value = T)
tmp <- melt(target.clin, id.vars = cols)
tmp$value <- gsub('^ +','',tmp$value)
colnames(tmp) <- gsub('x ','',colnames(tmp))
tmp$variable <- sub('^.* ','',tmp$variable)
tmp$value[tmp$value == "Not Applicable"] <- NA

########################################################################
# TCGA
l <- list.files(path = '~/Projects/toil-rnaseq-20k/data/metadata/tcga/clinical/gdc_download_20161214_200436/', pattern = ".txt", recursive = T, full.names = T)
l <- l[grep('patient',l)]

tcga <- data.frame()
for(i in 1:length(l)){
  x <- read.delim(file = l[i], skip = 1, check.names = F, stringsAsFactors = F)
  x <- x[-1,]
  if(i == 1){
    tcga <- x
  } else {
    tcga <- plyr::rbind.fill(tcga, x) 
  }
}
tcga.clin <- tcga
tcga.clin <- tcga.clin[which(tcga.clin$bcr_patient_barcode %in% sample.info$patient_barcode),]
tcga.clin <- tcga.clin[,colSums(is.na(tcga.clin))<nrow(tcga.clin)]
tcga.clin <- tcga.clin[,c('bcr_patient_barcode','gender','race','ethnicity',
                          'age_at_initial_pathologic_diagnosis',
                          'vital_status','days_to_death','pathologic_stage','tumor_tissue_site')]
colnames(tcga.clin) <- c('patient_barcode','Gender','Race','Ethnicity','Age at Diagnosis in Days','Vital Status','Overall Survival Time in Days','Stage','Tumor Tissue Site')
tcga.clin <- merge(sample.info[,c('analysis_id','patient_barcode','definition')], tcga.clin, by = 'patient_barcode')
write.table(tcga.clin, file = '~/Projects/toil-rnaseq-20k/data/clinical/tcga_clinical_info.txt', quote = F, row.names = F, sep = "\t")
