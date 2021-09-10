setwd('~/Projects/toil-rnaseq-20k/data/metadata/target/to_process/')
library(dplyr)

sra <- read.delim('data/SraRunTable_toprocess.txt')
sra <- unique(sra[,c('Sample_Name_s','Run_s','body_site_s','study_name_s')])
sra$study <- gsub('[-][A-Z0-9]+[-][0-9A-Z]{3}-[0-9A-Z]{3}','',sra$Sample_Name_s)

sra$study_name_s <- gsub('.*[(]|[)].*','',sra$study_name_s)
sra[grep('TARGET-15',sra$Sample_Name_s),'study_name_s'] <- 'MPAL'
sra[grep('TARGET-21',sra$Sample_Name_s),'study_name_s'] <- 'AML-IF'
sra$body_site_s_short <- gsub("(?<=[A-Z-])[^A-Z-]+", "", sra$body_site_s, perl = TRUE)
sra$body_site_s_short <- sub('-','_',sra$body_site_s_short)
meta <- sra[,c('Run_s','study_name_s','body_site_s','body_site_s_short','Sample_Name_s')]
meta$study <- 'TARGET'
colnames(meta) <- c('analysis_id','disease','definition','definition_short','barcode','study')
# meta <- meta[which(meta$disease != 'AML-IF'),]
ct <- plyr::count(meta, c('disease','definition_short'))

# add RT
all.tumors <- read.delim('data/SraRunTable_allTumors.txt', stringsAsFactors = F)
target.RT <- read.delim('~/Projects/Maris-lab/Kris/TARGET_RT/data/TARGET_RT_metadata.txt')
target.RT <- all.tumors[which(all.tumors$Run_s %in% target.RT$analysis_id),]
target.RT <- target.RT[,c('Run_s','Sample_Name_s','body_site_s')]
colnames(target.RT) <- c('analysis_id','barcode','definition')
target.RT$definition_short <- 'PST'
target.RT$study <- 'TARGET'
target.RT$disease <- 'RT'
target.RT <- target.RT[,colnames(meta)]

meta <- rbind(meta, target.RT)

# doesn't look like it matters if RNA-seq or WES/WGS etc
res <- all.tumors[which(all.tumors$Sample_Name_s %in% meta$barcode),]
res$duplicated[which(duplicated(res$Sample_Name_s))] <- "Yes"
res <- res[,c('Assay_Type_s','body_site_s','molecular_data_type_s','analyte_type_s','Sample_Name_s')]

# so continue
res <- meta[!duplicated(meta$barcode),]
res <- res[-which(res$disease %in% c("AML-IF","MPAL")),]
write.table(res, file = 'data/TARGET_samples_toadd.txt', quote = F, sep = "\t", row.names = F)

# new.meta <- meta
# new.meta$dups[which(duplicated(new.meta$barcode))] <- 'Yes'
# dups <- unique(as.character(new.meta[which(new.meta$dups == "Yes"),'barcode']))
# new.meta <- new.meta[which(new.meta$barcode %in% dups),]
# new.meta$dups <- NULL
# write.table(new.meta, file = '~/Projects/toil-rnaseq-20k/data/metadata/to_process/TARGET_duplicate_barcodes.txt', quote = F, sep = "\t", row.names = F)

# weird examples
# egs <- c('TARGET-50-CAAAAQ-01A-01R','TARGET-52-PADYCE-01A-01R','TARGET-40-PASKZZ-01A-01R')
# all.tumors.sub <- all.tumors[which(all.tumors$Sample_Name_s %in% egs),c('Sample_Name_s','Run_s','study_name_s','ReleaseDate_s','LoadDate_s','MBases_l','MBytes_l')]
# 
# # write out whole sra file
# dat <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/to_process/TARGET_duplicate_barcodes.txt', stringsAsFactors = F)
# dat <- dat[order(dat$barcode),]
# 
# sra <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_allTumors.txt')
# sra <- sra[which(sra$Sample_Name_s %in% dat$barcode),]
# sra <- sra[order(sra$Sample_Name_s),]
# write.table(sra, file = '~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_TARGET_duplicate_barcodes.txt', quote = F, sep = "\t", row.names = F)
# 
# library(dplyr)
# sra <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_allTumors.txt', stringsAsFactors = F)
# sra.dups <- sra[duplicated(sra$Sample_Name_s),'Sample_Name_s']
# sra.dups <- unique(sra.dups)
# sra.dups <- sra[which(sra$Sample_Name_s %in% sra.dups),]
# sra.dups <- sra.dups[order(sra.dups$Sample_Name_s),]
# write.table(sra.dups, file = '~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_allTumors_duplicates.txt', quote = F, sep = "\t", row.names = F)
# 
# unique(sra.dups$study_name_s)
# nbl <- sra.dups[grep('NBL',sra.dups$study_name_s),]
# nbl <- nbl[order(nbl$Sample_Name_s),]
# write.table(nbl, file = '~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_NBL_duplicates.txt', quote = F, sep = "\t", row.names = F)
# 
# egs <- c('TARGET-50-CAAAAQ-01A-01R','TARGET-52-PADYCE-01A-01R','TARGET-40-PASKZZ-01A-01R')
# sra.weird <- sra[which(sra$Sample_Name_s %in% egs),c('Sample_Name_s','Run_s','study_name_s','ReleaseDate_s','LoadDate_s','MBases_l','MBytes_l')]
# sra.weird <- sra.weird[order(sra.weird$Sample_Name_s),]
# write.table(sra.weird, file = '~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_duplicates_OtherExamples.txt', quote = F, sep = "\t", row.names = F)
# 
# # to resolve duplicates
# sra.tmp <- sra[,c('Sample_Name_s','LoadDate_s','MBases_l','MBytes_l')]
# ct <- plyr::count(sra,c('Sample_Name_s','LoadDate_s','MBases_l'))
# sra.tmp <- sra.tmp %>% group_by(Sample_Name_s) %>% filter(LoadDate_s == max(as.Date(LoadDate_s))) %>% 
#   group_by(Sample_Name_s, LoadDate_s) %>% 
#   filter(MBases_l == max(MBases_l)) %>% 
#   filter(MBytes_l == max(MBytes_l))
# xx <- plyr::count(sra.tmp, 'Sample_Name_s')
# 
# # conflict between Molecular Type and Assay Type
# sra <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_allTumors.txt', stringsAsFactors = F)
# sra.desc <- sra[which(sra$Assay_Type_s == "RNA-Seq" & sra$molecular_data_type_s != "RNA Seq (NGS)" & sra$molecular_data_type_s != "<not provided>"),]
# write.table(sra.desc, file = '~/Projects/toil-rnaseq-20k/data/metadata/to_process/SraRunTable_MolecularDataType_Discrepancy.txt', quote = F, sep = "\t", row.names = F)


# As plan B, I can resolve this duplication programmatically. These samples are unique based on the combination of 4 columns: Barcode, LoadDate_s, MBases_l and MBytes_l. What I can do is take the most recent ones based on LoadDate_s and then take the ones with max MBases_l and max MBytes_l. This will give me a unique set of TARGET Barcode-SRA ID map. 

