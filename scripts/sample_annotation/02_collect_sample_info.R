##############################################
# Author: Komal S Rathi
# Function: To get sample information for 
# TCGA, TARGET and GTEx from UCSC meta files
##############################################
setwd('~/Projects/toil-rnaseq-20k/data/metadata/')

# GTEx
# n = 7924 total
# 61 are improperly paired so we don't use them -> n = 7863
# 4 removing those with histological_type == <not provided> -> # not applicable anymore
gtex <- read.delim('gtex/gtex-manifest', header = F) 
gtex$V1 <- gsub('s3://cgl-rnaseq-recompute-fixed/gtex/|.tar.gz','',gtex$V1)
gtex$V2 <- 'PAIRED'
gtex[grep("IMPROPER",gtex$V1),'V2'] <- "IMPROPERLY PAIRED"
gtex <- gtex[-which(gtex$V2 == "IMPROPERLY PAIRED"),]
gtex.s <- read.delim('gtex/SraRunTable.txt', stringsAsFactors = F)
gtex.s <- gtex.s[which(gtex.s$Run_s %in% gtex$V1),]
gtex.s$group <- 'Normals'
gtex.s <- gtex.s[,c('Run_s','submitted_subject_id_s','submitted_sample_id_s','group','study_name_s','body_site_s','histological_type_s')]
gtex.s$definition <- "Normals"
colnames(gtex.s) <- c('analysis_id','patient_barcode','sample_barcode','group','study','tissue','subtissue','definition')
gtex.s$disease <- ''
gtex.s$disease_name <- ''
gtex.s$disease_subtype <- ''
write.table(gtex.s, file = '~/Projects/toil-rnaseq-20k/data/clinical/gtex_metadata_7863.txt', quote = F, sep = "\t", row.names = F)
rm(gtex)

# TCGA
# n = 11313
# n = 644 SINGLE ENDED, 6 IMPROPERLY PAIRED, 1 NO META = 651
# 10662 usable samples
tcga <- read.delim('tcga/tcga-manifest', stringsAsFactors = F, header = F)
tcga$V1 <- gsub('s3://cgl-rnaseq-recompute-fixed/tcga/|.tar.gz','',tcga$V1)
tcga$V2 <- gsub('[.].*','',tcga$V1)
tcga <- tcga[-which(tcga$V2 %in% c('SINGLE-END','IMPROPERLY_PAIRED')),]
tcga.s <- read.delim('tcga/unaligned.tsv')
tcga.s <- tcga.s[which(tcga.s$analysis_id %in% tcga$V1),]
tcga.t <- read.delim('target/LATEST_MANIFEST.tsv')
tcga.t <- tcga.t[which(tcga.t$analysis_id %in% tcga$V1),]
tcga.s <- tcga.s[-which(tcga.s$analysis_id %in% tcga.t$analysis_id),]
tcga.t <- tcga.t[,which(colnames(tcga.t) %in% colnames(tcga.s))]
tcga.all <- unique(rbind(tcga.s, tcga.t))
tcga.all <- tcga.all[which(tcga.all$analysis_id %in% tcga$V1),]
tcga.all$sample_type_name[tcga.all$sample_type_name == '1'] <- "Primary Solid Tumor"
rm(tcga.s, tcga.t, tcga)
tcga.all <- tcga.all[,c('analysis_id','barcode','disease','disease_name','sample_type_name')]
tcga.all$group <- "Tumors"
tcga.all$study <- "TCGA"
tcga.all$tissue <- ""
tcga.all$subtissue <- ""
tcga.all$patient_barcode <- gsub('-[0-9]{2}[A-Z]{1}-[0-9]{2}[A-Z]{1}|-[0-9A-Z]{3}-[0-9A-Z]{3}-[A-Z0-9]{4}-[0-9]{2}|-[0-9]{4}.*-SM-[0-9A-Z]{5}|-SM-[0-9A-Z]{5}','',tcga.all$barcode)
colnames(tcga.all)[c(2,5)] <- c('sample_barcode','definition')
write.table(tcga.all, file = '~/Projects/toil-rnaseq-20k/data/clinical/tcga_metadata_10662.txt', quote = F, sep = "\t", row.names = F)

# TARGET 
# n = 734
# n = 13 NO META 
# n = 721
target <- read.delim('target/target-manifest', stringsAsFactors = F, header = F)
target$V1 <- gsub('s3://cgl-rnaseq-recompute-fixed/target/|.tar.gz','',target$V1)
target.t <- read.delim('target/LATEST_MANIFEST.tsv', stringsAsFactors = F)
target.all <- target.t[which(target.t$analysis_id %in% target$V1),]

target_code <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/tcga_target_codetable.txt')
unique(target.all$sample_type)
target.all$sample_type <- gsub('TARGET_0|TARGET_','',target.all$sample_type)
target.all <- merge(target.all, target_code[,1:2], by.x = 'sample_type',by.y='code')
target.all$group <- "Tumors"
target.all$study <- "TARGET"
target.all$tissue <- ""
target.all$subtissue <- ""
target.all <- target.all[,c('analysis_id','barcode','disease','disease_name','definition','group','study','tissue','subtissue')]
target.all$patient_barcode <- gsub('-[0-9]{2}[A-Z]{1}-[0-9]{2}[A-Z]{1}|-[0-9A-Z]{3}-[0-9A-Z]{3}-[A-Z0-9]{4}-[0-9]{2}|-[0-9]{4}.*-SM-[0-9A-Z]{5}|-SM-[0-9A-Z]{5}','',target.all$barcode)
colnames(target.all)[2] <- 'sample_barcode'
target.all <- target.all[,colnames(tcga.all)]
target.all$disease_name <- stringi::stri_trans_totitle(target.all$disease_name)
target.all$disease_name[target.all$disease_name == "Acute Myeloid Leukemia (Non-Tcga)"] <- "Acute Myeloid Leukemia"
write.table(target.all, file = '~/Projects/toil-rnaseq-20k/data/clinical/target_metadata_721.txt', quote = F, sep = "\t", row.names = F)

# TARGET additional data
target.new <- read.delim('target/target_new_metadata.txt', stringsAsFactors = F)
target.new$sample_type <- gsub('TARGET-[0-9]{2}-[A-Z0-9]+-','',target.new$sample_id)
target.new$sample_type <- gsub("^0|[AB]$","",target.new$sample_type)
target.RT <- read.delim('target/SraRunTable_RT.txt')
target.RT$study <- "TARGET"
target.RT$barcode <- target.RT$Sample_Name_s
target.RT$disease <- "RT"
target.RT$disease_name <- "Rhabdoid Tumor"
target.RT$sample_type <- gsub('TARGET-[0-9]{2}-[A-Z0-9]+-','',target.RT$barcode)
target.RT$sample_type <- gsub("-.*","",target.RT$sample_type)
target.RT$sample_type <- gsub("^0|A$","",target.RT$sample_type)
target.RT$analysis_id <- target.RT$Run_s
target.RT$sample_id <- target.RT$Sample_Name_s
target.RT <- target.RT[,colnames(target.new)]
target.new <- unique(rbind(target.new, target.RT))

target_code <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/tcga_target_codetable.txt')
unique(target.new$sample_type)
target.new <- merge(target.new, target_code[,1:2], by.x = 'sample_type', by.y = 'code')
head(target.new)
target.new$tissue <- ""
target.new$subtissue <- ""
target.new$patient_barcode <- gsub('-[0-9]{2}[A-Z]{1}-[0-9]{2}[A-Z]{1}|-[0-9A-Z]{3}-[0-9A-Z]{3}-[A-Z0-9]{4}-[0-9]{2}|-[0-9]{4}.*-SM-[0-9A-Z]{5}|-SM-[0-9A-Z]{5}','',target.new$barcode)
target.new$sample_type <- NULL
target.new$group <- "Tumors"
colnames(target.new)[2] <- "sample_barcode"
target.new <- target.new[,colnames(tcga.all)]
rm(target.t, target)
target.new <- unique(target.new)
write.table(target.new, file = '~/Projects/toil-rnaseq-20k/data/clinical/target_metadata_1351.txt', quote = F, sep = "\t", row.names = F)

