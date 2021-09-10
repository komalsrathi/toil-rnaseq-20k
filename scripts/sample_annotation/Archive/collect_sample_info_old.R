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
write.table(gtex.s, file = '~/Projects/toil-rnaseq-20k/data/clinical/gtex_metadata_7863.txt', quote = F, sep = "\t", row.names = F)

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
target.all$sample_type_name <- target.all$definition
target.all$definition <- NULL
target.all <- target.all[,colnames(tcga.all)]
write.table(target.all, file = '~/Projects/toil-rnaseq-20k/data/clinical/target_metadata_721.txt', quote = F, sep = "\t", row.names = F)
