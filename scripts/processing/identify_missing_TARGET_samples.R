setwd('~/Projects/toil-rnaseq-20k/data/metadata/target/to_process/')

dat <- read.delim('data/SraRunTable_allTumors.txt', stringsAsFactors = F)
dat <- unique(dat[,c('Assay_Type_s','Center_Name_s','Run_s','Sample_Name_s','analyte_type_s','body_site_s','histological_type_s','is_tumor_s','molecular_data_type_s','study_name_s')])
dat <- dat[grep('Primary|Recurrent',dat$body_site_s),]

target <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/target/target_metadata.txt')
dat <- dat[-which(dat$Sample_Name_s %in% target$barcode),]
dat$link <- paste0('ln -s /mnt/isilon/diskin_lab/TARGET_RNA/fastq/',dat$Run_s,'_* ./')

RT <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/target/TARGET_RT_metadata.txt')
tmp <- dat[which(dat$Run_s %in% RT$analysis_id),]
dat <- dat[-which(dat$Run_s %in% RT$analysis_id),]

d_ply(dat, .variables = 'study_name_s', .fun = function(x) write.t(x))
write.t <- function(x){
  study <- unique(gsub('.*[(]|[)].*','',x$study_name_s))
  n <- nrow(x)
  name <- paste0('~/Projects/toil-rnaseq-20k/data/metadata/target/to_process/data/',study,'_',n,'.txt')
  write.table(x, file = name, row.names = F, quote = F, sep = "\t")
}

# TARGET-15 is actually MPAL
write.table(dat, file = 'data/SraRunTable_toprocess.txt', quote = F, row.names = F, sep = "\t")
