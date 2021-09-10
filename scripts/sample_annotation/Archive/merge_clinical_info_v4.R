###################################
# Author: Komal S Rathi
# Date: 10/16/2018
# Function: Merge all clinical data + new CBTTC
###################################

setwd('~/Projects/toil-rnaseq-20k/')

# new clinical data for cbttc
dat <- read.delim('data/clinical/clinical_info.txt', stringsAsFactors = F)
dat[is.na(dat)] <- ''
plyr::count(dat$study)
cbttc.meta <- dat[which(dat$study == "CBTTC"),]

# add cbttc file
cbttc <- read.delim('data/metadata/cbttc/cbttc_metadata_v3.txt', stringsAsFactors = F)
cbttc[is.na(cbttc)] <- ''

setdiff(colnames(cbttc), colnames(cbttc.meta))
setdiff(colnames(cbttc.meta), colnames(cbttc))
cbttc <- cbttc[,colnames(cbttc.meta)]
for(i in 1:ncol(cbttc)){
  print(colnames(cbttc)[i])
  set1 <- unique(cbttc[,i])
  set2 <- unique(cbttc.meta[,i])
  print(setdiff(set1, set2))
}
dat <- rbind(dat, cbttc)
write.table(dat, file = '~/Projects/toil-rnaseq-20k/data/clinical/final_clinical/clinical_info.txt', quote = F, sep = "\t", row.names = F)
