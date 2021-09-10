###################################
# Author: Komal S Rathi
# Date: 08/20/2018
# Function: Merge all clinical data + new CBTTC
###################################

setwd('~/Projects/toil-rnaseq-20k/')

# new target data added
clin <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/clinical_info_v2.txt', stringsAsFactors = F)

# new clinical data for cbttc
cbttc <- read.csv('~/Projects/toil-rnaseq-20k/data/metadata/cbttc/cbttc_metadata_v2.csv', stringsAsFactors = F)
clin <- rbind(clin, cbttc)
clin[is.na(clin)] <- ''
clin <- clin[!clin$study == "MB_TaylorLab",]
write.table(clin, file = '~/Projects/toil-rnaseq-20k/data/clinical/clinical_info.txt', quote = F, sep = "\t", row.names = F)
