#############################################################
# Author: Komal S Rathi
# Institute: CHOP
# Function: This script adds cancer subtypes to TCGA metadata
#############################################################

library(dplyr)

# read full datasets metafile
all.samples <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata_filtered/full_datasets.txt', stringsAsFactors = F)
all.samples <- all.samples[-which(all.samples$disease == "CNTL"),]

# add subtypes
allclin <- data.frame()
lfiles <- list.files(path = '~/Projects/toil-rnaseq-20k/data/metadata/tcga/clinical/',pattern = "_clinical_data.tsv",full.names = T)
for(i in 1:length(lfiles)){
  datp <- read.delim(lfiles[i])
  datp$Cancer.Studies <- toupper(sub('_tcga.*','',datp$Cancer.Studies))
  datp <- datp[,c('Sample.ID','Cancer.Type.Detailed','Cancer.Studies')]
  allclin <- rbind(allclin, datp)
}
allclin$Cancer.Studies[allclin$Cancer.Type.Detailed == "Rectal Adenocarcinoma"] <- "READ"
allclin$Cancer.Studies[allclin$Cancer.Type.Detailed == "Colon Adenocarcinoma"] <- "COAD"

allclin <- unique(allclin)
all.samples$newbarcode <- sub('[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[0-9]{2}','', all.samples$barcode)
dat <- merge(all.samples, allclin, by.x = c('newbarcode','disease'), by.y = c('Sample.ID','Cancer.Studies'), all.x = T)
colnames(dat)[13] <- 'cancer_subtype'
dat <- unique(dat)
dat$cancer_subtype <- as.character(dat$cancer_subtype)

# for all empty subtypes, replace with max sample subtypes
counts <- plyr::count(dat[which(dat$study == 'TCGA'),], c('disease','cancer_subtype'))
w <- counts %>% 
  group_by(disease) %>% 
  filter(freq == max(freq)) %>% 
  unique %>% as.data.frame()
w$freq <- NULL
for(i in 1:nrow(dat)){
  if(dat[i,'study'] == "TCGA" & is.na(dat[i,"cancer_subtype"])){
    disease <- dat[i,'disease']
    print(disease)
    subtype <- w[which(w$disease == disease),'cancer_subtype']
    print(subtype)
    dat[i,'cancer_subtype'] <- subtype
    print(subtype)
  }
}

# add subtypes for other studies (just use disease name where unavailable)
dat$cancer.subtype <- ifelse(dat$study %in% c('TARGET','MB','SCLC'), dat$disease_name, dat$cancer_subtype)
final <- dat[,c('analysis_id','disease','definition','study','group','disease_name','cancer_subtype')]
write.table(final, file = '~/Projects/toil-rnaseq-20k/data/metadata/metadata_with_cancersubtypes.txt', quote = F, sep = "\t", row.names = F)
