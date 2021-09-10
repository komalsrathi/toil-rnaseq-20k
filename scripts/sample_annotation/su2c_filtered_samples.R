##############################################################
# Author: Komal Rathi
# Institute: CHOP
# Function: Filter all the usable samples for SU2C project
##############################################################

library(dplyr)

code <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/tcga_target_codetable.txt', stringsAsFactors = F)
code[1:9,'code'] <- paste0('0',code[1:9,'code'])

# 7863 GTEx
gtex <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/gtex/gtex_metadata.txt', stringsAsFactors = F)
gtex$disease <- gtex$histological_type_s
gtex$disease_name <- gtex$disease
gtex <- gtex[,c('study','barcode','disease','disease_name','analysis_id')]
gtex$definition <- 'Normals'
gtex.ct <- plyr::count(gtex, c('disease','definition'))
gtex <- merge(gtex, gtex.ct, by = c('disease','definition'))
rm(gtex.ct)

# 721 -> 675 TARGET
target <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/target/target_metadata.txt', stringsAsFactors = F)
target$code <- sub('TARGET_','',target$sample_type)
target <- merge(target, code, by = 'code', all.x = T)
target.ct <- plyr::count(target, c('disease','definition'))
target <- merge(target, target.ct, by = c('disease','definition'))
target <- target[,colnames(gtex)]
target$disease_name <- tools::toTitleCase(tolower(target$disease_name))
rm(target.ct)

# 10662 -> TCGA
tcga <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/tcga/tcga_metadata.txt', stringsAsFactors = F)
tcga <- merge(tcga, code, by = 'sample_type') 
tcga.ct <- plyr::count(tcga, c('disease','definition'))
tcga <- merge(tcga, tcga.ct, by = c('disease','definition'))
tcga <- tcga[,colnames(gtex)]
tcga$disease_name <- tools::toTitleCase(tolower(tcga$disease_name))
rm(tcga.ct)

# SCLC 86 -> 79
sclc <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/sclc/SCLC_metadata.txt', stringsAsFactors = F)
sclc <- unique(sclc[,c('Assay.Name','Comment..ENA_RUN.','Comment..Sample_source_name.','Characteristics..disease.status.')])
sclc <- sclc[-which(sclc$Comment..Sample_source_name. == "normal"),]
sclc$Comment..Sample_source_name. <- 'Primary Solid Tumor'
colnames(sclc) <- c('barcode','analysis_id', 'definition','disease')
sclc$study <- 'SCLC_GSE60052'
sclc$freq <- 79
sclc$disease <- 'SCLC'
sclc$disease_name <- 'Small Cell Lung Cancer'
sclc <- sclc[,colnames(gtex)]

# Medullo 
medullo <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/medullo/MB_metadata.txt', stringsAsFactors = F)
medullo$disease <- 'MB'
colnames(medullo)[1] <- 'barcode'
medullo$disease_name <- 'Medulloblastoma'
medullo <- medullo[,colnames(gtex)]
medullo$study <- "MB_TaylorLab"

# full dataset
total <- rbind(gtex, medullo, sclc, tcga, target)
total$definition[total$definition == "Additional - New Primary"] <- "Primary Solid Tumor"
total$definition[total$definition == "Additional Metastatic"] <- "Metastatic"
total$group <- ifelse(total$study == "GTEx", "Normals", "Tumors")
colnames(total)[grep('freq',colnames(total))] <- 'total.by.definition'
total <- total %>% group_by(disease) %>% mutate(total.by.disease = n())
total <- total %>% group_by(study) %>% mutate(total.by.study = n())
total$def.short <- sub('-','_',gsub('[a-z]| ','',total$definition))
write.table(total, file = '~/Projects/toil-rnaseq-20k/data/metadata_filtered/full_datasets.txt', quote = F, row.names = F, sep = "\t")

# filtered datasets
filtered <- total %>%
  filter(total.by.definition > 1) %>%
  filter(!definition %in% c("Solid Tissue Normal","Metastatic")) %>%
  filter(disease != "AML-IF")
filtered <- filtered %>% group_by(disease) %>% mutate(total.by.disease = n())
filtered <- filtered %>% group_by(study) %>% mutate(total.by.study = n())
filtered$def.short <- sub('-','_',gsub('[a-z]| ','',filtered$definition))
filtered <- filtered[,c('analysis_id','barcode','study','total.by.study','disease','disease_name','total.by.disease','definition','def.short','total.by.definition','group')]
write.table(filtered, file = '~/Projects/toil-rnaseq-20k/data/metadata_filtered/filtered_datasets.txt', quote = F, row.names = F, sep = "\t")


