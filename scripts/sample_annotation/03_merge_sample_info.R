####################################################
# Author: Komal S Rathi
# Function: This script will merge 
# MB, SCLC to exising TCGA, TARGET and GTEx metadata
####################################################

# start merging sample information first
# cancer subtype to sample info
allclin <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/metadata_with_cancersubtypes.txt', stringsAsFactors = F)

# add individual infor
gtex <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/gtex_metadata_7863.txt', stringsAsFactors = F)

target <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/target_metadata_721.txt', stringsAsFactors = F)
target.new <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/target_metadata_1351.txt', stringsAsFactors = F)
tcga <- read.delim('~/Projects/toil-rnaseq-20k/data/clinical/tcga_metadata_10662.txt', stringsAsFactors = F)
tumors <- rbind(tcga, target, target.new)
tumors$disease_name <- stringi::stri_trans_totitle(tumors$disease_name)
tumors$disease_subtype <- ''
rm(target, tcga, target.new)

# merge tumors
alldata <- rbind(gtex, tumors)
rm(gtex, tumors)

# SCLC
sclc <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/sclc/SCLC_metadata.txt')
sclc <- sclc[which(sclc$Comment..Sample_source_name. == "tumor"),]
sclc <- unique(sclc[,c('Comment..ENA_RUN.','Comment..Sample_title.','Comment..Sample_title.')])
colnames(sclc) <- c("analysis_id","patient_barcode","sample_barcode")
sclc$group <- "Tumors"
sclc$study <- "SCLC_GSE60052"
sclc$disease <- "SCLC"
sclc$disease_name <- "Small Cell Lung Cancer"
sclc$disease_subtype <- "Small Cell Lung Cancer"
sclc$tissue <- ""
sclc$definition <- "Primary Solid Tumor"
sclc$subtissue <- ""

# merge SCLC
alldata <- rbind(alldata, sclc)
rm(sclc)

# Medulloblastoma
mb <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/medullo/MB_metadata.txt', stringsAsFactors = F)
colnames(mb)[1] <- "patient_barcode"
mb$sample_barcode <- mb$patient_barcode
mb$group <- "Tumors"
mb$study <- "MB_TaylorLab"
mb$disease <- "MB"
mb$disease_name <- "Medulloblastoma"
mb$disease_subtype <- mb$Subgroup
mb$tissue <- ''
mb$subtissue <- ''
mb <- mb[,colnames(alldata)]

# merge MB
alldata <- rbind(alldata, mb)
write.table(alldata, file = '~/Projects/toil-rnaseq-20k/data/clinical/sample_info.txt', quote = F, sep = "\t", row.names = F)
