setwd('~/Projects/toil-rnaseq-20k/data/metadata/')

library(RDiseaseXpress)
library(tidyverse)

# manifest from s3
tcga_manifest_from_s3 <- read.delim('tcga/tcga-manifest', stringsAsFactors = F, header = F, col.names = c("analysis_id"))
tcga_manifest_from_s3$analysis_id <- gsub('s3://cgl-rnaseq-recompute-fixed/tcga/|.tar.gz','',tcga_manifest_from_s3$analysis_id)
tcga_manifest_from_s3$type <- gsub('[.].*','',tcga_manifest_from_s3$analysis_id)
tcga_manifest_from_s3$type[!tcga_manifest_from_s3$type %in% c('SINGLE-END','IMPROPERLY_PAIRED')] <- "USABLE"
tcga_manifest_from_s3$analysis_id <- gsub("IMPROPERLY_PAIRED.|SINGLE-END.", "", tcga_manifest_from_s3$analysis_id)

# now add tcga barcodes + disease
tcga_file_one <- read.delim('tcga/unaligned.tsv')
tcga_file_one <- tcga_file_one %>%
  filter(analysis_id %in% tcga_manifest_from_s3$analysis_id) %>%
  dplyr::select(barcode, analysis_id, disease)
tcga_file_two <- read.delim('target/LATEST_MANIFEST.tsv')
tcga_file_two <- tcga_file_two %>%
  filter(analysis_id %in% tcga_manifest_from_s3$analysis_id) %>%
  dplyr::select(barcode, analysis_id, disease)
tcga_mapping_file <- unique(rbind(tcga_file_one, tcga_file_two))

# add back to manifest
final <- tcga_mapping_file %>%
  inner_join(tcga_manifest_from_s3, by = c("analysis_id"))
final$type[final$disease == "CNTL"] <- "CONTROL"

# which sample is missing from DE?
disease_express <- getSamples(myStudy = "TCGA")
setdiff(disease_express$sample_barcode, final$barcode)
final %>%
  filter(barcode %in% setdiff(final$barcode, disease_express$sample_barcode)) %>%
  group_by(type, disease) %>%
  summarise(n())

# looks all good!