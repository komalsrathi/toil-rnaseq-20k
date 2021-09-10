library(RDiseaseXpress)

# unaligned 
tcga_manifest_from_s3 <- read.delim('tcga/tcga-manifest', stringsAsFactors = F, header = F)
tcga_manifest_from_s3$V1 <- gsub('s3://cgl-rnaseq-recompute-fixed/tcga/|.tar.gz','',tcga_manifest_from_s3$V1)
tcga_manifest_from_s3$V2 <- gsub('[.].*','',tcga_manifest_from_s3$V1)
tcga_manifest_from_s3 <- tcga_manifest_from_s3[-which(tcga_manifest_from_s3$V2 %in% c('SINGLE-END','IMPROPERLY_PAIRED')),]

# from diseaseXpress
tcga_meta <- getSamples(myStudy = "TCGA")
tcga_meta <- tcga_meta %>%
  dplyr::select(sample_barcode, sample_id)
setdiff(tcga_manifest_from_s3$V2, tcga_meta$sample_id)

# 2 samples are missing, one is Control Analyte and the other one has no documentation so it was not included

# histology file from OpenPedCan repo
histology <- read.delim('~/Projects/PediatricOpenTargets/OpenPedCan-analysis/data/histologies.tsv')
histology <- histology %>%
  filter(cohort == "TCGA")

# tpm from diseaseXpress
tcga_tpm <- readRDS('~/Projects/toil-rnaseq-20k/data/results/expression/tumor_datasets_TPM.RDS')
tcga_tpm <- tcga_tpm %>%
  dplyr::select(tcga_meta$sample_id)
identical(colnames(tcga_tpm), tcga_meta$sample_id)
colnames(tcga_tpm) <- tcga_meta$sample_barcode
saveRDS(tcga_tpm, '~/Projects/toil-rnaseq-20k/data/results/expression/tcga-gene-expression-rsem-tpm-collapsed.rds')

# counts from diseaseXpress
tcga_counts <- readRDS('~/Projects/toil-rnaseq-20k/data/results/expression/tcga_target_counts.rds')
tcga_counts <- tcga_counts %>%
  dplyr::select(tcga_meta$sample_id)
identical(colnames(tcga_counts), tcga_meta$sample_id)
colnames(tcga_counts) <- tcga_meta$sample_barcode
saveRDS(tcga_counts, '~/Projects/toil-rnaseq-20k/data/results/expression/tcga-gene-counts-rsem-expected_count-collapsed.rds')

length(setdiff(colnames(de_tcga_exp), histology$Kids_First_Biospecimen_ID))
# 266 are in OpenPedCan but not in expression data

length(setdiff(histology$Kids_First_Biospecimen_ID, colnames(de_tcga_exp)))
# 681 are in histology but not in OpenPedCan


