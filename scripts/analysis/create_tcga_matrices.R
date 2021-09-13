library(RDiseaseXpress)

# unaligned 
tcga_dir <- file.path('~/Projects/toil-rnaseq-20k/', 'data', 'metadata', 'tcga')
tcga_manifest_from_s3 <- read.delim(file.path(tcga_dir, 'tcga-manifest'), stringsAsFactors = F, header = F)
tcga_manifest_from_s3$V1 <- gsub('s3://cgl-rnaseq-recompute-fixed/tcga/|.tar.gz','',tcga_manifest_from_s3$V1)
tcga_manifest_from_s3$V2 <- gsub('[.].*','',tcga_manifest_from_s3$V1)
tcga_manifest_from_s3 <- tcga_manifest_from_s3[-which(tcga_manifest_from_s3$V2 %in% c('SINGLE-END','IMPROPERLY_PAIRED')),]

# from diseaseXpress
tcga_meta <- getSamples(myStudy = "TCGA") %>%
  dplyr::select(sample_barcode, sample_id)
set.seed(100)
tcga_meta <- tcga_meta %>%
  group_by(sample_barcode) %>%
  arrange(sample_id) %>%
  distinct(sample_barcode, .keep_all = T)
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
write.table(colnames(tcga_counts), file = '~/Desktop/TCGA_QC/TCGA_usable_barcode_list.txt', quote = F, row.names = F, col.names = F)

# 266 are in OpenPedCan but not in expression data
x <- setdiff(colnames(tcga_counts), histology$Kids_First_Biospecimen_ID)
write.table(x, file = '~/Desktop/TCGA_QC/TCGA_in_diseaseXpress_not_in_OpenPedCan.txt', quote = F, row.names = F, col.names = F)

# 681 are in histology but not in OpenPedCan
y <- setdiff(histology$Kids_First_Biospecimen_ID, colnames(tcga_counts))
write.table(y, file = '~/Desktop/TCGA_QC/TCGA_in_OpenPedCan_not_in_diseaseXpress.txt', quote = F, row.names = F, col.names = F)

