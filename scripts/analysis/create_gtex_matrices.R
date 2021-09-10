# Author: Komal S. Rathi
# Date: 10/08/2018
# Function: Create matrices by tissues

library(data.table)
library(reshape2)
library(limma)
library(dplyr)
library(RDiseaseXpress)

# Step1: Create matrix
# gtex tissues
gtex.studies <- getSamples(myStudy = "GTEx")
gtex.studies <- gtex.studies[,c('sample_id','tissue')]
studies <- unique(gtex.studies$tissue)
studies <- gsub(' ','_',studies)

# annotation
ann <- data.table::fread('/mnt/toil_20k/data/annotation/gencode.v23.annotation.gi_ti_gs.txt')
ann$transcript_id <- NULL
ann <- unique(ann)
ann <- data.table(ann, key = 'gene_id')

# do gsva once per tissue
for(i in 1:length(studies))
{
  print(paste0("Reading...", studies[i]))
  fname <- paste0('/mnt/toil_20k/data/split_files/gtex_wtt/', studies[i], '_wtt.csv')
  gtex.study <- fread(fname)
  colnames(gtex.study) <- c("gene_id", "transcript_id.s.", "length", "effective_length", "expected_count", "TPM", "FPKM", "id", "tissue")
  gtex.study <- dcast(gtex.study, gene_id~id, value.var = 'FPKM')
  gtex.study <- data.table(gtex.study, key = 'gene_id')
  
  # merge res with ann
  gtex.study <- merge(ann, gtex.study)
  gtex.study$gene_id <- NULL
  gtex.study <- as.data.frame(gtex.study)
  
  # collapse by symbol
  gtex.study <- gtex.study %>% 
    unique() %>% 
    dplyr::mutate(means = rowMeans(.[2:ncol(gtex.study)])) %>% 
    arrange(desc(means)) %>% 
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    ungroup() %>% 
    dplyr::select(-means)
  
  # write down matrix
  write.table(gtex.study, file  = paste0('/mnt/toil_20k/data/gtex/',studies[i],'_FPKM_matrix.txt'), quote = F, sep = "\t", row.names = F)
  print(paste0("Done...", studies[i]))
}
