#############################################
# Author: Komal Rathi
# Institute: CHOP
# Function: Make plots of all usable samples
# Summarize data that can be used for SU2C
#############################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)

# 60498 Genes, 198619 Transcripts
annotation <- as.data.frame(data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation.txt'))
annotation <- unique(annotation[,c('gene_id','gene_symbol','transcript_id')])

# code table
code <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/tcga_target_codetable.txt', stringsAsFactors = F)
code[1:9,'code'] <- paste0('0',code[1:9,'code'])

# all normals from GTEx
# 31 Tissue types, 7859 Samples
gtex <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/gtex/gtex_metadata.txt', stringsAsFactors = F)
colnames(gtex)[4] <- 'disease'
gtex <- gtex[,c('study','disease','analysis_id')]
gtex <- gtex[-which(gtex$disease=="<not provided>"),]
gtex$definition <- 'Normals'
# write.table(gtex, file = '~/Projects/toil-rnaseq-20k/data/metadata_filtered/gtex/gtex_metadata.txt', quote = F, row.names = F, sep = "\t")
gtex.ct <- plyr::count(gtex, 'disease')
gtex.ct$study <- 'GTEx'
gtex.ct$definition <- 'Normals'
# write.table(gtex.ct, file = '~/Projects/toil-rnaseq-20k/data/metadata/gtex/GTEx_sample_count.txt', quote = F, row.names = F, sep = "\t")

# all primary tumors from TARGET
target <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/target/target_metadata.txt', stringsAsFactors = F)
target.ct <- plyr::count(target, c('disease'))
target.subsets.ct <- plyr::count(target, c('disease','sample_type'))
target.subsets.ct$code <- sub('TARGET_','',target.subsets.ct$sample_type)
target.subsets.ct <- merge(target.subsets.ct, code, by = 'code', all.x = T)
target.subsets.ct$study <- 'TARGET'
target.subsets.ct <- target.subsets.ct[,colnames(gtex.ct)]

# all primary tumors from TCGA 9307
# 10662 samples, 34 Tissue types
tcga <- read.delim('~/Projects/toil-rnaseq-20k/data/metadata/tcga/tcga_metadata.txt', stringsAsFactors = F)
tcga.ct <- plyr::count(tcga, 'disease')
tcga.subsets.ct <- plyr::count(tcga, c('disease','sample_type'))
tcga.subsets.ct <- merge(tcga.subsets.ct, code, by = 'sample_type', all.x = T)
tcga.subsets.ct$study <- 'TCGA'
tcga.subsets.ct <- tcga.subsets.ct[,colnames(gtex.ct)]

# Medullos: 97 Primary + 1 Relapse
medullo <- data.frame(disease = 'Medulloblastoma', study = 'Medulloblastoma', freq = 97, definition = 'Primary Solid Tumor')

# SCLC: 79 Primary
sclc <- data.frame(disease = 'SCLC', study = 'SCLC', freq = 79, definition = 'Primary Solid Tumor')

# rbind
total <- rbind(gtex.ct, target.subsets.ct, tcga.subsets.ct, sclc, medullo)

# total.ct <- total %>% group_by_(c('disease')) %>% summarise(total.by.disease = sum(freq))
# total <- merge(total, total.ct, by = 'disease')
# total.ct <- total %>% group_by_(c('study')) %>% summarise(total.by.study = sum(freq))
# total <- merge(total, total.ct, by = 'study')
# colnames(total)[3] <- 'total'
# total <- total[,c('study','disease','definition', "total.by.study",'total.by.disease','total')]
# total <- total[order(total$study, total$disease, total$definition),]
# write.table(total, file = '~/Projects/toil-rnaseq-20k/data/metadata/datasets_summary.txt', quote = F, sep = "\t", row.names = F)

# usable samples
all <- rbind(tcga, target)
total.sub <- merge(total, unique(all[,c('disease','disease_name')]), by = 'disease', all.x = T)
total.sub$definition[total.sub$definition == "Additional - New Primary"] <- "Primary Solid Tumor"
total.sub$disease_name[total.sub$study=='GTEx'] <- toupper(total.sub$disease[total.sub$study=='GTEx'])
total.sub$disease_name[total.sub$study=='Medulloblastoma'] <- 'MEDULLOBLASTOMA'
total.sub$disease_name[total.sub$study=='SCLC'] <- 'SMALL CELL LUNG CANCER'
total.sub$type <- ifelse(total.sub$study == "GTEx", "Normals", "Tumors")

total.sub <- total.sub %>% filter(freq > 1) %>% filter(!disease %in% c('CNTL','RT','AML-IF'))  %>% filter(!definition %in% "Metastatic")
total.sub <- total.sub[-which(total.sub$type == 'Tumors' & total.sub$definition == 'Solid Tissue Normal'),] # we are going to compare with GTEx
total.sub <- total.sub %>% group_by_(c('disease')) %>% mutate(total.by.disease = sum(freq))
total.sub <- total.sub %>% group_by_(c('study')) %>% mutate(total.by.study = sum(freq))
total.sub <- total.sub %>% mutate(total = freq) %>% select(study, disease, disease_name, definition, total.by.study, total.by.disease, total) %>% as.data.frame
total.sub <- total.sub[order(total.sub$study, total.sub$disease, total.sub$definition),]

write.table(total.sub, file = '~/Projects/toil-rnaseq-20k/data/metadata/datasets_summary_usable_samples.txt', quote = F, sep = "\t", row.names = F)

# plots
total <- total.sub
q <- plot_ly() %>%
  add_pie(data = total[which(total$study=='GTEx'),], labels = ~disease, values = ~total.by.disease,
          name = "GTEx", domain = list(x = c(0, 0.4), y = c(0.4, 1))) %>%
  add_pie(data = total[which(total$study!="GTEx"),], labels = ~disease, values = ~total.by.disease,
          name = "Tumors", domain = list(x = c(0.6, 1), y = c(0.4, 1))) %>%
  add_pie(data = total[which(total$study=="TARGET"),], labels = ~disease, values = ~total.by.disease,
          name = "Target", domain = list(x = c(2, -4), y = c(0, 0.6))) %>%
  layout(title = 'Toil 20k Summary', showlegend = FALSE,
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# barplot
total.sub <- total.sub[order(total.sub$study, total.sub$total.by.disease),]
total.sub$disease <- factor(total.sub$disease, levels = unique(as.character(total.sub$disease)))
pdf(file = '~/Projects/toil-rnaseq-20k/plots/Toil_SampleSummary.pdf', width = 25, height = 15)
p.not <- ggplot(total.sub, aes(x = disease, y = total, fill = study)) + 
  geom_bar(stat = 'identity') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ggtitle('GTEx Samples (n = 18245)') + ylab('Number of Samples')
p.not
dev.off()
p <- ggplotly(p.not)
p
