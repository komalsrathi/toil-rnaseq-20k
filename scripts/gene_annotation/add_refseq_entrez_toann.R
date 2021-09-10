############################################
# Author: Komal S Rathi
# Institute: CHOP
# Function: Add Entrez/Refseq ids to Gencode
############################################

setwd('~/Projects/toil-rnaseq-20k/')

library(plyr)
library(data.table)

# annotation file
ann <- data.table::fread('data/annotation/superset/gencode.v23.chr_patch_hapl_scaff.annotation.txt', stringsAsFactors = F)
ann <- unique(ann)

# entrez/transcript_id mapping
entrez <- data.table::fread('data/annotation/gencode.v23.metadata.EntrezGene.txt', header = F, stringsAsFactors = F)
colnames(entrez) <- c('transcript_id','entrez_id')
entrez <- unique(entrez)

# merge
ann.sub <- merge(ann, entrez, by = 'transcript_id', all.x = TRUE)
ann.sub <- ann.sub %>% group_by(transcript_id, chr, gene_start, gene_end, gene_strand, gene_id, gene_biotype, gene_symbol, transcript_start, transcript_end, transcript_strand, transcript_biotype) %>%
  summarise(entrez_id = toString(entrez_id)) %>%
  as.data.frame()

# refseq/transcript id mapping
refseq <- data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.metadata.RefSeq.txt', header = F, stringsAsFactors = F)
colnames(refseq) <- c('transcript_id','refseq_mrna_id','refseq_protein_id')
refseq <- unique(refseq)

# merge
ann.all <- merge(ann.sub, refseq, by = 'transcript_id', all.x = TRUE)
ann.all <- ann.all %>% group_by(transcript_id, chr, gene_start, gene_end, gene_strand, gene_id, gene_biotype, gene_symbol, transcript_start, transcript_end, transcript_strand, transcript_biotype, entrez_id) %>%
  summarise(refseq_mrna_id = toString(refseq_mrna_id),
            refseq_protein_id = toString(refseq_protein_id)) %>%
  as.data.frame()
final <- ann.all[,c('chr','gene_symbol',
           'gene_id','gene_start','gene_end','gene_strand','gene_biotype',
           'transcript_id','transcript_start','transcript_end','transcript_strand','transcript_biotype',
           'entrez_id','refseq_mrna_id','refseq_protein_id')]
write.table(final, file = '~/Projects/toil-rnaseq-20k/data/annotation/superset/gencode.v23.chr_patch_hapl_scaff.annotation_otherids.txt', quote = F, sep = "\t", row.names = F)

