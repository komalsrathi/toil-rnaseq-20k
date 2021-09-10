###################################
# Author: Komal S Rathi
# Date: 08/20/2018
# Function: Format new CBTTC
###################################

# new clinical data for cbttc
cbttc <- read.csv('~/Projects/toil-rnaseq-20k/data/metadata/cbttc/cbttc_metadata_v2.csv', stringsAsFactors = F)
cbttc$analysis_id.1 <- NULL
cbttc$study <- NULL
cbttc$study.1 <- NULL
cbttc[is.na(cbttc)] <- ''
cbttc$group <- 'tumors'
cbttc$age_normal_tissue_collected_in_years <- ''
cbttc$tags <- ''
cbttc$gender <- tolower(cbttc$gender)
cbttc$race <- tolower(cbttc$race)
cbttc$ethnicity <- tolower(cbttc$ethnicity)
cbttc$definition <- 'Primary Solid Tumor'
cbttc$study <- "CBTTC"

cbttc$disease <- cbttc$disease_name
cbttc$disease[cbttc$disease_name == "High-grade glioma/astrocytoma (WHO grade III/IV)"] <- "HGG"
cbttc$disease[cbttc$disease_name == "Atypical Teratoid Rhabdoid Tumor (ATRT)"] <- "ATRT"
cbttc$disease[cbttc$disease_name == "Dysplasia/Gliosis"] <- "FCD_Gliosis"
cbttc$disease[cbttc$disease_name == "Medulloblastoma"] <- "MB"
cbttc$disease[cbttc$disease_name == "Low-grade glioma/astrocytoma (WHO grade I/II)"] <- "LGG"
cbttc$disease[cbttc$disease_name == "Metastatic secondary tumors"] <- "Mets_Sec"
cbttc$disease[cbttc$disease_name == "Choroid plexus papilloma"] <- "CPP"
cbttc$disease[cbttc$disease_name == "Dysembryoplastic neuroepithelial tumor (DNET)"] <- "DNET"
cbttc$disease[cbttc$disease_name == "Malignant peripheral nerve sheath tumor (MPNST)"] <- "MPNST"
cbttc$disease[cbttc$disease_name == "Neurofibroma/Plexiform"] <- "NF_PF"
cbttc$disease[cbttc$disease_name == "Choroid plexus carcinoma"] <- "CPC"
cbttc$disease[cbttc$disease_name == "Langerhans Cell histiocytosis"] <- "LCH"
cbttc$disease[cbttc$disease_name == "Brainstem glioma- Diffuse intrinsic pontine glioma"] <- "DIPG"
cbttc$disease[cbttc$disease_name == "Primary CNS lymphoma"] <- "PCNSL"
cbttc$disease[cbttc$disease_name == "Non-germinomatous germ cell tumor"] <- "NGGCT"
cbttc$disease[cbttc$disease_name == "Supratentorial or Spinal Cord PNET"] <- "PNET"
cbttc$disease[cbttc$disease_name == "Non-germinomatous germ cell tumor"] <- "NGGCT"
cbttc[is.na(cbttc)]
cbttc$sample_barcode <- cbttc$analysis_id
write.csv(cbttc, file = '~/Projects/toil-rnaseq-20k/data/metadata/cbttc/cbttc_metadata_v2.csv', quote = F, row.names = F)

apply(cbttc[,c(4:16,21)], 2, FUN = function(x) unique(x))
