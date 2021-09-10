.. |date| date::

*******************
RNASeq 20k Analysis
*******************

:authors: Komal Rathi, Anthony Cros, Karthik Kalletla, Pichai Raman
:contact: rathik@email.chop.edu
:organization: DBHi, CHOP
:status: This is "work in progress"
:date: |date|

.. meta::
   :keywords: rnaseq, 20k, tcga, gtex, pnoc, target, 2018
   :description: RNAseq 20k samples.

Introduction
============

RNA-seq data Pre-processing, Processing and Analysis for GTEx, TCGA, TARGET, PNOC, CBTTC and other pediatric datasets.

Task List
=========

Processing:
"""""""""""

* Set up an EC2 instance with s3cmd. Reference: docs/getdata.rst
* Get and Put data from and to an S3 bucket. Reference: scripts/get-data.sh.
* Get metadata/sample info for each dataset.
* Concatenate - group by study and normalization type into 12 datasets. Reference: scripts/collect-data.sh and scripts/merge-data.R.
* Sort datasets by gene ID.
* For each dataset, split files by disease/tissue type.

Bioinformatics Analysis:
""""""""""""""""""""""""

* Gene and Isoform level Tumor vs Normal analysis - Primary and Relapse. Reference: scripts/tumorvsnormal_diffexpr_genes.R and scripts/tumorvsnormal_diffexpr_isoforms_64GB.R
* Cox Regression/Neuroblastoma Survival analysis.
* Tumor-Normal clustering.
* Use CD19 and HER2 as reference genes.
  
`Data Summary`_
"""""""""""""""

* Total samples: 19418; Usable samples: 18245.
* Adult Tumors: n = 9711 (53.22%)

    - TCGA: n = 9535
    - Medulloblastomas: n = 97
    - SCLC: n = 79
    - 99.52% Primary, 0.48% Recurrent

* Pediatric Tumors: n = 675 (3.70%)

    - TARGET: n = 675
    - 82.51% Primary, 17.48% Recurrent

* Normals: n = 7859 (43.07%)

    - GTEx: n = 7859
  
.. note::

	Metastatic samples are only available for TCGA and have a small sample size, usually n = 1 or 2. BRCA = 7, TCHA = 8, SKCM = 366. So we are not using them in the analysis.

.. _getdata: ./docs/getdata.rst
.. _Data Summary: ./data/metadata_filtered/filtered_datasets.txt