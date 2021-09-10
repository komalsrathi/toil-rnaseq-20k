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

.. code-block:: R

    # sample count by study

           study_id   n()
    1         CBTTC  1110
    2          GTEx  7863
    3          PNOC    27
    4 SCLC_GSE60052    79
    5        TARGET  1399
    6          TCGA 10661  

    # sample count by study + definition

            study_id                                              definition  n()
    1          CBTTC                                              Metastatic   26
    2          CBTTC                                     Primary Solid Tumor 1048
    3          CBTTC                                   Recurrent Solid Tumor   18
    4          CBTTC                                             unavailable   18
    5           GTEx                                                 Normals 7863
    6           PNOC                                              Metastatic    2
    7           PNOC                                     Primary Solid Tumor   25
    8  SCLC_GSE60052                                     Primary Solid Tumor   79
    9         TARGET      Blood Derived Cancer - Bone Marrow, Post-treatment   12
    10        TARGET Blood Derived Cancer - Peripheral Blood, Post-treatment    1
    11        TARGET                                              Metastatic    1
    12        TARGET              Primary Blood Derived Cancer - Bone Marrow  672
    13        TARGET         Primary Blood Derived Cancer - Peripheral Blood  129
    14        TARGET                                     Primary Solid Tumor  449
    15        TARGET            Recurrent Blood Derived Cancer - Bone Marrow  108
    16        TARGET       Recurrent Blood Derived Cancer - Peripheral Blood    3
    17        TARGET                                   Recurrent Solid Tumor   13
    18        TARGET                                     Solid Tissue Normal   11
    19          TCGA                                Additional - New Primary   11
    20          TCGA                                   Additional Metastatic    1
    21          TCGA                                              Metastatic  392
    22          TCGA         Primary Blood Derived Cancer - Peripheral Blood  172
    23          TCGA                                     Primary solid Tumor 9289
    24          TCGA                                     Primary Solid Tumor   18
    25          TCGA                                   Recurrent Solid Tumor   50
    26          TCGA                                     Solid Tissue Normal  728  
  
.. _getdata: ./docs/getdata.rst
.. _Data Summary: ./data/metadata_filtered/filtered_datasets.txt