# Cis-regulatory element pipeline for Heliconius

## Overview

This pipeline peforms a large scale annotation and characterisation of cis-regulatory promoters and enhancers across spatiotemporal conditions in Heliconius, alongside using sequence characteristics derived from cis-regulatory enhancers to identify enhancers in additional Heliconius species

### Project integrates:
- ATAC-seq (chromatin accessibility)
- Histone modification data (H3K27ac, H3K4me3)
- Motif and transcription factor analysis
- Differential accessibility/enrichment testing
- Grelu dilated convolutional neural network for  regulatory sequence prediction

---

## Cis-regulatory characterisation pipeline

Each stage of the main cis regulatory pipeline is modular and controlled via binary flags (0 = off / 1 = on)

## Pipeline structure

### 1. Data preprocessing and quality control
Handles raw data acquisition, mapping, and QC

#### Modules
- DOWNLOADFASTQ – Retrieve sequencing data from SRR accessions
- RUNMAKEINDEX – Build genome indices for alignment
- RUNBUSCO – Assess genome completeness (BUSCO/compleasm)
- RUNPREPROCESSING – Read trimming, alignment, filtering
- RUNCHECK – Replicate concordance checks
- RUNFRAGMENT – ATAC fragment size distribution analysis
- RUNPEAKCALLING – Peak calling using MACS3
- RUNQUALITYCHECK – Peak QC using deeptools

---

### 2. Cis regulatory element characterisation
Defines promoters and enhancers and links them to target genes

#### Modules
- RUNCONSENSUS - Generates consensus peak sets across datasets (ATAC, H3K27ac, H3K4me3)
- RUNINTERSECT - Defines CRE classes (Promoters: ATAC + H3K27ac + H3K4me3 + within 1kb of TSS / Enhancers: ATAC + H3K27ac + absence of H3K4me3 + outside 1kb of TSS)
- RUNREGMAP - Maps CREs to genes using cis-regulatory mapping
- RUNGENEANNOTATION - Assigns CREs to known wing-development genes via orthologous RBH search

---

### 3. Transcription factor and motif analysis
Identifies regulatory motifs and infers active transcription factors

#### Modules
- RUNMOTIF – Motif discovery (de novo + known) using HOMER
- RUNGC – GC content analysis of CREs vs genomic background
- RUNCLUSTER – Motif clustering (RSAT-based matrix clustering)
- RUNTOBIAS - Performs Tn5 bias correction, footprinting, and transcription factor activity inference via BINDetect)
- RUNTFMATCH - Motif to transcription factor mapping via TOMTOM (uses JASPAR 2024 insect database)

---

### 4. Differential analysis with DESeq2
Quantifies changes in regulatory activity across spatiotemporal conditions

#### Modules
- RUNDESEQ2 - Generates featureCount matrices, builds sample metadata, performs differential analysis of ATAC accessibility and histone enrichment

---

## Regulatory sequence prediction pipeline

## Pipeline structure

### 1. 01_prepare_labels
Prepare CNN labels and split into train/validation/test datasets

- matchbackground.py
- splitbeds.py

### 2. 02_train_model
Label setup and script to train gReLU CNN models

- prepare_grelu_data.py
- tune_grelu.py
- train_grelu.py

### 3. 03_predict_enhancers
Predict enhancer regions using gReLU CNN model and run ATAC overlap permutation assessment

- predict_grelu.py

### 4. 04_CNEE_overlap
Intersect predicted enhancers with conserved CNEEs
