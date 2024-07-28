# scRNAseq-and-Nanostring-Cosmx-SMI



## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
- [Issues](https://github.com/ebridge2/lol/issues)

# Overview
Single cell study on mouse melanoma tumour infiltrating immune cells and human melanoma peripheral blood mononuclear cells from patients post anti-PD1 therapy.Spatial transcriptomics analysis of mouse tumor tissue from Rag-/- B6 mice. Mice were reconstituted with B16 tumor and then at day 7 reconstituted with WT T effectors along with either WT T regulatory cells or PD1-/- T regulatory cells. Nanostring spatial analysis on tumor tissue was performed when tumor size reached >800mm3.

# Repo Contents

- [R](./R): `R` package code.
- [tests](./tests): A URL to the BD Rhapsody Github bootcamp training.

# System Requirements

## Hardware Requirements

The `R` package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3.3 GHz) and internet of speed 25 Mbps.

## Software Requirements

### OS Requirements

The R scripts were tested on *Mac OSX* operating system. The developmental version of the package has been tested on the following systems:
  
Mac OSX:  
Windows:  

Before running the listed 'R' scripts, users should have `R` version 4.4.1 or higher, and several packages set up from CRAN.

#### Installing R version 4.4.1 on Mac OSX
Installation of R requires Mac OSX with the copy of XQuartz which should install in about 2 hours.
The latest version of R can be installed by adding the latest repository to `terminal`
which should install in about 5 minutes.

# Installation Guide

## Package Installation

Users should install the following packages prior to running the script, from an `R` terminal:

```
install.packages("ggplot2, tidyverse, "Matrix", "RCurl", "scales", "data.table", "readxl", "BiocManager", "ggpubr", "Seurat")
BiocManager::install("ensembldb", "org.Hs.eg.db", "clusterProfiler", "biomaRt", "enrichplot", "AnnotationHub")
```

Each of which will install in about 1 minute on a machine with the recommended specs.

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/JXLIM99/scRNAseq-and-nanostring-Cosmx-SMI/issues). 

# Demo 
Demo files include a test run of BD Rhapsody single-cell RNA analysis bootcamp

# Results
In this study, we analysed the effect of blocking PD-1 on inhibitory coreceptor expression in T regulatory cells found within melanoma samples.
