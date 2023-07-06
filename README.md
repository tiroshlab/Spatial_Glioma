# Spatial_Glioma

This resource provides the R code to reproduce key results described in Greenwald, Galili Darnell, Hoefflin, et al. "Integrative spatial analysis reveals a multi-layered organization of glioblastoma".

The analyses are divided into 6 main modules:

Module 1: Per sample clustering, metaprogram generation, and spot annotation

Module 2: Spatial coherence and organizational zones

Module 3: Measures of spatial associations 

Module 4: defining consensus interactions

Module 5: Spatial CNA inference

Module 6: CODEX analysis

## Getting started
1. Clone Github repository.
2. Download and extract the data provided [Inputs.zip](https://drive.google.com/file/d/19YULZZbQpXnLwCnOCos5xm2PRFeHbuCe/view?usp=sharing)
3. Set the working directory to Inputs.
4. Run one of the 6 code modules in R.

## General notes

Please note results of modules 1-3 might slightly differ dependending on the version of R/R packages used.

Each code module can be run independently.

The code uploaded here is the working code being used throughout the work on the project. We are currently working on generating a more user-friendly, readable and easy to run code. 

The following metaprogram naming in metadata and some modules matches a pervious version of the text and will be updated to reflect current naming shortly:
metabolism = Prolif-Metab
inflammatory.resp = Inflammatory-Mac

Requirements
R (tested in version 4.1.1).
