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
2. Download and extract the data provided in [Inputs.zip](https://drive.google.com/file/d/1YFOwEDFLpyxG6WX243FNA7lnrVSaiVL-/view?usp=sharing)
3. Set the working directory to Inputs.
4. Run one of the 6 code modules in R.

## General notes

Please note results of modules 1-3 might slightly differ dependending on the version of R/R packages used.

Each code module can be run independently.

The code uploaded here is the working code being used throughout the work on the project. 

The MP generation approach is based on our earlier work described in Kinker et al. 2020 and Gavish et al. 2023.  Further documentation can be found [here](https://github.com/gabrielakinker/CCLE_heterogeneity) and [here](https://github.com/tiroshlab/3ca).

We have updated Module 2 to include the code used to generate Figures 4A-C and clarified which code was used for different panels.   

The Inputs file contains an additional README regarding alignment of Visium and CODEX samples. 

You can find the Visium H&E [here](https://drive.google.com/file/d/19ROs5wKtDH-RLqELFFBeEexKScW20FJL/view?usp=drive_link). 

Requirements
R (tested in version 4.1.1).
