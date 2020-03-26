# **Developing and testing an elastic-net logistic regression immune cell classifier**

This repository supplies the code developed in the study of Torang _et al._ **_"An elastic-net logistic regression approach to generate classifiers and gene signatures for types of immune cells and T helper cell subsets"_**. The corresponding pre-print can be found on BMC Bioinformatics ( https://doi.org/10.1186/s12859-019-2994-z ). It can be used to reproduce the results of the study and investigate the methodology to be used for other datasets.

## **Requirements**

* R version 3.5.3 -- "Joy in Playing".
* R libraries: dplyr, annotables, EDASeq, glmnet, DescTools, matlab, tidyverse, plotrix, MASS, plyr, rafalib, factoextra, NbClust, Matrix, pROC, xlsx, MESS, ggplot2

## **Data**

All necessary data is provided in the "Files" folder of the repository or can be found in Gene Expression Omnibus repository [https://www.ncbi.nlm.nih.gov] with the following GEO accession numbers: GSE60424, GSE64655, GSE36952, GSE84697, GSE74246, GSE70106, GSE55536, GSE71645, GSE66261, GSE96538, GSE75688, GSE72056. Sourcing for all data can be found in Torang _et al._

## **Quick start**

To reproduce the results, download the relevant script and load the corresponding data into the R workspace. Running the script will generate all relevant figures and data tables for the given portion of the study.

# General notes

The code provided in this repository reproduces the main results of the study of Torang _et al._ **_"An elastic-net logistic regression approach to generate classifiers and gene signatures for types of immune cells and T helper cell subsets"_** but it is not meant as a self-contained module for use in further analysis.

## Citation
To download the citation please go to  https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2994-z.ris

Torang, A., Gupta, P. & Klinke, D.J. **_"An elastic-net logistic regression approach to generate classifiers and gene signatures for types of immune cells and T helper cell subsets." BMC Bioinformatics 20, 433 (2019)."_** https://doi.org/10.1186/s12859-019-2994-z
