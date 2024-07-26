# Multimodal Approach using Polygenic Risk Scores for Prostate Cancer: An End-to-End Methodology

## Introduction

This project leverages data from the All of Us Research Program, along with supplementary data from the study _Trans-ancestry genome-wide association meta-analysis of prostate cancer (PCa) identifies new susceptibility loci and informs genetic risk prediction_ by Conti et al. (2021). The project focuses on developing prostate cancer polygenic risk scores (PRSs) using the betas for 269 susceptibility loci. These PRSs are calculated for approximately 26,000 prostate cancer cases and 2,000 controls.The focus of this study is to evaluate the the combined predictive power of social determinants of health data (SDoH) with genetic and health data.

## Overview

The methodology presented in this project serves as an end-to-end solution for generating polygenic risk scores for any disease. By inputting your research cohort data—comprising survey information, genetic data, and relevant health data specific to the disease of interest—you can utilize the provided effect alleles and weights to develop accurate PRSs.

## Key Features

- **Comprehensive Data Integration**: Combines data from large-scale research programs and meta-analyses for robust PRS development.
- **Scalability**: Adaptable to various disease states beyond prostate cancer, with customizable input requirements.
- **Genetic Risk Prediction**: Enhances understanding of genetic predisposition by identifying susceptibility loci and calculating associated risk scores.

This analysis was adapted from a demonstration project in the All of Us Research Program v6, which outlined the development of polygenic risk scores using the program's data. The original project included several steps, such as:

1. **Liftover Weights Files**: Transition weights files from hg19 to hg38.
2. **Build Test and Train Cohorts**: Select ancestry-balanced train and test cohorts, subset to sites in weights files, and create new datasets. 
3. **Repartition**: Reduce cohort dataset into partitions to improve efficiency in later steps.
4. **Score Cohorts**: Score samples using Hail, based on [Hail's polygenic score calculation guide](https://hail.is/docs/0.2/guides/genetics.html#polygenic-score-calculation), with adjustments for multiallelic sites.
5. **Fit Model**: Fit models and evaluate performance.

This work acknowledges the All of Us Research Program for contributing to the methodology used in this analysis.


