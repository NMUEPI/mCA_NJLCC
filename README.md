# **Mosaic-chromosomal-alteration**

This repository contains custom scripts used for processing raw IDAT files from Illumina GSA (Infinium Global Screening Array) and detecting mosaic chromosomal alterations (mCAs). 

## **Overview**

This repository has three directories:

+ `Array_Data_Processing` included the code used for converting the raw intensity data (IDAT) files from the Illumina Infinium Global Screening Array into VCF files.

+ `mCAs_Detection` included the code used for detecting mosaic chromosomal alterations (mCAs).

    + Step1: Phasing pipeline

    + Step2: Imputation pipeline

    + Step3: Chromosomal alterations pipeline

+ `mCAs_Summary` included mCAs summary data and the code used for association analysis in the Nanjing Lung Cancer Cohort (NJLCC) study.

The mCA call sets from the UK Biobank can be obtained from dataset Return 3094 from UKB application 19808 at https://biobank.ndph.ox.ac.uk/ukb/dset.cgi?id=3094.

