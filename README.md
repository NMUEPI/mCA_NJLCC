# **Mosaic-chromosomal-alteration**

This repository contains custom scripts used for processing the raw IDAT files from the Illumina GSA and detecting mosaic chromosomal alterations (mCAs). 

## **Overview**

This repository has two directories:

+ `Array_Data_Processing` is to convert the raw intensity data (IDAT) files from the Illumina Infinium Global Screening Array into VCF files. 

+ `mCAs_Detection` is to detect mosaic chromosomal alterations (mCAs). 

    + Step 1 : Phasing

    + Step 2 : Imputation

    + Step 3 : MoChA Detection

The mCA call sets can be obtained from dataset Return 3094 from UKB application 19808 at https://biobank.ndph.ox.ac.uk/ukb/dset.cgi?id=3094.

