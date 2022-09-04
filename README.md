# Mosaic-chromosomal-alteration

This repository contains custom scripts used for processing the raw IDAT files from the Illumina GSA and detecting mosaic chromosomal alterations (mCAs). 

## **Overview**

This repository has two directories:

+ `Process_SNP_Genotyped_Array_Data` is to convert the raw intensity data (IDAT) files from the Illumina Infinium Global Screening Array into VCF files. 

+ `Detect_mosaic_chromosomal_alterations` is to detect mosaic chromosomal alterations (mCAs). 

    + Step 1 : Phasing

    + Step 2 : Imputation

    + Step 3 : MoChA Detection

You can find the required UK Biobank data at https://biobank.ndph.ox.ac.uk/ukb/dset.cgi?id=3094.
