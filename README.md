# Ofav_SHratiosNH4

* R_Script_Ofav-SH-NH4.R: Script to clean data and run preliminary analysis. It uses:
    - data: folder that contains the raw qPCR data in .csv files
    - STEPoneFunction.R: Ross Cunning function to calculate the S/H cell ratios using uotput files from StepOne machines
    - TGFB-Final_AP: Contains information about each O.fav core, as treatments and notes
    - ToRem1: File that contains the duplicated samples I decided to remove from the analysis
    
* Ofav-SH-NH4_Notes.html and .Rmd ducument the code in R_Script_Ofav-SH-NH4.R and show the results

* DuplicatesA.csv: is generated each time that R_Script_Ofav-SH-NH4.R is run and containes teh samples that were ran multiple times.    
