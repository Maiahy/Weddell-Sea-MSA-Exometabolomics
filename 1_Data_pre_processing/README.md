# Data preprocessing

This folder contains the R scripts used for data cleaning, quality control, 
normalization, blank correction, and averaging of technical and environmental 
replicates prior to statistical analysis:

1. Exploration of sample and pooled-QC intensities
2. Assessment of analytical reproducibility using QC injections (RMAD + RSD)
3. Normalization according to the volume of seawater extracted
4. Assessment of procedural blank contribution
5. Averaging of injection replicates
6. Correction of procedural blank contribution
7. Averaging of onboard technical replicates
8. Export of the preprocessed dataset


This code requires as input the *SKYLINE_Peak_areas.txt* feature table and the 
*Metadata_final_with_operators.txt* metadata table, both available on MetaboLights 
among the supplementary files under the study accession **MTBLS15522**.


>***NOTE*** 
Code written with R version 4.5.2

