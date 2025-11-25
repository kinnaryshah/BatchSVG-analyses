# BatchSVG-analyses

This repository contains the analyses for the BatchSVG project.

There are two folders for the two datasets used for the project: `humanHippocampus2024/` and `spatialLIBD_DLPFC_12_3_7_12_expanded/`. Each of those two folders contains a `run_analysis.nf` Nextflow file detailing the analysis. The scripts used in the Nextflow file are located in the `scripts/spatial/` folder. Each folder also contains `plots/` and `results/` folders for the intermediate results.

The files `HPC_analysis.R` and `DLPFC_analysis.R` contain code to recreate the figures in the manuscript.
