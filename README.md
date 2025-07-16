# Code for: "Interaction legacies and functional redistribution restructure the FPOM and DOM transformation process in aquatic ecosystems"
This repository contains code used for statistical analyses and figure generation in the above manuscript.

## Data availability
All datasets used in the analysis are available on both:
- GitHub repository: https://github.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025
- Zenodo (DOI): https://doi.org/10.5281/zenodo.15959018
  
The `read.csv()` calls in the analysis scripts refer to the GitHub raw files, which are archived under the same Zenodo DOI for reproducibility.

## Raw sequence data
The raw amplicon sequencing data have been deposited in the DNA Data Bank of Japan (DDBJ) under the BioProject accession number: **PRJDB18899**.

You can access the raw FASTQ files at:
[https://ddbj.nig.ac.jp/resource/bioproject/PRJDB18899](https://ddbj.nig.ac.jp/resource/bioproject/PRJDB18899)

Sample metadata and mapping files used for analysis are available in the `data/` folder of this repository.

## Folder structure
- `scripts/`: R scripts used for data processing, modeling, and visualization
- `data/`: Sample input data
- `scripts/setup.R`: Installs all required R packages

## Requirements
- R version 4.3.0 or later
- Required R packages are listed and installed via scripts/setup.R

## How to run
1. Clone this repository
Download or clone the repository to your local environment:
```
git clone https://github.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025.git 
cd Interaction-legacies-and-functional-redistribution-2025
```

3. Install R and required packages
Ensure you have R (version ≥ 4.3.0) installed.
Then, run the following in R or RStudio to install all required packages (including CRAN and Bioconductor packages):
```
source("scripts/setup.R")
```

3. Run analysis scripts
You can execute each script in the scripts/ folder depending on your interest (e.g., dbRDA, diversity metrics, pathway analysis).
For example:
```
source("scripts/Main/01.fpom_analysis.R")
```

## Output
All plots and tables will be displayed directly in the R session.
No files will be saved to disk by default.
If you wish to export figures, please modify the relevant scripts to include ggsave() or similar output functions.


## Preprocessing (FASTQ to ASV CSV)
Raw FASTQ files were processed using QIIME2 and DADA2 to generate ASV and taxonomy tables.  
All commands used for QIIME2 preprocessing are described in [`scripts/qiime_dada2_processing.md`](scripts/qiime_dada2_processing.md).


---


### License
This repository is shared under the MIT License.
