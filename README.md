# R4scSHARP

## Purpose:
  This R package was created as a supplementary package for a larger cell classification project, scSHARP
  (find that repository here: https://github.com/mperozek11/scSHARP_tool). R4scSHARP contains a number of helper 
  functions, but the main function, run_r4scsharp, is used to run up to five cell classification tools on a
  RNA gene expression matrix, then return the tool predictions in a succinct file.
  
## Installation:
  For the time being, this package is only available on GitHub. It can be installed by calling the following
  in an R command line or script.
  ```
  library(devtools)
  install_github("W-Holtz/R4scSHARP")
  library(R4scSHARP)
  ```
  We have plans to publish to CRAN soon (a week or so), so stay tuned!

### Installing Suggested Packages:
  Our package contains 4 suggested dependencies that need to be installed on their own.

  #### SingleR:
  Install this package with the command:
  ```
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("SingleR")
  ```

  #### scSorter:
  Install this package with the command:
  ```
  install.packages('scSorter')
  ```

  #### scPred:
  Install this package with the command:
  ```
  devtools::install_github("immunogenomics/harmony")
  devtools::install_github("powellgenomicslab/scPred")
  ```

  #### SCINA
  ```
  install.packages('SCINA')
  ```
  
## Usage:
  More documentation is coming! For now, we've included an examples file containing an R script showcasing
  the basic functionality on a simulated dataset.
  
### Contact:
  Emails: 
  willholtz2001@gmail.com,
  d_lewinsohn@coloradocollege.edu 
