# Colorectal Cancer Meta-Analysis

TODO

## Requirements

- This project requires R version __3.5.1__ -- "Feather Spray"

- Please make sure that you have all required packages installed. You can
easily check whether that is the case by running:

    ```bash
    Rscript requirements.R
    ```
- This project relies heavily on the R package `SIAMCAT`.  
    For the case that you are using an older version of the package, it will
    undoubtedly throw some errors. Thus, please make sure that you have the
    version __1.1.0__ installed. You can install the correct version of  
    `SIAMCAT` via the
    [gitlab repository](https://git.embl.de/grp-zeller/SIAMCAT/tags/v1.1.0).
    Download the tar-ball of the package and then install the package by typing
    in R:

    ```R
    install.packages("SIAMCAT-v1.1.0.tar",
        repos = NULL, type="source")
    ```
    - _you may have to change the path to the path of the downloaded file_
    - _in Windows, you will need to have the package
    [rtools](https://cran.r-project.org/bin/windows/Rtools/) installed_
        - _i guess... never tested on a Windows system_

## Workflow

To reproduce the results, please follow these instructions.

### Setup

First, you will have to download the taxonomic and functional profiles and the
metadata needed for this project. They are stored on the ftp server of the
Zeller team and will be automatically downloaded and cleaned by running
this script.

```bash
cd ./src
Rscript prepare_data.R
```
