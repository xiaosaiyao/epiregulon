![plot](inst/epiregulon_logo.png)<br>

# Introduction

![epiregulon_logo-01](https://github.com/xiaosaiyao/epiregulon/assets/70164860/67b4ce37-0cf7-4832-8343-4dfd9c96892a)


Gene regulatory networks model the underlying gene regulation hierarchies that drive gene expression and cell states. The main function of the epiregulon package is to construct gene regulatory networks and infer transcription factor (TF) activity in single cells by integration of scATAC-seq and scRNA-seq data and incorporation of public bulk TF ChIP-seq data.

For full documentation, please refer to the epiregulon [book](https://xiaosaiyao.github.io/epiregulon.book/).

![plot](inst/epiregulon_schematics.svg)   
There are three related packages. The core epiregulon package supports input in the form of `SingleCellExperiment` objects. If the users would like to start from `ArchR` projects, they may choose to use `epiregulon.archr` package, which allows for the seamless integration with `ArchR` package through accepting its output to be used in the downstream workflow. Moreover, we provide a suite of tools in `epiregulon.extra` package for the enrichment analysis, visualization, and network analysis which can be run on the `epireglon` or `epiregulon.archr` output.

# Installation

```
# install devtools
if(!require(devtools)) install.packages("devtools")

# install basic epiregulon package
devtools::install_github(repo='xiaosaiyao/epiregulon')

# install extended version of epiregulon
devtools::install_github(repo='xiaosaiyao/epiregulon.archr')

# install extended version of epiregulon
devtools::install_github(repo='xiaosaiyao/epiregulon.extra')
```

Example data included in the tutorial are available from [scMultiome](https://bioconductor.org/packages/release/data/experiment/html/scMultiome.html) 

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scMultiome")
```

# Functions
Functions in the suite of Epiregulon packages
![plot](inst/epiregulon_functions.png)


Contact: [Xiaosai Yao](mailto:yao.xiaosai@gene.com), Genentech Inc.


