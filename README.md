# sumo_analysis
Analysis and figures for SUMO manuscript.

## Set-up instructions 
1. Clone this repository and its submodule
```
git clone --recurse-submodules https://github.com/ratan-lab/sumo_analysis.git
```
2. Install SUMO from command line (please note that the package require python3.6+):
```
python3 -m pip install --upgrade pip
python3 -m pip install python-sumo
```
3. Install packages from R console:
```R
install.packages(c('PMA', 'PINSPlus', 'R.matlab', 'devtools', 'Matrix', 'rticulate', 'cluster', 'survival', 
                   'tidyverse', 'ggsci', 'ggpubr', 'cowplot', 'gridExtra'))

library("devtools")
install_github("danro9685/CIMLR", ref = 'R')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("iClusterPlus", "DESeq2"))
```
4. Download and install source packages of [SNFtool](https://cran.r-project.org/src/contrib/Archive/SNFtool/SNFtool_2.0.3.tar.gz) and [LRAcluster](http://bioinfo.au.tsinghua.edu.cn/member/jgu/lracluster/LRAcluster_1.0.tgz):
```R
install.packages('LRAcluster_1.0.tgz', repos = NULL, type="source")
install.packages('SNFtool_2.0.3.tar.gz', repos = NULL, type="source")
```

## Multi-omic benchmark

1. Download benchmark data from http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html. 
2. Extract all .zip files into the *sumo_analysis/benchmark/data* directory.
3. From *sumo_analysis/benchmark* directory run **run_benchmark.R** script.
4. To create Fig2 run **plot_benchmark.R** script.

## Benchmark evaluation


## Noisy simulation


## Subsampling simulation
