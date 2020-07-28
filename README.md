# MethylPipe.1.0 Single Cell DNA Methylation Data Processing Pipeline

MethylPipe is an R package for processing single cell DNA methylation data. It accepts fastq files as input, performs demultiplexing, adapter trimmming, mapping, quantification, dimensionality reduction and differential methylation analysis for single cell DNA methylation datasets.

## Installation
To install MAPLE, type the following commands in R command prompt:
```R
library(devtools)
install_github("yasin-uzun/MethylPipe.1.0")
```
Once you have installed the MethylPipe package, verify that it is installed correctly as follows:

```R
library(MethylPipe)
test()
```
If MethylPipe is installed without problems, you should see the following message:
> \[1\] "MethylPipe works ok."


