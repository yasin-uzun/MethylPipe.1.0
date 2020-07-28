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
If MethylPipe is installed without any problems, you should see the following message:
> \[1\] "MethylPipe works ok."

## Dependencies 

To run MethylPipe, you need to have the underlying software:
* Adapter Trimmer: Cutadapt or 
* Aligner: Bismark (with Bowtie) or BSMAP or BS3
* Duplicate removal: samtools or picard
* Demultiplexer: demultiplex_fastq.pl perl script (see below).

You also need genomic sequence and annotated genomic regions for quantification of methylation calls. Note that you only need the tools you will use to be installed, i.e, you don't need BSMAP or BS3 if you will only use Bismark as the aligner.

You can install these tools by yourself. For convenience, we provide the binaries in [here](https://chopri.box.com/s/l8o4v6ko8aeabo3fsdtfan8gxjxzg39h) . Please cite the specific tool when you use it, in adition to MethylPipe.

You can download demultiplex_fastq.pl from here[https://chopri.box.com/s/vplpxht3r7u6i0fcnio803wlnezuc5o3] .



