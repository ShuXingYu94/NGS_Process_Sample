# Sample code for NGS processing
> This programme takes advantage of Stacks to analyze NGS data.

* [Required environment](#required-environment)
  + [Shell](#shell)
  + [R](#r)
  + [Python](#python)
* [Installing Prerequisites](#installing-prerequisites)
  + [For shell package(s)](#for-shell-package-s-)
  + [For R package(s)](#for-r-package-s-)
  + [For Python package(s)](#for-python-package-s-)
* [Quick Guide](#quick-guide)
  + [1. Working directory](#1-working-directory)
  + [2. Download shell scripts and configuration](#2-download-shell-scripts-and-configuration)
    - [Download scripts](#download-scripts)
    - [Required configurations](#required-configurations)
    - [Optional configurations](#optional-configurations)
  + [3. Executing script](#3-executing-script)
  + [Ex. To rerun the programme](#ex-to-rerun-the-programme)
* [Documentation](#documentation)
* [Interpreting Results](#interpreting-results)
  + [1. Introduction of the folders](#1-introduction-of-the-folders)
  + [2. Statistical results preview](#2-statistical-results-preview)
  + [3. About numerical data](#3-about-numerical-data)
  + [4. About figures](#4-about-figures)
    - [SNPs Distribution Figure](#snps-distribution-figure)
    - [SNP Depth Figure](#snp-depth-figure)




## Required environment

### Shell
- stacks
- vcftools
- samtools
- trimmomatic
- java
- bwa
### R
- CMplot
### Python
- matplotlib
- pandas
- scipy
- numpy
- etc.

## Installing Prerequisites

### For shell package(s)

    Please refer to the homepage of each package for more information.

### For R package(s)

    Run `install.packages("package_needed")` in R.

### For Python package(s)

    Run the following script in shell. This script will download requirements.txt and automatically install.

```
if command -v wget >/dev/null 2>&1; then
  wget -O ./requirements.txt https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/requirements.txt
else
  if command -v curl >/dev/null 2>&1; then
    curl -o ./requirements.txt https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/requirements.txt
  else
    echo "Wget and Curl are not installed. Please install."
    exit 0
  fi
fi

if command -v pip3 >/dev/null 2>&1; then
  pip3 install -r requirements.txt
else
  if command -v pip >/dev/null 2>&1; then
    pip install -r requirements.txt
  else
    echo "Something wrong with pip/pip3. Please install python packages with other software like conda."
    exit 0
  fi
fi

rm ./requirements.txt
```

## Quick Guide

### 1. Working directory
Please start with the following folder structure.
```
work_dir
├── BaseCall
│   ├── sample1_R1_001.fastq.gz
│   ├── sample1_R2_001.fastq.gz
│   ├── sample2_R1_001.fastq.gz
│   └── sample2_R2_001.fastq.gz
├── MIGadapter.fasta → (Optional)
└── popmap.txt → (Optional)
```
Normally, putting NGS data files in `./work_dir/BaseCall` is all you need to do.

>- For trimming of the original data files, an adapter.fasta file is needed. By default, a `MIGadapter.fasta` file will be downloaded.
>
>- In case of multiple population analysis, you can put a `popmap.txt` file in `./work_dir/`. By default, a `popmap.txt` file will be automatically generated with all the samples recognized in the same population.

### 2. Download shell scripts and configuration

#### Download scripts
Download `main.sh` and `statistics.sh` files from the Master branch to your working directory.

#### Required configurations

Input the required configurations in `main.sh` as follows.

    # Required Settings
    
    workdir=/Users/.../work_dir
    mapping_db=/Users/.../reference_genome.fa
    trimmomatic_dir=/Users/.../Trimmomatic-0.39/trimmomatic-0.39.jar
    adapter=${workdir}/MIGadapter.fasta
    popmap_dir=${workdir}/popmap.txt
    
    # Read file names
    
    files="sample1 sample2"

#### Optional configurations
Basically, there's no need to change the optional configurations.

However, in case of analyzing more than 2 samples, it is recommended to set `stacks_r = 0.8` rather than `stacks_r = 1.0` by default.

 → For more information, please refer to the documentation of Trimmomatic and Stacks.

### 3. Executing script
Please execute the following commands in the same directory as `main.sh` and `statistics.sh`,

```
bash main.sh
bash statistics.sh
```

or simply copy all the commands and execute them.

> **Note**
>
> **Currently, copy all the commands and execute is recommended.**


### Ex. To rerun the programme

In case you might want to delete all the generated files to rerun the programme:

```
rm -rf ${workdir}/aligned
rm -rf ${workdir}/log
rm ${workdir}/popmap.txt
rm -rf ${workdir}/stacks
rm -rf ${workdir}/trimmed
rm -rf ${workdir}/statistics
rm -rf ${workdir}/*.py
```

This will remove all the generated and downloaded files and return to the original [folder structure](#1-working-directory).

> **Warning**
> 
> This is not a part of the analysis !
> 
> Check the directory and files you are deleting !

## Documentation

        **...still in progress, please follow the steps in [Quick Guide](#quick-guide).**

## Interpreting Results

### 1. Introduction of the folders
After execution of `main.sh` and `statistics.sh`, the folder should look something like this:

```
work_dir
├── BaseCall
├── aligned
├── log
├── stacks
├── statistics
│   └── Coverage
└── trimmed

# the files inside is not shown
```
Here, the original `fastq.gz` data files are in folder `./BaseCall`. After trimming with given `adapter.txt` file, the reads are saved in folder `./trimmed`. The trimmed files are subsequently mapped to reference genome and saved in `.bam` format in folder `./algined`. Then, Stacks were used to analyze mapped reads, carry out SNPs calling and population analysis. These results are saved in `./stacks`.

Eventually, the statistical analysis were carried out with `statistics.sh` using the above files. The generated result files are kept in folder `./statistics`, which should look like below.
```
statistics
├── Coverage
│   ├── sample1.txt
│   └── sample2.txt
├── SNPs_Distribution_0.svg
├── SNPs_Distribution_10.svg
├── compare_output.csv
├── compare_statistic.csv
├── file_names.txt
└── results.txt
```
Basically, all the statistical results and figures we need is in this folder.

### 2. Statistical results preview
Currently, the following results are available.

- **Data in table 1**
- [x]  Reads Count of Raw Data
- [x]  Reads Count of Trimmed Data
- [x]  Mapped Length of Trimmed Data 
- [x]  Mapped Coverage of Trimmed Data
- [x]  Average Depth of Mapped Reads
- [x]  Consensus Length
- [x]  Consensus Coverage
- **Data in table 2**
- [x]  SNPs Count
- [x]  Average SNP Depth
- [x]  Mapped Length / SNPs
- [x]  Consensus Length / SNPs
- **Figures**
- [x]  SNPs Distribution Based on Chromosomes
    - [x]  Python - SNP distribution
    - [ ]  R - SNP depth and distribution

All the numerical data are saved in `./statistics/results.txt`. It is recommended to open the file with Excel or any software that supports `.tsv` format.

### 3. About numerical data

File `./statistics/results.txt` contains two independent tables.

- Table 1 contains information from **Reads Count of Raw Data** to **Consensus Coverage**.
- Table 2 contains information from **SNPs Count** to **Consensus Length / SNPs**

### 4. About figures

#### SNPs Distribution Figure 

SNPs Distribution Figure shows the distribution of SNPs on the chromosomes using reference genome information.

- File `./statistics/SNPs_Distribution_0.svg` shows distribution of all the SNPs from `./stacks/populations.snps.vcf`
- File `./statistics/SNPs_Distribution_10.svg` shows distribution of SNPs the depth of which are > 10 on both of the samples.

#### SNP Depth Figure

SNP Depth Figure shows the depth of each SNPs on the chromosomes.

        **...still in progress.**


