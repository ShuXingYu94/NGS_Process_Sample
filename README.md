# Sample code for NGS processing
> This programme takes advantage of Stacks to analyze NGS BaseCall data.

* [Required environment](#required-environment)
  + [Shell](#shell)
  + [R](#r)
  + [Python](#python)
* [Installing Prerequisites](#installing-prerequisites)
  + [For shell package(s)](#for-shell-package-s--)
  + [For R package(s)](#for-r-package-s--)
  + [For Python package(s)](#for-python-package-s--)
* [Quick Guide](#quick-guide)
  1. [Working directory](#1-working-directory)
  2. [Download shell scripts and configuration](#2-download-shell-scripts-and-configuration)
  3. [Executing script](#3-executing-script)
* [Standard Procedure](#standard-procedure)
  1. [Working directory](#1-working-directory-1)
* [Interpreting Results](#interpreting-results)

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
***
### For shell package(s):

    Please refer to the homepage of each package for more information.

### For R package(s):

    Run `install.packages("package_needed")` in R.

### For Python package(s):

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

***

### 1. Working directory
Please start with the following folder structure.
```
work_dir/
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

- Download `main.sh` and `statistics.sh` files from the Master branch to your working directory.

- Input the required configurations in `main.sh` as follows.
    ```
    # Required Settings
    
    workdir=/Users/.../work_dir
    mapping_db=/Users/.../reference_genome.fa
    trimmomatic_dir=/Users/.../Trimmomatic-0.39/trimmomatic-0.39.jar
    adapter=${workdir}/MIGadapter.fasta
    popmap_dir=${workdir}/popmap.txt
    
    # Read file names
    
    files="d1-LE
    d2-LE"
    ```
- Optional configurations

    In case of analyzing more than 2 samples, it is recommended to set `stacks_r = 0.8` rather than `stacks_r = 1.0` by default.
> For more information, please refer to the documentation of Trimmomatic and Stacks.

### 3. Executing script
Please execute the following commands in the same directory as `main.sh` and `statistics.sh`,

```
./main.sh
./statistics.sh
```

or simply copy all the commands and execute them.

## Standard Procedure

***

### 1. Working directory
Please start with the following folder structure.

```
work_dir/
├── BaseCall
│   ├── sample1_R1_001.fastq.gz
│   ├── sample1_R2_001.fastq.gz
│   ├── sample2_R1_001.fastq.gz
│   └── sample2_R2_001.fastq.gz
├── MIGadapter.fasta → (Optional)
└── popmap.txt → (Optional)
```
Put NGS data files in `./work_dir/BaseCall` is all you need to do.

For trimming of the original data files, an adapter.fasta file is needed. By default, a `MIGadapter.fasta` file will be downloaded.

In case of multiple population analysis, you can put a `popmap.txt` file in `./work_dir/`. By default, a `popmap.txt` file will be automatically generated with all the samples recognized in the same population.

>**...still in progress, please follow the steps in Quick Guide.**

## Interpreting Results

***

