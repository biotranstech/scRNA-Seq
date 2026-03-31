scRNA-Seq RNA Editing Workflow
======

## scRNA-Seq overview
This project provides an automated workflow for single-cell RNA editing analysis, built on the Python DAGflow framework.
It enables a complete pipeline from raw FASTQ data to RNA editing results, including:
1. Single-cell RNA-seq processing
2. Expression matrix generation
3. RNA editing detection
4. Statistical analysis
5. Result output and visualization

## Why Use This Workflow
Single-cell RNA editing analysis is typically complex and fragmented:
```
| Challenge            | Traditional Approach                |
| -------------------- | ----------------------------------- |
| Multiple tools       | cellranger + custom scripts + R     |
| Complex workflow     | Multi-step manual execution         |
| Error-prone          | Path / format / parameter issues    |
| Poor reproducibility | Different users → different results |
```
✅ Advantages of This Workflow

This pipeline provides:

1. Fully automated workflow
2. Standardized processing
3. Integrated RNA editing analysis
4. Structured outputs
5. Restartable execution (checkpoint support)

## Workflow Overview
```
FASTQ
 ↓
Cell Ranger (scRNA processing)
 ↓
Expression Matrix
 ↓
RNA Editing Detection
 ↓
Statistical Analysis
 ↓
Results
```

## Project Structure
```
├── up_analysis.py        # Stage 1: scRNA processing
├── after_anno.py        # Stage 2: RNA editing analysis
├── config.py            # Configuration file (MUST EDIT)
├── common.py            # Utility modules
├── dagflow/             # DAGflow framework
├── scripts/             # Analysis scripts
├── sample.list          # Sample input file
└── human_common_marker.xlsx
```

## Installation

### Install Conda (if not installed)
Recommended:
```
Miniconda
Anaconda
```

### Create Python Environment
```
conda create -n scrna python=3.9 -y
conda activate scrna
```

### Install R Environment (Required)
Create environment.yml:
```
name: scrna_r
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base=4.2.3
  - r-seurat
  - r-tidyverse
  - r-data.table
  - r-openxlsx
  - r-optparse
  - bioconductor-annotationdbi
```
Run:
```
conda env create -f environment.yml
conda activate scrna_r
```

### Install Required Software
🔹 Single-cell processing
Install [Cell Ranger](（https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger）)
Download from official 10x Genomics website.



### What is RNA editing
RNA editing is a dynamic post-transcriptional modification with significant implications for gene regulation and disease mechanisms.

### Introduction
This is a workflow for processing single-cell RNA editing analysis, developed based on the Python package DAGflow. The workflow is divided into two main stages:
1. **up_analysis.py**: Complete the whole process from raw data processing to scRNA analysis
2. **after_anno.py**: RNA editing was performed after processing

## Workflow file directory structure

```
├── up_analysis.py             # Main workflow definition file
├── after_anno.py              # Second workflow definition file
├── common.py                  # Some workflow file processing modules
├── config.py                  # Pipeline software configuration file
├── dagflow/                   # Python DAGflow package
├── sample.list                # Sample list file
├── scripts/                   # Script directory
└── human_common_marker.xlsx   # Human marker information
```

## Installation
* REQUIREMENT
   * [Cromwell](https://github.com/broadinstitute/cromwell/releases)
   * [Cromwell](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
   * [R](https://www.r-project.org)
### Create and activate R environment
For scRNA-Seq, the R version need is over 4.2. You can create a YAML file called environment.yml with the following content:
```
name: r_environment
dependencies:
  - r-base=4.2.3
  - bioconda::bioconductor-annotationdbi
  - r-corrplot
  - r-corrplot
  - r-optparse
  - r-homologene
  - r-openxlsx
  - r-Seurat
  - r-tidyverse
  - r-data.table
  - r-plyr
  - r-ggplot
```

## Usage

### First：From raw data to scRNA analysis (up_analysis.py)

```bash
python ./up_analysis.py -i sample.xls --work_dir ./01_work --out_dir ./02_result 
```

### Second：RNA editing processing (after_anno.py)

After the first phase is complete, the second phase of processing is performed：

```bash
python ./after_anno.py --work_dir ./01_work --out_dir ./02_result -s sample.xls
```

### sample.xls
Sample information file. The first column contains sample IDs, the second and third columns contain data paths, and the fourth column contains sample grouping information: 

```json
sample1 /../demo/sample1/R22071372-J22120174-J22120174_combined_R1.fastq.gz    /../demo/sample1/R22071372-J22120174-J22120174_combined_R2.fastq.gz    CC1
sample2 /../demo/sample2/R22071373-J22120175-J22120175_combined_R1.fastq.gz    /../demo/sample2/R22071373-J22120175-J22120175_combined_R2.fastq.gz    CC1
sample5 /../demo/sample5/R22071383-J22120191-J22120191_combined_R1.fastq.gz    /../demo/sample5/R22071383-J22120191-J22120191_combined_R2.fastq.gz    CC2
sample6 /../demo/sample6/R22071384-J22120192-J22120192_combined_R1.fastq.gz    /./demo/sample6/R22071384-J22120192-J22120192_combined_R2.fastq.gz    CC2
```


## Status tracking and restart

When submitting a job using nohup, you can monitor the workflow’s progress in real time by checking the nohup.out file. When submitting the job directly in the command line, you can monitor the workflow’s progress through the messages printed on the screen. If the workflow interruption needs to be restarted, the system will automatically check the file, skip the completed task, and continue to execute from the breakpoint, improving work efficiency.

